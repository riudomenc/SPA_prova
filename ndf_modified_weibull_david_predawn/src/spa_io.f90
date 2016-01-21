! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_io

  !! coordinates calling of all io routines, including    !!
  !! reading user-config and opening/reading input files. !!

  use scale_declarations, only: fname_length, user_opts
  !use spa_io_netcdf,      only: nc_met_file

  implicit none

  ! Not sure this is the right place to put this..
  !type(nc_met_file),pointer,save :: met_nc

  ! --- private to this module ---

  ! variables.. 
  character(fname_length)           :: spa_config_filename = "SPA.config" !default value
  integer                           :: prev_step_of_day

  private :: spa_config_filename, prev_step_of_day
  
contains
  !
  !-------------------------------------------------------------------
  !
  subroutine check_enough_timeslices( met , flag )

    ! Check that there are enough timeslices in the driver  !
    ! file to run the model for the number of days the user !
    ! wants to run for.                                     !

    use log_tools
    use scale_declarations, only : met_drivers
    use obsdrivers,         only : all_baseobs

    implicit none

    ! arguments..
    type(met_drivers),intent(in) :: met
    logical :: flag ! if TRUE, check time slices of met file, else of obs data

    ! local variables..
    integer :: nos_of_timeslices

    ! nos of available timeslices = length of time dimension.
    ! nos of required timeslices  = user_run_length in days * nos of steps per day (we read at every step)

    if ( flag ) then
      nos_of_timeslices = size( met%temp )
    else
      nos_of_timeslices = size( all_baseobs , 2 )
    endif

    if ( nos_of_timeslices .lt. user_opts%nos_of_years * user_opts%days_in_year * user_opts%timesteps_per_day ) then
      if (flag ) then
        write(message,*)"There are not enough timeslices in the met driver file to run the model "//&
                      "for the number of days declared in the config file (default is 365)"
      else
        write(message,*)"There are not enough timeslices in the obs data file to run the model "//&
                      "for the number of days declared in the config file (default is 365)"
      endif
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

  end subroutine check_enough_timeslices
  !
  !-------------------------------------------------------------------
  !
  subroutine start_spa( initialwater , A , B )

    ! Read the user config, then call the relevant routines !
    ! to open files and do any initialisation required.     !

    use config_tools,       only: ConfigSection
    use log_tools
    use scale_declarations, only: met, user_opts
    use spa_cmd_line,       only: parse_spa_cmd_line
    use spa_config,         only: read_user_config, update_parameters_from_user_config
    use spa_initialise
    use spa_io_csv
    !use spa_io_netcdf
    use enkf
    use obsdrivers,         only: baseobs

    implicit none

    ! arguments.. 
    real,intent(out)   :: initialwater
    real,dimension(ndim,nrens),intent(inout) :: A , B

    ! local variables..
    type(ConfigSection), pointer :: section_names
    real,dimension(ndim)         :: initmean ! initial mean values of ensemble state vector
    real,dimension(ndim)         :: initerr ! initial model error (or rather uncertainty)

    real :: dummy

    ! Find out the config filename..
    call parse_spa_cmd_line( spa_config_filename )

    ! Read user's config file..
    !allocate( met_nc, met_nc%header )
    call read_user_config( spa_config_filename , user_opts , section_names , dummy) !met_nc )

    ! Open the meteorology file..(but don't load any data from it yet)
    if ( user_opts%met_file_is_nc ) then

       ! load basic data associated with file...
       call write_log( "opening the NC meteorology input file" , msg_info , __FILE__ , __LINE__ )
       !call open_nc_in( met_nc%header )

       ! ! Check latitude and longitude desired by user are inside bounds of this file..
       ! call write_log( "checking it contains the desired lat-lon" , msg_info , __FILE__ , __LINE__ )
       ! call check_LatLon( met_nc%header, met_nc%grid )

       ! Load all the meteorological data into the pointer..
       !call load_met_nc_data( met_nc )
       call write_log_div

    else

       ! Load all the meteorological data into the pointer..
       call write_log( "opening the CSV meteorology input file" , msg_info , __FILE__ , __LINE__ )
       call read_met_csv( user_opts%met_filename , met )
       call write_log_div

    endif

    ! check the met-driver file has enough timeslices to perform the run..
    call write_log( "checking met file contains enough time-slices for run length" , msg_info , __FILE__ , __LINE__ )
    call check_enough_timeslices( met , .TRUE. )

    ! Read the soils file..    
    call write_log("opening the soils file" , msg_info , __FILE__ , __LINE__ )
    call read_soils_csv( user_opts%soils_filename, user_opts%veg_is_deciduous )

    ! Read the veg file..
    call write_log("opening the veg file" , msg_info , __FILE__ , __LINE__ )
    call read_veg_csv( user_opts%veg_filename , user_opts%veg_is_deciduous )

    ! Read the EnKF file..
    call write_log("opening the EnKF file" , msg_info , __FILE__ , __LINE__ )
    call read_enkf_csv( user_opts%enkf_filename , initmean , initerr )

    ! Read the observations data file for assimilation
    allocate(baseobs(maxobs))
    call write_log("opening the observations data file" , msg_info , __FILE__ , __LINE__ )
    call read_obs_csv( user_opts%obs_filename )
    call write_log( "checking obs file contains enough time-slices for run length" , msg_info , __FILE__ , __LINE__ )
    call check_enough_timeslices( met , .FALSE. )

    ! Read LAI phenology file
    call write_log("opening the phenology data file" , msg_info , __FILE__ , __LINE__ )
    call read_phen_csv( user_opts%phen_filename )

    ! Initialise the soils..
    call write_log( "Initialising the soils" , msg_info , __FILE__ , __LINE__ )
    call initialise_soils( initialwater )

    ! Initialise the ensemble A
    call initialise_A ( A , B , initmean , initerr , user_opts%initials_filename )

    ! Initialise the veg.. (On day 1 set leaf WP in each layer)
    call initialise_veg( user_opts%timesteps_per_day , A )

    ! Make sure parameters read in from the user-config overwrite any in the files/default values..
    call write_log( "Reading user-defined params from user-config" , msg_info , __FILE__ , __LINE__ )
    call update_parameters_from_user_config( section_names )

    ! Open output files..
    call write_log( "Opening output files" , msg_info , __FILE__ , __LINE__ )
    call handle_output( 0 )

  end subroutine start_spa
  !
  !-------------------------------------------------------------------
  !
  subroutine update_drivers( time )

    ! wrapper to call appropriate nc/csv routine !
    ! for reading next slice of meteorology data !

    use clim,               only: atmos_press, coa, daypar, dayppt, par_top, ppt, sw_rad, temp_bot, &
                                  temp_top, vpd_bot, vpd_top, wdbot, wdtop, wetev, wind_spd
    use hourscale,          only: discharge, runoff
    use log_tools
    use scale_declarations, only: met, time_holder
    use veg,                only: flux, gppt, respt, totass, totevap, totres, transt
    use obsdrivers,         only: all_baseobs , baseobs , sapobs
    use enkf,               only: Prades , Vallcebre

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time 

    ! local variables
    integer :: tm

    ! check we haven't been called already this step...
    if ( time%step .ne. prev_step_of_day ) then

      ! On first hour of day, set counters to zero..
      if ( time%step .eq. 1 ) then
        daypar  = 0. ; dayppt  = 0. ; discharge = 0. ; flux   = 0.
        gppt    = 0. ; respt   = 0. ; runoff    = 0. ; totass = 0.
        totevap = 0. ; totres  = 0. ; transt    = 0. ; wetev  = 0. 
      endif

      ! load this time-step's driver into each var..
      tm = time%steps_count - (( time%year -1 ) * user_opts%days_in_year * user_opts%timesteps_per_day)
      atmos_press = met%sfc_pressure( tm ) !( time%steps_count )
      coa         =  met%ambient_co2( tm ) !( time%steps_count )
      par_top     =          met%par( tm ) !( time%steps_count )
      ppt         =       met%precip( tm ) !( time%steps_count )
      sw_rad      =       met%sw_rad( tm ) !( time%steps_count )
      temp_bot    =         met%temp( tm ) !( time%steps_count )
      temp_top    = temp_bot
      vpd_bot     =          met%vpd( tm ) !( time%steps_count ) DELETE!!!
      vpd_top     = vpd_bot
      wind_spd    =     met%wind_spd( tm ) !( time%steps_count )
      ! and for obs data
      if (Prades) then
        sapobs      = all_baseobs( : , tm )
      else
        baseobs     = all_baseobs( : , tm ) !time%steps_count )
      endif

      ! Calculate the absolute water deficit (g m-3)
      wdtop = vpd_top * 217. / ( 0.1 * ( temp_top + 273.4 ) )
      wdbot = vpd_bot * 217. / ( 0.1 * ( temp_bot + 273.4 ) )

      !--- The following are only used for output files.. ---

      ! keep a (daily) running total of the precip..
      dayppt = dayppt + ppt

      ! sum the daily energy (PAR multiplied by nos of seconds per timestep)..
      daypar = daypar + par_top * time%seconds_per_step

      prev_step_of_day = time%step

      call write_log( "Loaded met and obs data for step" , msg_info , __FILE__ , __LINE__ )

    endif

  end subroutine update_drivers
  !
  !-------------------------------------------------------------------
  !
  subroutine handle_output( flag , time , iwater , output_data )

    ! wrapper to the various possible writes. !

    use clim,               only: gdd, mint
    use scale_declarations, only: time_holder
    use spa_io_csv
    !use spa_io_netcdf
    use log_tools
    use enkf

    implicit none

    integer,intent(in)      :: flag
    ! flag takes the following values...
    ! 0   = open output files
    ! 1   = perform a standard write
    ! 2-7 = specific writes for subroutines within SPA
    ! 8 = spare
    ! 9   = close output files
    type(time_holder),optional,intent(in) :: time
    real,optional,intent(inout)           :: iwater
    real,dimension(:),optional            :: output_data

    if ( ( flag .ge. 3 ) .and. ( flag .le. 6) &
         .and. ( .not. present(output_data) ) ) then
      write(message,*)"When handle_output is called with flag=",flag,&
                         " you MUST supply output_data."
      call write_log( trim(message), msg_fatal, __FILE__ , __LINE__ )
    endif

    select case (flag)
    case (0) !open output files and write headers..
      if ( user_opts%std_csv_output ) &
            call open_output_csv( user_opts%output_directory, user_opts%veg_is_deciduous )
    case (1)
      if ( user_opts%std_csv_output ) then
        if ( present( iwater ) ) then
          !call write_output_csv( time, iwater, user_opts%veg_is_deciduous )
       else
          write(message,*)"When handle_output is called with flag=1 you MUST supply iwater."
          call write_log( trim(message), msg_fatal, __FILE__ , __LINE__ )
        endif
      endif
    case (2)
      ! moved writes from predictor (allocate_carbon.f90) to here 
      call write_to_file( 31,(/gdd, mint/) )
    case (3)
      call write_assimilate_output_csv( time, output_data )
    case (4)
      call write_solar_output_csv( output_data )
    case (5)
      call write_soilday_output_csv( 1 , output_data )
    case (6)
      call write_soilday_output_csv( 2 , output_data )
    case (9)
      if (user_opts%std_csv_output)   call close_output_csv
    case default
      write(message,*)"flag supplied to handle_output:",flag," was not recognised"
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    end select

  end subroutine handle_output
  !
  !----------------------------------------------------------------------------
  !
end module spa_io
!
!----------------------------------------------------------------------------
