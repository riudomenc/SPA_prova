! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_config

  !! This module reads the user-provided  !!
  !! runtime-configuration file for SPA.  !!

  use scale_declarations, only: fname_length
  use spa_io_netcdf ! for the nc_met/nc_soils/nc_veg type definitions

  implicit none

  private
  public :: read_user_config, update_parameters_from_user_config

contains
  !
  !----------------------------------------------------------------------------
  !
  subroutine read_user_config( config_filename , config_options , section_names , met_nc )

    ! reads the user config, looking for the filenames !
    ! output dir and SPA configuration flags.          !

    use config_tools
    use netcdf_tools
    use scale_declarations, only: fname_length, met, user_config_holder, time
    use log_tools

    ! arguments..
    character(fname_length),intent(in) :: config_filename !(input)  where to look for the config file
    type(user_config_holder)           :: config_options  !(output) holds the user's config choices.
    type(ConfigSection), pointer       :: section_names   !(output) structure holding sections of the configuration file   
    type(nc_met_file),pointer          :: met_nc          !(output)

    ! local variables..
    type(ConfigSection), pointer :: section
    type(nc_output), pointer     :: output
    logical                      :: exists, file_nc, ios
    integer                      :: i
    character(len=100)           :: outputdir1, outputdir2

    ! Get the section names..
    call ConfigRead( config_filename, section_names )

    ! print the user's config to the logfile
    call write_log( "The contents of the configuration file are.." , msg_info , __FILE__ , __LINE__ )
    call PrintConfig( section_names, get_logunit() )
    call write_log_div

    ! Get the General options...
    call GetSection(section_names,section,'Options')
    if ( associated(section) ) then
      call GetValue(section,'veg_is_deciduous',config_options%veg_is_deciduous )
      if (config_options%veg_is_deciduous) then
        call write_log( "Deciduous simulation." , msg_info , __FILE__ , __LINE__ )
      else
        call write_log( "Evergreen simulation." , msg_info , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'latitude_deg_north',config_options%latitude         )
      call GetValue(section,'longitude_deg_east',config_options%longitude        )
      call GetValue(section,'days_per_year',     config_options%days_in_year     )
      call GetValue(section,'number_of_years',   config_options%nos_of_years     )
    else
      write(message,*)"Could not find [Options] section in the config file.  Cannot continue without it!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Get the initialisation input files...
    call GetSection(section_names,section,'Initialisation')
    if ( associated(section) ) then
      ! soils
      call GetValue(section,'soils',config_options%soils_filename)
      if (config_options%soils_filename=='') then
        write(message,*)'Soils filename must be supplied in the user-config file.'
        call write_log( message , msg_fatal, __FILE__ , __LINE__ )
      endif
      ! vegetation
      call GetValue(section,'vegetation',config_options%veg_filename)
      if (config_options%veg_filename=='') then
        write(message,*)'Vegetation filename must be supplied in the user-config file.'
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      endif
      ! enkf
      call GetValue(section,'enkf',config_options%enkf_filename)
      if (config_options%enkf_filename=='') then
        write(message,*)'EnKF filename must be supplied in the user-config file.'
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      endif
      ! enkf initials
      call GetValue(section,'initials',config_options%initials_filename)
      if (config_options%initials_filename=='') then
        write(message,*)'EnKF initials filename must be supplied in the user-config file.'
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      endif
      ! observation data for DA
      call GetValue(section,'obs',config_options%obs_filename)
      if (config_options%obs_filename=='') then
        write(message,*)'Observation data (for EnKF) filename must be supplied in the user-config file.'
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      endif
      ! phenology data
      call GetValue(section,'phenology',config_options%phen_filename)
      if (config_options%phen_filename=='') then
        write(message,*)'Phenology filename must be supplied in the user-config file.'
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      endif
    else
      write(message,*)"Could not find [Initialisation] section in the config file."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Meteorology driver
    call GetSection( section_names , section , 'Meteorology driver' )
    if ( associated(section) ) then
      call GetValue( section , 'meteorology' , config_options%met_filename )
      if (config_options%met_filename=='') then
        write(message,*) "Meteorology filename must be supplied in the user-config file."
        call write_log( message , msg_fatal , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'file_is_netcdf',config_options%met_file_is_nc)
      if ( config_options%met_file_is_nc ) then
        allocate( met_nc%header )
        met_nc%header%name = config_options%met_filename
        call get_met_nc_varnames_from_config( section , met_nc )
      end if
      call GetValue( section , 'file_has_co2' ,  config_options%co2_in_met_file  )
      if ( .not. config_options%co2_in_met_file ) then 
        call write_log( "Using a constant CO2 concentration." , msg_info , __FILE__ , __LINE__ )
        call GetValue( section , 'constant_CO2_concentration' , met%const_ambient_co2 )
      endif
      call GetValue( section , 'file_has_par' ,  config_options%par_in_met_file  )
      if ( .not. config_options%par_in_met_file ) then
        call write_log( "PAR will be calculated from SW-radiation." , msg_info , __FILE__ , __LINE__ )
      endif
      call GetValue(section,'file_has_sfc_pressure',  config_options%sfc_press_in_met_file  )
    else
      write(message,*)"Could not find [Meteorology driver] section in the config file."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Get the general output detail..
    call GetSection(section_names,section,'Output')
    if (associated(section)) then
      call GetValue(section,'directory',config_options%output_directory)
      call GetValue(section,'standard_csv_output',config_options%std_csv_output)
    end if

  end subroutine read_user_config
  !
  !----------------------------------------------------------------------------
  !
  ! a couple of tools for handling the output linked-list..
  !
  !----------------------------------------------------------------------------
  !
  function process_output(section, output)

    use config_tools
    use log_tools
    use netcdf_tools

    implicit none

    ! arguments..
    type(ConfigSection), pointer :: section
    type(nc_output), pointer    :: output

    ! local variables..
    type(nc_output), pointer :: process_output

    process_output=>add_output(output)

    ! get filename
    allocate(process_output%header)
    call GetValue(section,'filename',process_output%header%name)
    call GetValue(section,'frequency',process_output%frequency)

    if (process_output%header%name(1:1).eq.' ') then
      write(message,*)'Error, no filename specified [netCDF output]'
       call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    end if  

  end function process_output
  !
  !----------------------------------------------------------------------------
  !
  function add_output(oc)

    ! add new element to linked list !

    use netcdf_tools

    implicit none

    ! arguments..
    type(nc_output), pointer :: add_output ! function
    type(nc_output), pointer :: oc

    allocate(add_output)

    if (associated(oc)) then
       add_output%previous => oc
       if (associated(oc%next)) then
          add_output%next => oc%next
          oc%next%previous => add_output
       end if
       oc%next => add_output
    end if

  end function add_output
  !
  !----------------------------------------------------------------------------
  !
  !
  !----------------------------------------------------------------------------
  !
  subroutine get_met_nc_varnames_from_config( section , met_nc )

    ! If the input met file is netcdf we need to check the !
    !  [Meteorology driver] section for details of what    !
    !  the variables are called in the met input file.     !

    use config_tools
    use log_tools
    use scale_declarations, only: met, user_opts
    use spa_io_netcdf,      only: nc_met_file

    implicit none

    ! arguments..
    type(ConfigSection),pointer :: section
    type(nc_met_file),pointer   :: met_nc

    ! define default names for dimensions, and then go see if the user
    ! has given different ones...
    ! (I'm sure there should be a better way of doing this...)
    allocate( met_nc%header%time   ) ; met_nc%header%time%name = 'time'
    allocate( met_nc%header%lat    ) ; met_nc%header%lat%name  = 'lat'
    allocate( met_nc%header%lon    ) ; met_nc%header%lon%name  = 'lon'
    call GetValue( section, 'time',      met_nc%header%time%name   )
    call GetValue( section, 'latitude',  met_nc%header%lat%name    )
    call GetValue( section, 'longitude', met_nc%header%lon%name    )

    ! define default names for variables.. nd then go see if the user
    ! has given different ones...
    ! (I'm sure there should be a better way of doing this...)
    allocate( met_nc%coa    ) ; met_nc%coa%name    = 'co2'
    allocate( met_nc%par    ) ; met_nc%par%name    = 'par'
    allocate( met_nc%ppt    ) ; met_nc%ppt%name    = 'ppt'
    allocate( met_nc%sat    ) ; met_nc%sat%name    = 'temp'
    allocate( met_nc%sfc_p  ) ; met_nc%sfc_p%name  = 'pressure'
    allocate( met_nc%swrad  ) ; met_nc%swrad%name  = 'sw'
    allocate( met_nc%vpd    ) ; met_nc%vpd%name    = 'vpd'
    allocate( met_nc%windsp ) ; met_nc%windsp%name = 'windsp'
    call GetValue( section, 'air_temperature',         met_nc%sat%name      )
    call GetValue( section, 'temp_in_kelvin',          met_nc%sat_in_kelvin )
    if ( user_opts%co2_in_met_file ) then
      call GetValue( section, 'carbon_dioxide',        met_nc%coa%name      )
    else
      met_nc%coa%name = ''
    endif
    if ( user_opts%par_in_met_file ) then
      call GetValue( section, 'photo_active_rad',       met_nc%par%name     )
    else
      met_nc%par%name = ''
    endif
    call GetValue( section, 'precipitation',           met_nc%ppt%name      )
    call GetValue( section, 'precip_is_rate',          met_nc%ppt_is_rate   )
    if ( user_opts%sfc_press_in_met_file ) then
      call GetValue( section, 'surface_pressure',      met_nc%sfc_p%name    )
    else
      met_nc%sfc_p%name = ''
      call GetValue(section,'constant_sfc_pressure',   met%const_sfc_pressure )
    endif
    call GetValue( section, 'sw_radiation',            met_nc%swrad%name    )
    call GetValue( section, 'vapour_pressure_deficit', met_nc%vpd%name      )
    call GetValue( section, 'wind_speed',              met_nc%windsp%name   )

    call write_log("The variable names SPA is expecting to find in met nc file are:", &
                     msg_info, __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%header%time%name),' ',trim(met_nc%header%lat%name),' ',trim(met_nc%header%lon%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%coa%name),' ',trim(met_nc%par%name),' ',  &
                     trim(met_nc%ppt%name),' ',trim(met_nc%sat%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    write(message,*) trim(met_nc%sfc_p%name),' ',trim(met_nc%swrad%name),' ', &
                     trim(met_nc%vpd%name),' ',trim(met_nc%windsp%name)
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )

  end subroutine get_met_nc_varnames_from_config
  !
  !----------------------------------------------------------------------------
  !
  subroutine update_parameters_from_user_config( section_names )

    ! make sure any parameters provided by the user overwrite !
    ! the values taken from the input files.                  !

    use config_tools
    use snow_info,          only: snowheight,snowweight
    use soil_structure,     only: rootrad
    use scale_declarations, only: met, user_opts

    implicit none

    ! arguments..
    type(ConfigSection), pointer :: section_names ! holds structure of the configuration file   

    ! local variables..
    type(ConfigSection), pointer :: section

    ! For each section, we look for possible parameters.
    ! GetValue only changes the third argument IF it
    ! finds something.

    call GetSection( section_names , section , 'Met Parameters' )
    if (associated(section)) then
    endif

    call GetSection(section_names,section,'Soil Parameters')
    if (associated(section)) then
       call GetValue(section,'snowweight',snowweight)
       call GetValue(section,'snowheight',snowheight)
    endif

    call GetSection(section_names,section,'Veg Parameters')
    if (associated(section)) then
       call GetValue(section,'rootrad',rootrad)
    endif

  end subroutine update_parameters_from_user_config
  !
  !-------------------------------------------------------------------
  !
end module spa_config
