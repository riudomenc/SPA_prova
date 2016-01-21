! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


program main_spa

  !! Soil Plant Atmosphere model coupled to DALEC and EnKF !!

  use allocate_carbon
  use canopy
  use log_tools
  use scale_declarations
  use spa_io
  use spa_io_csv
  use enkf
  use soil_structure , only: rootresist,SWCobs
  use soil_functions
  use veg !, only: gplant , minlwp , kappac , lafrac
  use clim
  use snow_info
  use hourscale
  use irradiance_sunshade
  use metab
  use meteo
  !use omp_lib
  use math_tools

  implicit none

  integer            :: yr, day, step, memb
  logical            :: frstrn = .true. , update = .true.
  real               :: iwater
  real               :: A(ndim,nrens), covmat(ndim,ndim), cormat(ndim,ndim)
  real               :: B(ndim,nrens), hbase(ndim,maxobs) !ensemble, covariance and correlation matrices
  real               :: SWC(core),LWPout(nos_canopy_layers),LWPout_sd(nos_canopy_layers)
  real               :: av(ndim),as(ndim),runoff_out,BGtot(core),BG1(core),BG2(core)
  real,dimension(nos_canopy_layers) :: CSR,rsoilout,rplantout
  character(len=200) :: header

  ! open SWC obs file
  OPEN(UNIT=2810,FILE='input/Prades_SWC_gapfilled.csv',status='old')
  !read(2810,*) header

  ! start logging..
  call open_log( unit=100 , fname="spa.log" )

  ! read user config, open files & initialise (spa_io.f90)
  call start_spa( iwater , A , B )
  do memb = 1 , nrens
    call update_ensemble_variables( memb , .false. , frstrn , SWC , LWPout , &
      LWPout_sd, update, runoff_out, BGtot, BG1, BG2 , CSR, rsoilout , rplantout )
  enddo

  ! set up observation operator for EnKF
  call obsop( hbase )

  ! for each year..
  do yr = 1 , user_opts%nos_of_years

     call increment_time( 'year' , time )

     ! for each day...
     do day = 1 , user_opts%days_in_year
        write(*,'(i4)',advance='no') day

        call increment_time( 'day' , time)

        totlai = all_phenology( day )

        ! for each sub-daily slice...
        do step = 1 , user_opts%timesteps_per_day

           read(2810,*) SWCobs
           !write(*,*) SWCobs

           call increment_time( 'step' , time )

           ! get met driver and obs data for step (spa_io.f90)
           call write_log( 'Get next chunk of met and obs data' , msg_info , __FILE__ , __LINE__ )
           call update_drivers( time )

           ! transfer carbon and water (canopy.f90)
           call write_log('Entering the main calculations (roots and hour)')

           ! for each ensemble member
           do memb = 1 , nrens
             call update_ensemble_variables( memb , .true. , frstrn , SWC , &
               LWPout , LWPout_sd , update , runoff_out, BGtot, BG1, BG2 , CSR , rsoilout , rplantout )
             call transform_parameters( A , .false. , time , memb )
             iota       = A( 24 , memb )
             gplant     = A( 25 , memb )
             capac      = A( 26 , memb )
             root_leaf_ratio  = A( 27 , memb )
             call roots( A , memb )
             call timestep_calcs( A , B , memb , time )
             call transform_parameters( A , .true. , time , memb )
             call update_ensemble_variables( memb , .false. , frstrn , SWC , &
               LWPout , LWPout_sd , update , runoff_out, BGtot, BG1, BG2 , CSR , rsoilout , rplantout )
           enddo ! ends loop over ensemble members

           frstrn = .false.
           do memb = 1 , ndim
             av(memb) = mean(A(memb,:))
             as(memb) = stddev(A(memb,:))
           enddo

           ! execute enkf routines (e.g. analysis + smoother) after timestep_calcs are finished
           !call enkf_calcs ( A , B , time , hbase )

           ! write output if needed (spa_io.f90)
           call write_log('Dealing with any output (if needed)' , msg_info , __FILE__ , __LINE__ )
           !call handle_output( 1 , time , iwater=iwater )
           !call write_output_csv( time, iwater, user_opts%veg_is_deciduous , A )
           call write_enkf_output ( A , B , SWC , LWPout , LWPout_sd , runoff_out , &
             BGtot, BG1, BG2 , CSR , rsoilout , rplantout )
           call write_log_div

        enddo ! ends loop over sub-daily slices

        call write_log('              ' , msg_info , __FILE__ , __LINE__ )
        call write_log_div

     enddo ! ends loop over days
  enddo ! ends loop over years

  if ( smooth .eq. 1 ) then
    !call smoothout()
  endif

  ! finish logging..
  call close_log( 0 )

contains
  !
  !-------------------------------------------
  !
  subroutine increment_time( type , time )

    ! Method to update the time-holder variable. !
  
    use scale_declarations,only: time_holder, user_opts
    use log_tools,         only: write_log

    implicit none

    ! arguments..
    !integer,intent(in)              :: value
    character(*),intent(in)         :: type
    type(time_holder),intent(inout) :: time

    select case (type)
    case ('year')
      ! new year, reset day & step..
      time%year = time%year + 1
      time%day  = 0
      time%step = 0
      write(message,*)'year is: ',time%year
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('day')
      ! new day, reset step...
      time%day  = time%day + 1
      time%step = 0
      write(message,*)'day is: ',time%day
      call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case ('step')
      ! next step
      ! Also update the count of total-steps, and the daytime..
      time%step = time%step + 1
      write(message,*)'step-of-day is: ',time%step
      call write_log( trim(message), msg_info , __FILE__ , __LINE__ )
      time%steps_count = time%steps_count + 1
      write(message,*)'number of steps so far is: ',time%steps_count
      call write_log( trim(message), msg_info , __FILE__ , __LINE__ )

    ! Re-calculate the current time in units of days..
    ! Consider step 1 to be at 00:00 on the day, and the
    ! last step of the day to be at or before 23.59...
    time%daytime = ( time%year -1 ) * user_opts%days_in_year    &
                    + time%day                                  &
                     + ( real(time%step) - 1 ) / real(user_opts%timesteps_per_day)
    write(message,*)'daytime is: ',time%daytime
    call write_log( trim(message) , msg_info , __FILE__ , __LINE__ )
    case default
      write(message,*)'Do not recognise type of update: ',type
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    endselect

  end subroutine increment_time
  !
  !-------------------------------------------
  !

end program main_spa
