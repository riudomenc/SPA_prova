! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_initialise

  !! This module gives initial values to variables at the !!
  !!  start of a 'from-scratch' SPA simulation.           !!
  !! The alternative is to start SPA from some previous   !!
  !!  state, which does not require these steps.          !!

  ! procedure private to this module..
  private :: open_file

contains
  !
  !-------------------------------------------------------------------
  !
  subroutine initialise_soils( initialwater )

    ! Go through the steps necessary to initialise the soils from scratch !

    use hourscale,          only: canopy_store, hourts
    use meteo,              only: gi
    use scale_declarations, only: nos_soil_layers
    use soil_air,           only: saxton_parameters, soil_porosity, water_retention
    use soil_structure,     only: drythick, porosity, soiltemp, thickness, &
                                  waterfrac, wettingbot

    implicit none

    ! arguments..
    real,intent(out) :: initialwater

    ! local variables..
    integer          :: i

    call saxton_parameters ! calculate the conductivity parameters
    call water_retention

    canopy_store  = 0.
    drythick      = 0.1
    hourts        = soiltemp(1)  ! initial soil temp estimate for light routine
    wettingbot(1) = thickness(1) ! top layer is saturated

    call soil_porosity      ! calculate the soil porosity

    ! v large mesophyll conductance - ACi curves have not been Epron adjusted
    initialwater  = 0.
    do i = 1 , nos_soil_layers
      initialwater = initialwater + 1e3 * ( waterfrac(i) * thickness(i) )
    enddo

  end subroutine initialise_soils
  !
  !-------------------------------------------------------------------
  !
  subroutine initialise_veg( tsteps , A )

    ! Calculate the initial leaf-water potential, which requires !
    ! knowing the soil water potential, which in turn requires   !
    ! knowing about the conductance and resistance of the soil.  !

    use meteo,              only: head
    use scale_declarations, only: core, nos_canopy_layers
    use soil_air,           only: slayer, soil_conductivity, soil_resistance, &
                                  soil_water_potential, water_uptake_layer
    use soil_structure,     only: conduc, rootl, rootlength, rootmass, SWP, waterfrac, weighted_SWP
    use veg,                only: ess, flux, gppt, lafrac, lai, layerht, LMA, &
                                  LWPstore, respt, soiletmm, transt
    use enkf
    use math_tools

    implicit none

    ! arguments..
    integer,intent(in) :: tsteps
    real,dimension(ndim,nrens),intent(inout) :: A

    ! local variables..
    real    :: dummy_output
    integer :: i , memb

    !allocate( ess(tsteps)   , flux(tsteps,nos_canopy_layers) , gppt(tsteps) , &
    !          respt(tsteps) , soiletmm(tsteps)               , transt(tsteps) )

    !if (Vallcebre) call soil_water_potential() ! first calculate SWP, as soil_conductivity for Vallcebre
      ! is calculated as a function of SWP
      ! SWP is a function of waterfrac in soil_water_potential


    ! Calculate soil-conductance. This requires the saxton parameters
    ! to have been calculated already (which happens in init_soils)
    do slayer = 1 , core   ! loop over layers..
      conduc(slayer) = soil_conductivity( waterfrac(slayer) )
    enddo

    ! Need to provide some (made-up!) initial values..
    rootlength = 0.1
    rootmass   = 0.1
    rootl      = core

    do memb = 1 , nrens
      if ( A( 8 , memb ) .eq. 0. ) A( 8 , memb ) = 0.1
    enddo
    lai        = lafrac * mean( A( 8 , : ) ) / LMA

    ! Calculate resistance to flow (we need soilR)
    call soil_resistance( A , -999 )
    ! and the soil water potential..
    call soil_water_potential()
    ! ..use to calculate the weighted soil-water-potential..
    call water_uptake_layer( dummy_output )
    ! And finally used to get the initial leaf-water potential..
    do i = 1 , nos_canopy_layers     
       LWPstore(i) = weighted_SWP - head * layerht(i)
    enddo
 
  end subroutine initialise_veg
  !
  !-------------------------------------------------------------------

  subroutine initialise_A ( A , B , initmean , initerr , filename )

  ! subroutine to initialise the model ensemble (mean and spread)
  ! will be executed only once before main calculations begin

    use enkf
    use math_tools
    use config_tools,       only: ConfigSection
    use log_tools
    use scale_declarations, only: met, user_opts
    use spa_cmd_line,       only: parse_spa_cmd_line
    use spa_config,         only: read_user_config, update_parameters_from_user_config
    use spa_io_csv

    implicit none

    ! arguments
    real,dimension(ndim,nrens),intent(inout) :: A , B
    real,dimension(ndim),intent(inout)       :: initmean ! initial mean values of ensemble state vector
    real,dimension(ndim),intent(inout)       :: initerr ! initial model error (or rather uncertainty)

    ! local variables
    real                    :: check,minA,maxA
    real,save               :: modvarflag = 0.
    real,dimension(nrens)   :: R
    real,dimension(5,nrens) :: initial_parameters
    integer,dimension(4)    :: iseed
    integer                 :: i ,j
    character( len = 10)    :: dummy
    character(len=*),intent(in) :: filename


    if (modvarswitch) modvarflag = 1.
    if (read_initial_distribution) then
        call open_file( filename , 24 , readonly=.true. )
        read(24,*) dummy
        do i = 1 , nrens
          read(24,*) initial_parameters(:,i)
        enddo
    endif

    iseed = (/188,2541,3560,4017/)

    do j = 1 , ndim
      !call slarnv( 3, iseed, nrens, R )   ! random numbers
      !R = R - mean( R )
      !if ( stddev( R ) .ne. 0) R = R / stddev( R )
      do i = 1 , nrens
        !B( j , i ) = R( i ) ! B is a matrix with zero mean and unit variance in each row,
          ! used for analysing spurious correlations and estimating size of inflation parameter
        if( j .le. nsv ) then
          A( j , i ) = initmean( j ) !+ R( i ) * initerr( j ) * initmean( j ) * modvarflag    ! state variables
        else
          if ((read_initial_distribution) .and. (j .lt. 29)) initmean( j ) = initial_parameters((j-nsv),i)
          check = initmean( j ) !+ R( i ) * initerr( j ) * initmean( j ) * modvarflag ! parameters
          !A( j , i ) = initmean( j ) + R( i ) * initerr( j ) * initmean( j ) * modvarflag ! parameters
          if ( check .lt. lo( j ) ) check = lo( j )  ! parameters to be within preset limits (lo/hi)
          if ( check .gt. hi( j ) ) check = hi( j )
          !if ( A( j , i ) .lt. lo( j ) ) A( j , i ) = lo( j )  ! parameters to be within preset limits (lo/hi)
          !if ( A( j , i ) .gt. hi( j ) ) A( j , i ) = hi( j )
          A( j , i ) = asin( ( ( check - lo( j ) ) / ( hi( j ) - lo( j ) ) ) ** 0.5 )
        endif
      enddo
      !if ( j .gt. nsv ) then
        !minA = minval( A( j , : ) )
        !if ( minA .lt. lo( j ) ) then
        !  do i = 1 , nrens
            !A( j , i ) = A( j , i ) + lo( j ) - minA
        !  enddo
        !endif
        !maxA = maxval( A( j , : ) )
        !if ( maxA .gt. hi( j ) ) then
        !  do i = 1 , nrens
        !    A( j , i ) = A( j , i ) + hi( j ) - maxA
        !  enddo
        !endif
        !do i = 1 , nrens
        !  A( j , i ) = asin( sqrt( ( A( j , i ) - lo( j ) ) / ( hi( j ) - lo( j ) ) ) )
        !enddo
      !endif
    enddo ! end loop over state variables

    j = 0
    do i = nsv+1 , ndim
      j = j + 1
      minerr( j ) = abs(0.25 * initerr( i ) * initmean( i )) ! calculate minimum error standard deviation of parameters, see Aksoy06
    enddo

  end subroutine initialise_A
  !

!=========================================================================

  subroutine open_file( filename , unitnos , readonly , header )

    ! This s/r is a wrapper to open(). It performs !
    ! sanity checks when opening files.  Control   !
    ! of write-permission is available through the !
    ! optional readonly argument, and if header is !
    ! provided it will be written to the file.     !

    use log_tools

    implicit none

    ! arguments..
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: unitnos
    logical,optional, intent(in) :: readonly
    character(len=*),optional, &
                      intent(in) :: header

    ! local variables..
    integer           :: ios
    logical           :: read_mode
    character(len=20) :: file_status, write_status

    ios = 0

    ! determine whether reading or writing..
    ! (base assumption is that it's okay to write)
    read_mode    = .false.
    file_status  = 'replace'
    write_status = 'readwrite'
    if ( present( readonly ) ) then
      if ( readonly ) then
        read_mode    = .true.
        file_status  = 'old'
        write_status = 'read'
      endif
    endif

    ! open file..
    open( unit=unitnos , file=trim(filename) , iostat=ios ,    &
          status=trim(file_status) , action=trim(write_status) )

    ! check opened file okay..
    if ( ios .ne. 0 ) then
      write(message,*)"Problem opening file : ",trim(filename)
      call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! write header if provided..
    if ( present( header ) ) then
      if ( read_mode ) then
        write(message,*)"You cannot write a header to a file that "&
                      //"you have opened in read-only mode!"
        call write_log( message , msg_warning , __FILE__ , __LINE__ )
      else
        write(unitnos,*)'Time (days),',trim(header)
      endif
    endif

  end subroutine open_file
!
!----------------------------------------------------------------------------

end module spa_initialise

