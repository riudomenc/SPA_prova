! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

module spa_io_netcdf

  !! This module declares derived types describing SPA-specific  !!
  !! input/output files, and procedures for reading those files. !!

  use scale_declarations, only: fname_length
  use netcdf_tools

  implicit none

  ! Required for finding particular spot on grid that we will use..
  type nc_grid
     integer,dimension(2) :: indices     = (1,1)   ! bottom-left corner of box of grid-pts
     ! that encompass user-desired location.
     logical           :: no_bilinear = .false. ! bilinear flag (assume false=>do calcs)
  end type nc_grid

  ! Input Meteorolgy file structure..
  type nc_met_file
     type(nc_header),pointer   :: header => NULL() ! information about the file, such as name & dims
     type(nc_variable),pointer :: sat    => NULL() ! surface air temperature
     logical                   :: sat_in_kelvin = .False. ! assume SAT in Celcius
     type(nc_variable),pointer :: coa    => NULL() ! carbon dioxide atmospheric concentration
     type(nc_variable),pointer :: par    => NULL() ! photosyntheticaly active radiation
     type(nc_variable),pointer :: ppt    => NULL() ! precipitation
     logical                   :: ppt_is_rate = .False. ! assume ppt is volume (mm) not rate (mm/sec)
     type(nc_variable),pointer :: swrad  => NULL() ! shortwave radiation
     type(nc_variable),pointer :: sfc_p  => NULL() ! surface atmospheric pressure
     type(nc_variable),pointer :: vpd    => NULL() ! vapour pressure deficit
     type(nc_variable),pointer :: windsp => NULL() ! wind speed
     type(nc_grid),pointer     :: grid   => NULL() ! point on grid that we will use
  end type nc_met_file

  ! --- private to this module ---
  private :: is_value_within_dim

contains
  !
  !-------------------------------------------------------------------
  ! reading routines...
  !-------------------------------------------------------------------
  !
  subroutine open_nc_in( header )

    ! Open a netcdf input file, and fill out the header with !
    ! information on handles to the file, and the dimensions !
    ! of time, latitude and longitude.                       !

    use netcdf_tools

    implicit none

    ! arguments..
    type(nc_header),pointer :: header

    ! open netCDF file
    call open_nc_file( header )

    ! populate the time,lat & lon dims with their dim/var handles,
    !  dim length and the actual dim values...
    call get_dim_info( header%id, header%time )
    call get_dim_info( header%id, header%lat  )
    call get_dim_info( header%id, header%lon  )

  end subroutine open_nc_in
  !
  !-------------------------------------------------------------------
  !
  subroutine check_LatLon( header , grid )

    ! Check the lat/lon defined in the user-config falls within !
    ! that of the input file. If they are, retrieve the indices !
    ! of the points that encompass the user's desired location. !

    use log_tools
    use netcdf_tools
    use scale_declarations, only: user_opts

    implicit none

    ! arguments..
    type(nc_header),pointer :: header
    type(nc_grid),pointer   :: grid

    ! local variables..

    logical :: status 
    integer :: tmp(1),grid_index
    real    :: value,spacing

    status = is_value_within_dim( header%lat, user_opts%latitude )
    if ( .not. status ) then
      write(message,*) "Requested latitude lies outside of bounds of input file latitudes!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    status = is_value_within_dim( header%lon, user_opts%longitude )
    if ( .not. status ) then
      write(message,*)"Requested longitude lies outside of bounds of input file longitudes!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! If we are still going, the required lat/lon must lie within the bounds of the dimensions, so find it..

    ! Check if the user has picked a point on the longitude grid...
    ! Find the point closest..
    tmp = minloc( header%lon%values - user_opts%longitude ) 
    grid_index = tmp(1)
    value = header%lon%values( grid_index )
    spacing = header%lon%values(1) - header%lon%values(0) ! assume constant grid spacing
    if ( value .eq. user_opts%longitude ) then
       ! set this as our index..
       grid%indices(1) = grid_index
       grid%no_bilinear = .true.
    else
       ! if not, find the nearest point just less than it...
       tmp = minloc( abs( header%lon%values - user_opts%longitude + 0.5*spacing ) )
       grid%indices(1) = tmp(1)
    endif

    ! Check if the user has picked a point on the latitude grid...
    tmp = minloc( header%lat%values - user_opts%latitude )
    grid_index = tmp(1)
    value = header%lat%values( grid_index )
    spacing = header%lat%values(1) - header%lat%values(0) ! assume constant grid spacing
    if ( ( value - user_opts%latitude ) .eq. user_opts%latitude ) then
       grid%indices(2) = grid_index
    else
       ! we only avoid bi-linear if both are true...
       grid%no_bilinear = .false.
       tmp = minloc( abs( header%lat%values - user_opts%latitude + 0.5*spacing ) )
       grid%indices(2) = tmp(1)
    endif

  end subroutine check_LatLon
  !
  !-------------------------------------------------------------------
  !
  subroutine load_met_nc_data( ncfile )

    ! Read the met input file and load data to pointers !
    ! Unlike the csv file, which is read line-by-line as!
    ! it is needed, the netcdf file contents are read   !
    ! all in one go.                                    !

    use log_tools
    use scale_declarations, only: met, time, user_opts

    implicit none

    ! arguments..
    type(nc_met_file)  :: ncfile

    ! Carbon dioxide..
    allocate( met%ambient_co2(ncfile%header%time%length) )
    if ( user_opts%co2_in_met_file ) then      
      call write_log("Try to load Ambient CO2 concentration from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%coa )
      met%ambient_co2 = reshape( ncfile%coa%data , (/ ncfile%header%time%length /) )
      call write_log("..successful.")
    else
      ! use a default value..
      met%ambient_co2 = met%const_ambient_co2
      call write_log("Ambient CO2 concentration set to default value (330ppm)")
    endif

    ! Precipitation
    allocate( met%precip(ncfile%header%time%length) )
    call write_log("Try to load Precipitation from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%ppt )
    call write_log ("..successful.")
    met%precip = reshape( ncfile%ppt%data , (/ ncfile%header%time%length /) )
    if ( ncfile%ppt_is_rate ) then
      call write_log("Converting precipitation from rate to volume by "&
                   //"multiplying by nos of seconds per timestep.")
      met%precip = met%precip * time%seconds_per_step
    endif

    ! Short wave radiation
    allocate( met%sw_rad(ncfile%header%time%length) )
    call write_log("Try to load SW radiation from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%swrad )
    call write_log("..successful.")
    met%sw_rad = reshape( ncfile%swrad%data  , (/ ncfile%header%time%length /) )

    ! Surface air temperature
    allocate( met%temp(ncfile%header%time%length) )
    call write_log("Try to load Surface air temperature from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%sat )
    ! Adjust to Kelvin, if needed, and also check that units seem sane..
    if ( ncfile%sat_in_kelvin ) then
      if ( minval( ncfile%sat%data ) .lt. 100. ) then
        write(message,*)"Units of temperature appear to be wrong. You "//&
             "specified Kelvin, but minimum is below 100K; is this correct??"
        call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
      endif
      ncfile%sat%data = ncfile%sat%data - 273.15 ! convert T to Celcius
    else
      if ( maxval( ncfile%sat%data ) .gt. 150. ) then
        write(message,*)"Units of temperature appear to be wrong. You "//&
             "specified Celcius, but maximum is above 150C; is this correct??"
        call write_log( trim(message) , msg_fatal , __FILE__ , __LINE__ )
      endif
    endif
    call write_log("..successful.")
    met%temp = reshape( ncfile%sat%data , (/ ncfile%header%time%length /) )

    ! Surface pressure
    allocate( met%sfc_pressure(ncfile%header%time%length) )
    if ( user_opts%sfc_press_in_met_file ) then
      call write_log("Try to load Surface pressure from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%sfc_p )
      call write_log("..successful.")
      met%sfc_pressure = reshape( ncfile%sfc_p%data , (/ ncfile%header%time%length /) )
    else
      ! use a default value..
      met%sfc_pressure = met%const_sfc_pressure
      call write_log("Surface pressure set to default value (990mb)")
    endif

    ! Vapour pressure deficit
    allocate( met%vpd(ncfile%header%time%length) )
    call write_log("Try to load vapour pressure deficit from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%vpd )
    call write_log("..successful.")
    met%vpd = reshape( ncfile%vpd%data , (/ ncfile%header%time%length /) )

    ! Wind speed
    allocate( met%wind_spd(ncfile%header%time%length) )
    call write_log("Try to load wind speed from met netcdf file..")
    call get_nc_var( ncfile%header%id, ncfile%windsp )
    call write_log("..successful.")
    met%wind_spd = reshape( ncfile%windsp%data , (/ ncfile%header%time%length /) )
    ! adjust zero windspeeds to ensure always some (small) turbulence..
    where ( met%wind_spd .lt. 0.2 )  met%wind_spd = 0.2

    ! Photosynthetically active radiation
    allocate( met%par(ncfile%header%time%length) )
    if ( user_opts%par_in_met_file ) then
      call write_log("Try to load photosynthetically active radiation from met netcdf file..")
      call get_nc_var( ncfile%header%id, ncfile%par )
      call write_log("..successful.")
      met%par = reshape( ncfile%par%data , (/ ncfile%header%time%length /) )
    else
      ! PAR isn't available in most climate models, so
      ! instead we use a loose approximation...
      met%par = 2.3 * met%sw_rad
    endif

  end subroutine load_met_nc_data
  !
  !-------------------------------------------------------------------
  !
  ! PRIVATE PROCEDURES...
  !
  !-------------------------------------------------------------------
  !
  logical function is_value_within_dim( dim , value )

    ! check whether a value is within the bounds of a dimension !

    use log_tools

    ! arguments..
    type(nc_dimension) :: dim
    real,   intent(in) :: value

    ! check if value lies within it...
    if ( ( minval(dim%values) .lt. value ) .and. ( maxval(dim%values) .gt. value ) ) then
       is_value_within_dim = .true.
    else 
       is_value_within_dim = .false.
    endif

  end function is_value_within_dim
  !
  !-------------------------------------------------------------------
  !
end module spa_io_netcdf
