! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module scale_declarations

  !! These declarations are separated from the bulk of !!
  !!  the others because they define the size of the   !!
  !!  soil and canopy profiles, & the time resolution. !!
  !! They also declare constructs that are needed at   !!
  !!  all levels, such as the user's config options,   !!
  !!  the time holder, and the met-drivers.            !!
  !! Physical parameters, such as pi, are also here.   !!

  implicit none


  integer,parameter :: core              = 21,  & ! number of soil layers + 1
                       fname_length      = 100, & ! length of filename variables
                       nos_soil_layers   = 20,  & ! soil layers
                       nmax              = 2,   & ! number of iterations in math-loops
                       nos_canopy_layers = 10,  & ! number of canopy layers
                       steps             = 96     ! number of timesteps per day

  ! All the time information in one place..
  type time_holder
    integer :: year = 0
    integer :: day  = 0
    integer :: step = 0
    integer :: steps_count = 0        ! count of number of steps completed so far.
    real    :: seconds_per_step = 24./steps*60.*60.  ! timesteps (3600=>fixed at 24 steps per day)
    real    :: daytime = 0.              ! date+time in days, e.g. noon on day 30 of yr 2
                                         ! (if yr=365days) == 365 + 30 + 0.5 = 395.5
  end type time_holder
  type(time_holder),save :: time


  ! meteorological inputs/drivers
  type met_drivers
    real,allocatable,dimension(:) :: ambient_co2  ! (ppm) Ambient Carbon Dioxide concentration
    real,allocatable,dimension(:) :: par          ! (units?) Photosynthetically active radiation
    real,allocatable,dimension(:) :: precip       ! (units?) precipitation
    real,allocatable,dimension(:) :: sfc_pressure ! (pa)  Atmospheric surface pressure
    real,allocatable,dimension(:) :: sw_rad       ! (units?) surface downward short-wave radiation
    real,allocatable,dimension(:) :: temp         ! (C) temperature
    real,allocatable,dimension(:) :: vpd          ! (units?) vapour pressure deficit
    real,allocatable,dimension(:) :: wind_spd     ! (ms-1) wind strength
    ! If user does not provide variables, but still wants to alter co2/pressure from defaults
    ! then they can supply these in the 'met parameters' section of the config file..
    real :: const_ambient_co2  = 360    ! (ppm) Ambient Carbon Dioxide concentration
    real :: const_sfc_pressure = 100000 ! (pa)  Atmospheric surface pressure
  end type met_drivers
  type(met_drivers),save :: met


  ! All of the configuration info for SPA...
  type user_config_holder
     ! Input files --------------------------------------------------------
     character(fname_length) :: met_filename   = '' ! Met file name/path
     logical                 :: met_file_is_nc   = .false. ! determines whether
     character(fname_length) :: soils_filename = '' ! Soils file name/path
     character(fname_length) :: veg_filename   = '' ! Veg file name/path
     character(fname_length) :: enkf_filename  = '' ! EnKF file name/path
     character(fname_length) :: initials_filename  = '' ! EnKF file name/path
     character(fname_length) :: obs_filename  = '' ! observation data (for EnKF) file name/path
     character(fname_length) :: phen_filename  = '' ! phenology (LAI seasonality) data file name/path
     ! Options ----------------------------------------------------------
     logical :: veg_is_deciduous  = .False.
     real    :: latitude          = 50.00  ! +ve == nth, -ve == sth
     real    :: longitude         = 00.00  ! 0--360. -ve not recognised.
     integer :: timesteps_per_day = steps     ! User can control whether there are more/less steps per day..
     integer :: days_in_year      = 365    ! number of (whole!) days per year
     integer :: nos_of_years      = 1      ! number of years to simulate
     logical :: co2_in_met_file   = .False.! whether or not CO2 is in the input met file
     logical :: par_in_met_file   = .False.! likewise for PAR.
     logical :: sfc_press_in_met_file = .False.! likewise for src pressure
     ! Output ----------------------------------------------------------
     character(fname_length) :: output_directory =  ''    ! directory to put output data in
     logical                 :: std_csv_output   = .true. ! Set if you want to do calcs for deciduous veg
  end type user_config_holder
  type(user_config_holder),save  :: user_opts

  real,allocatable,dimension(:) :: all_phenology

  ! General physical constants..
  real,parameter :: boltz = 5.670400e-8, & ! (W m-2 K-4)
                    pi    = 3.14159265     ! (-)

end module scale_declarations
