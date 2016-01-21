! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


!!                                                  !!
!!              ARCHAIC SYSTEM !!                   !!
!!                                                  !!
!! Used to contain all global-variable declarations !!
!! (except for scale_declarations).  I have slowly  !!
!! started to migrate these to where they belong -  !!
!! modules that contain code that uses them!        !!
!!                                                  !!


module clim

  use scale_declarations, only: nos_canopy_layers, steps

  implicit none

  ! met-drivers used/updated each timestep..
  real :: atmos_press, & ! Surface atmospheric pressure (Pa)
          coa,         & ! Ambient CO2 concentration (ppm)
          par_top,     & ! PAR at top canopy layer
          ppt,         & ! Precipitation (mm)
          sw_rad,      & ! Incident short wave radiation (Wm-2)
          temp_bot,    & ! Surface temperature at bottom and..
          temp_top,    & !  ..at top of canopy
          vpd_bot,     & ! Vapour Pressure Deficit at bottom..
          vpd_top,     & !  ..and at top of canopy.
          wind_spd       ! Wind speed (m/s)

  ! other meteorology-related variables.. 
  real    :: avtemp,     & ! Average daily temperature (Celcius); used in deciduous
             daypar,     & ! Accumulation of PAR over a day (umol.m-2.s-1)
             dayppt,     & ! Accumulation of precipitation over a day (mm.d-1)
             daytempsum, & ! ?
 gbw(nos_canopy_layers), & ! Boundary layer conductance for each canopy layer (m.s-1)
             gdd,        & ! ?
             max_fol,    & ! Deciduous switch to indicate that maxfoliage has been reached
             mint,       & ! Minimum daily temperature
             multtf,     & ! ?
             multtl,     & ! Deciduous switch to indicate labile turnover has become
             rnet,       & ! ?
             wdbot,      & ! wind speed at bottom of canopy (m.s-1)
             wdtop,      & ! wind speed at top of canopy (m.s-1)
             wetev(steps)  ! Wet evaporation from canopy (mm.t-1) per time step

  real,parameter :: par_ratio = 4.6  ! parameter, PAR energy convertion ratio with SW radiation (umol.J-1)

end module clim
!
!-------------------------------------------------------------------
!
module hourscale

  implicit none

  integer :: hour            ! ?
  real    :: canopy_store, & ! ?
                  dayevap, & ! soil moisture flux '=-Qe'
                discharge, & ! ?
               evap_store, & ! ?
                      gaw, & ! boundary layer conductance (m.s-1), calculated each time step 
                      gws, & ! soil conductance to water vapour 
                  hourppt, & ! Met -> precip
                hourpress, & ! Met -> surface pressure
                  hourrad, & ! Soil net radiation balance including both long and shortwave (W.m-2), used in long wave determination
                 hourrnet, & ! Met -> sw rad
                 hourtemp, & ! Met -> temperature (oC) converted to (K)
                 hourtime, & ! ?
                   hourts, & ! ?
                  hourvpd, & ! Met -> VPD
                 hourwind, & ! Met -> wind speed (m.s-1)
                 overflow, & ! over flow from water input into the soil layer. Set as a fixed proportion of surface_watermm (mm)
                       Qc, & ! ?
                       Qe, & ! Latent energy flux (soil) within the model structure (W.m-2)
                       Qh, & ! ?
                       Qn, & ! Net LW emissions (W.m-2) - based upon the emissivity and the total radation penetration
                   runoff, & ! ?
          surface_watermm, & ! timestep calculated surface water. i.e. the water that cannot be infiltrated within the given timestep or surface layer saturation (mm) 
                    totet, & ! ?
                underflow, & ! ?
            unintercepted    ! ?

  real,parameter :: freeze = 273.15 ! parameter, freezing point of water in Kelvin

end module hourscale
!
!-------------------------------------------------------------------
!
module Hydrol

  use scale_declarations, only : core

  implicit none

  real,dimension(core) :: claypc, & ! Percentage of soil that is clay.
                          sandpc    ! Percentage of soil that is sand.

end module Hydrol
!
!-------------------------------------------------------------------
!
module Irradiance_Sunshade

  implicit none

  real ::check
  real :: soilnet  ! Net radiation balance at the soil surface (W.m-2)
   ! (used in light, soil_funcs, and spa_io_csv)

end module Irradiance_Sunshade
!
!-------------------------------------------------------------------
!
module metab

  use scale_declarations, only : nos_canopy_layers

  implicit none

  real :: an, & ! GPP for given canopy layer in the shade or light leaf area loops (umolC.m-2.timestep-1)
          ci, & ! intra leaf CO2 concentration (ppm)
          et, & ! Evapotranspiration through stomata,for a given canopy layer. Both (g.m-2.s-1) and (mol.m-2.s-1) are used
          ht, & ! Height of specific canopy layer, used in spa_canopy.F loop (m)
 layer_capac, & ! Canopy layer specific canopy capacitance, based on canopy area
        resp, & ! Canopy layer respiration (umolC.m-2)
          rn, & ! Respiration constant for leaf area under photosynthetic analysis (umol CO2.gN.m-2 at 10oC)
      rplant, & ! Plant hydraulic resistance for each canopy layer (MPa.s-1.mmol-1)
       rsoil, & ! Root and soil resistence for a given canopoy layer (MPa.s-1.mmol-1)
         vcm, & ! Maximum rate of carboxylation (umol CO2.gN.m-2.s-1); based on kappaV coefficient
         vjm    ! Maximum rate of electrob transport (umol.m-2.s-1)
  real,dimension(nos_canopy_layers) :: rplant_canopy , rsoil_canopy

  real,parameter  :: opttemp    = 30.0   ! parameter, metabolic temperature optimum
  real,parameter  :: propresp   = 1.0    ! parameter, proportional respiration constant based on N content and assumed temperature base of 10 oC
  real,parameter  :: Vckurtosis = 0.143  ! parameter, kurtosis of Vcmax temp response
  real,parameter  :: Vjkurtosis = 0.172  ! parameter, kurtosis of Jmax temp response

end module metab
!
!-------------------------------------------------------------------
!
module meteo

  implicit none

  real :: gbb, & !
        hflux, & !
           la, & !
          nit, & !
          par, & !
         psil, & !
         psip, & !
         psis, & !
          rad, & !
         temp, & !
         wdef    !

  real,parameter :: cp   = 0.001012  ! parameter, spec. heat of air, (KJ g-1 K-1)
  real,parameter :: gi   = 1.0       ! parameter, mesophyll conductance (m s-1) set high, and thus ignored, /0.0025/ old value
                                     ! v large mesophyll conductance - ACi curves have not been Epron adjusted
  real,parameter :: head = 0.009807  ! parameter, head of pressure  (MPa/m)
  real,parameter :: Rcon = 8.3144    ! parameter, gas constant

end module meteo
!
!-------------------------------------------------------------------
!
module soil_structure

  use scale_declarations, only : core, nos_soil_layers

  implicit none

  logical :: SaxtonNew
  integer ::  pass,  & ! pass indicator for soil layer in the Saxton equations
             rootl     ! number of root layers penetrated into the soil layers
  real :: drythick,  & !
         max_depth,  & !
       max_storage,  & !
          resprate,  & !
       rootbiomass,  & ! root biomass (g Biomass.m-2) in total is determined as twice the root C (gC.m-2)
           rootrad,  & !
         rootreach,  & !
   model_rootreach,  & !
  root_distr_param,  & !
            root_k,  & !
              snow,  & !
       surfbiomass,  & !
      through_fall,  & !
      weighted_SWP,  & !
           theta33,  & !
   SWCgravelfactor,  & !
              psie,  & !
              KbKs,  & !
            thetas,  & !
     lambda_Saxton,  & !
            SWCobs,  & !
     weighted_soilR    !
  real           :: abovebelow = 1.    ! saxton water retention equation are off by default
  real           :: rootdens   = 0.5e6 ! root density (g biomass m-3 root). Set as a constant value
  real           :: rootresist = 400.  ! number of m of roots in a layer to give resistance of 1 MPa m2 s mmol-1
  real           :: thermal    = 1.58  ! conversion rate of thermal conductivity from (W m-1 K-1) to (J m-1 K-1 h-1)
  real :: watericemm(nos_soil_layers) !
  real,dimension(core) :: conduc,  & !
                          cond1,   & !
                          cond2,   & !
                          cond3,   & !
                     draincheck,   & !
                 field_capacity,   & !
                      gravelvol,   & !
                fraction_uptake,   & ! fraction of evapotranspiration (i.e. root draw) from each layer. Reset at each timestep
                        iceprop,   & !
                    layer_depth,   & !
                    mineralfrac,   & !
                    organicfrac,   & !
                       porosity,   & !
                           potA,   & !
                           potB,   & !
                        pptgain,   & !
                       rootfrac,   & !
                     rootlength,   & !
                       rootmass,   & ! root mass (g biomass) per soil layer are determined every timestep in root call.
                       soiltemp,   & !
                soiltemp_nplus1,   & !
                            SWP,   & !
                          soilR,   & !
                         soilR1,   & !
                         soilR2,   & !
                      thickness,   & !
                      waterfrac,   & !
                      waterloss,   & ! Water lost (mm) from a given soil layer; reset with each timestep
                      watergain      !
  real,dimension(10) :: wettingtop = 0. ! surface layer wetting fronts
  real,dimension(10) :: wettingbot = 0. ! surface layer wetting fronts

end module soil_structure
!
!-------------------------------------------------------------------
!
module Snow_Info

  implicit none

  real            :: delta_t,   & !
                snow_watermm,   & !
                    snowheat,   & !
                  snowheight,   & !
                  snowweight,   & !
                    soilflux,   & !
                 totsoilflux      !
  real,dimension(3) :: f_I,     & !
                   ro_snow,     & !
                     snowh,     & !
                    snowht,     & !
                     snoww,     & !
                  totwaterflux, & !
                     Tsnow,     & !
                     waterflux, & !
                     W_L          !

end module Snow_Info
!
!-------------------------------------------------------------------
!
module veg

  use scale_declarations, only : nos_canopy_layers,steps

  implicit none

  integer :: nlink        = 1  !
  integer :: conductivity = 1  ! 1=>conductivity set, 0=>conductance set
  real         :: altitude,  & !
                       avN,  & !
                     canht,  & !
                     capac,  & !
                    co2amb,  & !
                    gplant,  & !
                      iota,  & !
                    kappac,  & !
                    kappaj,  & !
                       lat,  & !
                       LMA,  & !
                    minlwp,  & !
                   modRnet,  & !
                     prevC,  & !
                     sensh,  & !
                    totass,  & !
                   totevap,  & !
                     totlai, & !
                      totn,  & !
                    totres,  & !
                    towerht, & !
                     rbelow, & !
                    root_leaf_ratio
  real            :: theta(16) !A(22),  & !
                   !theta(16)   !
  real      :: dimen    = 0.08 ! HF leaf dimension
  real,dimension(nos_canopy_layers) :: LWPstore = 0., predawn_lwp = 0. ! initial LWP=0.
  real,dimension(nos_canopy_layers) :: &
               canopy_soil_resistance,  & !
                                 CSR1,  & !
                                 CSR2,  & !
                               lafrac,  & !
                                  lai,  & !
                              layerht,  & !
                                nfrac,  & !
                                  Nla,  & !
                          stom_conduct, &
                                 PWPstore !
  real,dimension(steps)         :: ess, & ! before: real,allocatable,dimension(:)
                                  gppt, & !
                                 respt, & !
                              soiletmm, & !
                                transt    !
  real,dimension(steps,nos_canopy_layers) :: flux ! before: real,allocatable,dimension(:,:)

end module veg

!-------------------------------------------------------------------
