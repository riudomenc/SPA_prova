! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module soil_air

  !! > this module needs a summary < !!

  implicit none

  integer :: slayer        ! Soil layer variable used when passing from a subroutine
                           !  to a function which is acting on a specific soil layer.
  real    :: drainlayer, & ! The field capacity of the specific layer being worked
                           !  on in the s/r soil_balance
                 liquid, & ! Liquid fraction of the water present within a given
                           !  soil layer, used in s/r soil_balance
                soilpor, & ! Soil porosity of a given soil layer, used in s/r 
                           ! soil_balance
                  unsat    ! Unsaturated volume of the soil layer below the current
                           ! one being modelled (m3 m-2). Calculated each timestep.

  SAVE

contains
  !
  !-------------------------------------------
  !
  subroutine boundary()

    ! determine boundary layer conductances at each canopy layer !

    use clim,               only: atmos_press, gbw, temp_bot, temp_top, wind_spd
    use scale_declarations, only: nos_canopy_layers
    use veg,                only: dimen, layerht, towerht

    implicit none

    ! local variables..
    integer :: i
    real    :: alpha, Dwv, mult, roughl, store, tempx, thick, u(nos_canopy_layers), xx

    alpha  = 4.
    roughl = 0.075 * towerht  ! tower height is canht+4
    do i = 1 , nos_canopy_layers
      gbw(i) = 0.
    enddo

    xx    = 0.
    store = 1.0
    do i = 1 , nos_canopy_layers
      xx    = xx + 1.
      tempx = temp_top - 0.1111 * ( xx - 1 ) * ( temp_top - temp_bot )
      Dwv   = .0000242 * ( ( ( tempx + 273.15 ) / 293.15 )**1.75 ) * 101300. / atmos_press    !Jones p 51
      ! windspeed in layer - no decay in windspeed in this open stand
      !             u(i,t)=windsp(t)
      mult   = exp( alpha * ( layerht( i ) / towerht - 1. ) )
      u( i ) = wind_spd * mult
      thick  = 0.004 * ( dimen / u( i ) )**0.5    ! boundary layer thickness
      ! conductance to water vapour (m s-1 - not mol m-2 s-1 as in Amthor p.333)
      ! i.e. remove P/RT
      gbw( i ) = Dwv / thick
    enddo

  end subroutine boundary
  !
  !-------------------------------------------
  !
  subroutine saxton_parameters()

    ! Key Saxton equations !

    use hydrol,             only: claypc, sandpc
    use scale_declarations, only: core
    use soil_structure,     only: cond1, cond2, cond3, potA, potB, porosity, SWCgravelfactor, SaxtonNew, theta33, psie, thetas, KbKs, gravelvol, lambda_Saxton
    use enkf,               only: Prades, Vallcebre

    implicit none

    integer :: i
    real    :: A, B, CC, D, E, F, G, H, J, K, mult1, mult2, mult3, P, Q, R, T, U, V
    real    :: organicpc, theta1500t, theta1500, theta33t, thetas33t, thetas33, density, S, C ! new Saxton parameters
    real    :: alpha, lambda, Ks, Ktheta, densityfactor=1., thetasdf, gravelpc=0., thetasold ! new Saxton parameters
    real    :: psiet, bulkdensity ! new Saxton parameters
    logical :: densityswitch

    SaxtonNew = .true. ! use new Saxton equations to calculate SWP and SHC !bulky
    densityfactor = 1. ! density adjustment factor (0.9-1.3), <1 if density larger than estimated (e.g. compacted agricultural soil), and vice versa
    densityswitch = .false.    ! if densityfactor NE 1, then porosity and sand adjusted saturation have to be recalculated
    gravelpc = 0. ! volumetric decimal percentage of gravel
    organicpc = 0. ! volumetric percentage of organic material
    if (Vallcebre) then
      gravelpc = 0.19
      organicpc = 0.5 ! absolute percentage (0-100)!!! 2.08 is average in Saxton soils, but Josep says its prob. less than 1% for VC and PR
      densityfactor = 1.
      densityswitch = .false.
    endif
    if (Prades) then
      gravelpc = 0.46 ! derived from Prades total and fine bulk density data
      organicpc = 0.5 ! absolute percentage (0-100)!!! 2.08 is average in Saxton soils
      densityfactor = 1.
      densityswitch = .false.
    endif
    mult1 = 100.
    mult2 = 2.778e-6
    mult3 = 1000.0
    A = -4.396    ;  B = -0.0715   ; CC = -4.880e-4 ; D = -4.285e-5
    E = -3.140    ;  F = -2.22e-3  ;  G = -3.484e-5 ; H =  0.332
    J = -7.251e-4 ;  K =  0.1276   ;  P = 12.012    ; Q = -7.55e-2
    R = -3.895    ;  T =  3.671e-2 ;  U = -0.1103   ; V =  8.7546e-4

    do i = 1 , core

       ! old Saxton
       potA(i)  = exp( A + B * claypc(i) + CC * sandpc(i) * sandpc(i) + D * sandpc(i) * sandpc(i) * claypc(i) ) * 100.
       potB(i)  = E + F * claypc(i) * claypc(i) + G * sandpc(i) * sandpc(i) * claypc(i)
       cond1(i) = mult2
       cond2(i) = P + Q * sandpc(i)
       cond3(i) = R + T * sandpc(i) + U * claypc(i) + V * claypc(i) * claypc(i)

       ! new saxton
       if (SaxtonNew) then
         S = sandpc(i)/100.; C = claypc(i)/100. ! here, soil texture has to be in decimal percentage
         theta1500t = -0.024 * S + 0.487 * C + 0.006 * organicpc + &
            0.005 * ( S * organicpc ) - 0.013 * ( C * organicpc ) + &
            0.068 * ( S * C ) + 0.031
         theta1500 =  theta1500t + ( 0.14 * theta1500t - 0.02 ) ! 1500 kPa moisture, %v
         theta33t = -0.251 * S + 0.195 * C + 0.011 * organicpc + &
            0.006 * ( S * organicpc ) - 0.027 * ( C * organicpc ) + &
            0.452 * ( S * C ) + 0.299
         theta33 = theta33t + ( 1.283 * theta33t * theta33t - & ! 33 kPa moisture, %v
            0.374 * theta33t - 0.015 )
         potB(i) = -( log( 1500. ) - log( 33. ) ) / & ! moisture-tension coefficient
            ( log( theta33 ) - log( theta1500 ) )
         potA(i) = exp( log( 33. ) + (-potB(i) * log( theta33 ) ) ) ! moisture-tension coefficient
         thetas33t = 0.278 * S + 0.034 * C + 0.022 * organicpc - &
            0.018 * ( S * organicpc ) - 0.027 * (C * organicpc ) - &
            0.584 * ( S * C ) + 0.078
         thetas33 = thetas33t + ( 0.636 * thetas33t - 0.107 ) ! saturation 33 kPa moisture, %v
         psiet = -21.67 * S - 27.93 * C - 81.97 * thetas33 + &
            71.12 * S * thetas33 + 8.29 * C * thetas33 + 14.05 * S * C + 27.16
         psie = psiet + ( 0.02 * psiet**2 - 0.113 * psiet - 0.7)
         thetas = theta33 + thetas33 - 0.097 * S + 0.043 ! sand adjusted saturation, i.e. porosity
         density = ( 1. - thetas ) * 2.65
         if (densityswitch) then ! density adjusted by densityfactor (0.9-1.3) to account for local variations of soil density by structure or management
           density = density * densityfactor
           thetasold = thetas
           thetas = 1. - ( density / 2.65 )
           theta33 = theta33 - 0.2 * ( thetasold - thetas )
           thetas33 = thetas - theta33
         endif
         thetasdf = 1.-(density/2.65) ! density adjusted porosity
         porosity(i) = thetasdf
         lambda = 1. / (-potB(i))
         lambda_Saxton = lambda
         alpha = density / 2.65
         gravelvol(i) = ( alpha * gravelpc ) / ( 1. - gravelpc * ( 1. - alpha ) )
         SWCgravelfactor = 1. - gravelvol(i)
         bulkdensity = density * ( 1. - gravelvol(i) ) + ( gravelvol(i) * 2.65 )
         KbKs = ( 1. - gravelpc ) / ( 1. - gravelpc * ( 1. - 3. * alpha / 2. ) )
         cond1(i) = 1930. * ( thetas - theta33 )**( 3. - lambda ) ! saturated conductivity of matric soil
         if (gravelpc > 0.) cond1(i) = cond1(i) * KbKs
         cond2(i) = thetas ! saturated moisture, or porosity
       else
         if (Prades) gravelvol(i) = 0.34
       endif ! end of if condition for SaxtonNew

    enddo ! end of loop over soil layers
 
  end subroutine saxton_parameters
  !
  !-------------------------------------------------------------------------------------!
  !
  real function soil_conductivity( wf )

    ! Used in the soil drainage integrator. !
    ! Returns a single-point value.         !
    ! 'slayer' is a module variable that    !
    !  provides the soil-layer number.      !

    use soil_structure, only: cond1, cond2, cond3, SWP, SaxtonNew, KbKs, lambda_Saxton
    use enkf, only: Vallcebre, Prades

    implicit none

    ! arguments..
    real, intent(in) :: wf ! fraction of water in soils
    real :: alpha ! new Saxton parameter
    real :: wf_adj !adjusted SWC!

    wf_adj = wf
    !wf_adj = ( wf - 0.129 ) / 0.139 !wf * 0.139 + 0.129 ! bulky

    if ( wf_adj .lt. 0.001 ) then    ! Avoid floating-underflow fortran error
      soil_conductivity = 1e-30
    else
      soil_conductivity = cond1(slayer) * exp( cond2(slayer) + cond3(slayer) / wf_adj ) ! Soil conductivity (m s-1 )
      if (SaxtonNew) soil_conductivity = ( cond1(slayer) * ( wf_adj / cond2(slayer) )**( 3. + ( 2. / lambda_Saxton ) ) ) / 1000. / 3600.
      !if (Vallcebre) soil_conductivity = ( 10.**(-6) ) * ( (-SWP(slayer) )**(-0.092) )
    endif

  end function soil_conductivity
  !
  ! --------------------------------------------------------------------
  !
  subroutine soil_porosity()

   ! Porosity is estimated from Saxton equations. !

    use hydrol,             only: claypc, sandpc
    use scale_declarations, only: core
    use soil_structure,     only: porosity, SaxtonNew, gravelvol
    use enkf,               only: Vallcebre, Prades

    implicit none

    ! local variables..
    integer :: i
    real    :: H, J, K
    real,dimension(21) :: dummy

    ! saxton params relevant to porosity..
    H = 0.332  ;  J = -7.251e-4  ;  K = 0.1276

    if ( SaxtonNew .neqv. .true.) then ! porosity is calculated in SR saxton_parameters if new Saxton equations are used
      ! loop over soil layers..
      do i = 1 , core
        porosity(i) = H + J * sandpc(i) + K * log10( claypc(i) )
      end do
    endif
    
    !if (Vallcebre) porosity = 0.35 !/ (1. - gravelvol ) ! Vallcebre specific value!!! 0.49 and 0.465 = 0.4775
    if (Prades) then
       porosity = 0.65 !0.53 !/ ( 1. - gravelvol ) ! gravelvol = 0.34; 0.204 is observed max SWC for bulk (i.e. gravel corrected) soil
       !porosity( 15 : core ) = 0.
    endif

  end subroutine soil_porosity

  !
  ! --------------------------------------------------------------------
  !

  subroutine soil_resistance( A , m )

    !> subroutine description <!

    use meteo,              only: head
    use scale_declarations, only: pi,core
    use soil_structure,     only: abovebelow, conduc, rootl, rootlength, rootmass, &
                                  rootrad, rootresist, soilR, soilR1, soilR2, thickness
    use enkf
    use math_tools

    implicit none

    ! arguments
    real,dimension(ndim,nrens),intent(inout) :: A
    integer,intent(in)                       :: m

    ! local variables..
    integer :: i
    real    :: Lsoil, rs, rs2

    soilR = 0. ; soilR1 = 0. ; soilR2 = 0.

    ! Calculate soil-root hydraulic resistance
    do i = 1 , rootl
      Lsoil = conduc(i) / head    !converts from ms-1 to m2 s-1 MPa-1
      if ( Lsoil .lt. 1e-35 ) then    !prevent floating point error
        soilR(i) = 1e35
      else 
        rs  = sqrt( 1. / (rootlength(i) * pi ) )
        rs2 = log( rs / rootrad ) / ( 2.0 * pi * rootlength(i) * Lsoil * thickness(i) )
        soilR1(i) = rs2 * 1E-6 * 18. * 0.001    ! convert from MPa s m2 m-3 to MPa s m2 mmol-1
        !second component of below ground resistance related to root hydraulics
        !rrcheck=rootresist/(rootmass(i)*thickness(i)/abovebelow)
        soilR2(i) = rootresist / ( rootmass(i) * thickness(i) / abovebelow )
        soilR(i)  = soilR1(i) + soilR2(i)
      endif
    enddo

  end subroutine soil_resistance
  !
  ! --------------------------------------------------------------------
  !
  subroutine soil_water_potential()

    ! Find SWP without updating waterfrac yet (we do that in !
    ! waterthermal). Waterfrac is m3 m-3, soilwp is MPa.     !

    use soil_structure, only: potA, potB, rootl, SWP, waterfrac, SaxtonNew, theta33, psie, thetas
    use enkf, only: Vallcebre, Prades, defoliated, nondefoliated

    implicit none

    ! local variables..
    integer :: i
    real    :: soil_WP,dummy1,dummy2

    do i = 1 , rootl
      if ( waterfrac(i) .gt. 0. ) then
        soil_WP = -0.001 * potA(i) * waterfrac(i)**potB(i)   !  Soil water potential (MPa)
        if ( ( SaxtonNew ) .and. ( soil_WP > -0.033 ) ) soil_WP = -0.001 * ( 33. - ( ( waterfrac(i) - theta33 ) * ( 33. - psie ) / ( thetas - theta33 ) ) ) ! linear equation for high SWPs
        if ( ( SaxtonNew ) .and. ( soil_WP > ( psie * -0.001 ) ) ) soil_WP = 0.
        !if (Vallcebre) soil_WP = -0.001 * 0.0005 * waterfrac(i)**(-8.525)
        !if ( ( Prades ) .and. ( nondefoliated ) ) soil_WP = -1. * exp( 1. ) ** ( -1. * log( waterfrac( i ) ) * 0.759 - 2.006 ) ! bulky
        !write (*,*) soil_WP
        !stop
        !if ( ( Prades ) .and. ( defoliated ) )    soil_WP = -1. * exp( 1. ) ** ( -1. * log( waterfrac( i ) ) * 0.759 - 2.006 - 0.128 )
      else
        soil_WP = -9999.
      endif
      SWP(i) = soil_WP
    end do

  end subroutine soil_water_potential
  !
  !--------------------------------------------------------------------
  !
  subroutine water_uptake_layer( total_est_evap )

    ! Find from which layer water is withdrawn !

    use log_tools
    use scale_declarations, only: core, nos_canopy_layers
    use soil_structure,     only: fraction_uptake, iceprop, rootl, soilR, soilR1, soilR2, SWP, weighted_SWP
    use veg,                only: canopy_soil_resistance, lai, minlwp, CSR1, CSR2

    implicit none

    ! arguments..
    real,intent(out) :: total_est_evap

    ! local variables..
    integer              :: i, j, p
    real                 :: est_evap(core), frac
    integer, allocatable :: AR2(:) ! put values in array

    ! -- calculations begin below --

    j = size( shape( est_evap ) ) ! Get number of dimensions in array
    allocate( AR2(j) )            ! Allocate AR1 to number of dimensions in array

    total_est_evap = 0.  ;  weighted_SWP    = 0.
    est_evap       = 0.  ;  fraction_uptake = 0.

    do i = 1 , rootl
      ! estimate max transpiration from gradient-gravity / soil resistance.
      est_evap(i) = ( SWP(i) - minlwp ) / soilR(i)
      est_evap(i) = max( 0. , est_evap(i) )       ! no negative 
      if ( iceprop(i) .gt. 0. ) est_evap(i) = 0. ! no uptake from frozen soils
    enddo
    total_est_evap = sum( est_evap )

    ! weighted soil water potential
    if ( total_est_evap .gt. 0. ) then
      ! Water was evaporated from some layers..
      do i = 1 , rootl
        weighted_SWP = weighted_SWP + swp(i) * est_evap(i)
        ! fraction of total et taken from layer i...
        fraction_uptake(i) = est_evap(i) / total_est_evap
      enddo
      weighted_SWP = weighted_SWP / total_est_evap
    else
      ! No water was evaporated..
      fraction_uptake = 1. / rootl
    endif

    if ( nint ( sum( fraction_uptake ) ) .ne. 1. ) then
       call write_log( 'The sum of uptake fraction is not (nearly) equal to 1 '&
                     //' in water_uptake_layer' , msg_warning , __FILE__ , __LINE__ )
    endif
    if ( ( fraction_uptake(1) .gt. 1. ) .or. ( fraction_uptake(1) .lt. 0. ) ) then
       call write_log( 'Problem with the uptake fraction (either >1 or 0<)' , &
                       msg_warning , __FILE__ , __LINE__ )
    endif

    canopy_soil_resistance = 0.    ! reset
    CSR1 = 0. ; CSR2 = 0.
    do p = 1 , nos_canopy_layers
      frac = lai(p) / sum(lai)
      do i = 1 , rootl
        if ( frac .gt. 0. ) then 
          ! soil resistance for each canopy layer is related to leaf area
          canopy_soil_resistance(p) = canopy_soil_resistance(p) + 1. / ( soilR(i)  / frac )
          CSR2(p) = CSR2(p) + 1. / ( soilR2(i) / frac )
        else
          canopy_soil_resistance(p) = 0.001
        endif
      enddo
      canopy_soil_resistance(p) = 1. / canopy_soil_resistance(p)
      CSR2(p) = 1. / CSR2(p)
    enddo
    CSR1 = canopy_soil_resistance - CSR2

  end subroutine water_uptake_layer
  !
  !--------------------------------------------------------------------
  !
  subroutine water_retention()

    ! field capacity calculations for saxton eqns !

    use math_tools, only: zbrent
    use soil_structure

    implicit none

    ! local variables..
    integer        :: i
    real           :: x1, x2
    real,parameter :: xacc = 0.0001

    do i = 1 , core
      x1   = 0.1    ! low guess
      x2   = 0.7    ! high guess
      pass = i
      ! field capacity is water content at which SWP = -10 kPa
      field_capacity( i ) = zbrent( water_retention_saxton_eqns , x1 , x2 , xacc )
    enddo
    
  end subroutine water_retention
  !
  !--------------------------------------------------------------------
  !
  real function water_retention_saxton_eqns( xin )

    ! field capacity calculations for saxton eqns !

    use soil_structure, only: pass, potA, potB
    use enkf, only: Vallcebre, Prades

    implicit none

    ! arguments..
    real, intent(in) :: xin

    ! local variables..
    real ::soil_wp

    soil_WP = -0.001 * potA( pass ) * xin**potB( pass )   !  Soil water potential (MPa)
    !if (Vallcebre) soil_WP = -0.001 * 0.0005 * xin**(-8.525)

    water_retention_saxton_eqns = 1000. * soil_wp + 10.        ! 10 kPa represents air-entry swp

  end function water_retention_saxton_eqns
  !
  !--------------------------------------------------------
  !
end module soil_air
!
!--------------------------------------------------------
