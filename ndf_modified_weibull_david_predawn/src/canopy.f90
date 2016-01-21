! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !

module canopy

  !! >this module requires a summary< !!

  implicit none

contains

  subroutine timestep_calcs( A, B , m , time )

    ! CANOPY models (previously called day_sunshade_soilwater)
    !  runs first by timestep, and within each step then runs through each layer
    !soilwater version 18 feb 1999
    !  {calculate daily carbon gain for one layer}
    ! new parallel calculations of deltaLWP for sun and shade fractions (29 June 1998).

    use allocate_carbon,     only: predictor , phen_sim
    use clim,                only: coa, gbw, par_ratio, temp_bot, temp_top, wdbot, wdtop
    use hourscale,           only: totet
    use leaf,                only: assimilate
    use light,               only: checkpar, leaf_sun, longem, nirabs, parabs_shade, parabs_sun, solar
    use meteo,               only: gbb, la, nit, par, psil, psis, rad, temp, wdef
    use scale_declarations,  only: nos_canopy_layers, time_holder, steps
    use soil_air,            only: boundary
    use soil_functions,      only: soil_day
    use soil_structure,      only: weighted_SWP
    use veg,                 only: avn, co2amb, gppt, lai, lafrac, lma, lwpstore, predawn_lwp, nfrac, &
         nla, respt, totevap, totn, transt , totlai , flux, stom_conduct
    use enkf

    implicit none

    ! arguments..
    integer,intent(in) :: m ! ensemble member subscript
    real,dimension(ndim,nrens),intent(inout) :: A , B
    type(time_holder),intent(in) :: time

    ! local variables..
    integer :: i
    real    :: agr,frac_sh,frac_sun,ggppt,gsm,hold_psil,la_shade,la_sun,lambda, &
         modet,res,shade_psil,sun_psil,totestevap,weighted_gs(nos_canopy_layers)
    real    :: dayint ! function to be called
    real    :: dummy1
    logical :: dummy2
    external dayint
    real,dimension(96,10) :: delete
    logical :: sun 

    ! variables for calculating predawn LWP
    !real :: predawn_lwp
    logical :: set_predawn_lwp
    integer :: predawn

    predawn = steps / 24 * 6
    A( 22 , m ) = 0. ; A( 23 , m ) = 0. ! reset transpiration and capacitance for each timestep
    ! (fluxes are cumulative over (non-)shaded leaf area and canopy layers)

    if ( time%step .eq. 1 ) then ! reset values at first timestep of given day
       gppt   = 0.
       respt  = 0.
       transt = 0.
    endif
    set_predawn_lwp = .false.
    if ((time%step .eq. predawn) .or. (time%steps_count .lt. predawn)) then
       set_predawn_lwp = .true.
    endif
    lai    = lafrac * A( 8 , m ) / LMA   ! determine LAI distribution from Cf
    dummy1 = sum(lai)
    dummy1 = LMA
    if (phen_sim .eqv. .false.) lai = lafrac * totlai
    dummy1 = totlai
    dummy1 = sum(lai)
    dummy2 = phen_sim
    totn   = avN * sum(lai)        ! total N
    Nla    = nfrac * totn          ! N per layer
    call boundary()
    call solar( time )
    call soil_day( time , totestevap , A , m )
    psis = weighted_SWP                !weighted soil water potential
    ! weighted stomatal conductance for each layer, based on sun and shaded fractions - a diagnostic variable
    weighted_gs = 0.
    if (m == 1) stom_conduct = 0.
    do i = 1 , nos_canopy_layers
       psil      = LWPstore( i )     ! load layer's LWP from previous timestep
       !if (set_predawn_lwp) predawn_lwp = psil
       co2amb    = coa   ! load current CO2 concentration
       hold_psil = psil              ! store LWP 
       frac_sun  = leaf_sun( i )
       la_sun    = frac_sun * lai(i) ! first do sunlit leaf area
       la        = la_sun
       sun = .true.
       if ( ( la .gt. 0. ) .and. ( totestevap .gt. 0. ) ) then    !if there is active leaf area
          nit = frac_sun * Nla(i)
          par = ( parabs_sun(i) + parabs_shade(i) * frac_sun ) / la    ! sunlit leaf area has all beam rad plus frac of diffuse rad
          ! net radiation= shortwave + longwave radiation
          rad = 0.001 * ( nirabs(i) * frac_sun / la + longem(i) * frac_sun / la + par / par_ratio )
          temp = temp_top - 0.1111 * ( i - 1 ) * ( temp_top - temp_bot )
          wdef = wdtop - 0.1111 * ( i - 1 ) * ( wdtop - wdbot )
          gbb  = gbw(i)                    ! boundary layer
          call set_leaf( i , frac_sun , A , m , sun , predawn_lwp(i) )    ! see below
          call assimilate( i , time , modet , agr , res , gsm , A , m )    ! see leaf.f90
          gppt(time%step)   = gppt(time%step)   + agr     ! umol CO2 assimilation m-2 ground area s-1
          respt(time%step)  = respt(time%step)  + res     ! umol CO2 respiration m-2 ground area s-1
          transt(time%step) = transt(time%step) + modet   ! evapotranspiration in W m-2 ground area
          checkpar(i)     = checkpar(i) + par * la
          weighted_gs(i) = gsm * la / lai(i)
       endif
       sun_psil = psil         ! record LWP in sun leaves
       psil     = hold_psil    ! reset LWP for shade leaf calculation

       ! then do shaded leaf area
       frac_sh = ( 1. - leaf_sun(i) )
       la_shade = frac_sh * lai(i)
       la = la_shade
       sun = .false.
       if ( ( la .gt. 0. ) .and. ( totestevap .gt. 0. ) ) then
          nit = frac_sh * Nla(i)
          par = parabs_shade(i) * frac_sh / la
          ! net radiation= shortwave + longwave radiation
          rad = 0.001 * ( nirabs(i) * frac_sh / la + longem(i) * frac_sh / la + par / par_ratio )
          temp = temp_top - 0.1111 * ( i - 1 ) * ( temp_top - temp_bot )
          wdef = wdtop - 0.1111 * ( i - 1 ) * ( wdtop - wdbot )
          ! boundary layer
          gbb = gbw(i)
          call set_leaf( i , frac_sh , A , m , sun , predawn_lwp(i) )                    ! see below
          call assimilate( i , time , modet , agr , res , gsm , A , m )    ! see leaf.f90
          gppt(time%step)   = gppt(time%step)   + agr     ! umol CO2 assimilation m-2 ground area s-1
          respt(time%step)  = respt(time%step)  + res     ! umol CO2 respiration m-2 ground area s-1
          transt(time%step) = transt(time%step) + modet   ! evapotranspiration in W m-2 ground area
          checkpar(i)    = checkpar(i) + par * la
          weighted_gs(i) = weighted_gs(i) + gsm * la / lai(i)
       endif
       shade_psil = psil    ! record LWP in shade leaves
       psil = frac_sun * sun_psil + frac_sh * shade_psil    ! calculate final LWP
       if (set_predawn_lwp) predawn_lwp(i) = psil
       LWPstore(i) = psil
    enddo	! end loop over canopy layers
    ! ensemble average stomatal conductance per canopy layer
    stom_conduct = stom_conduct + weighted_gs
    if (m == nrens) stom_conduct = stom_conduct / nrens
    !total evapotranspiration during this time period
    !A( 22 , m ) = A( 22 , m ) * time%seconds_per_step / 1000. ! converted from g m-2 GA s-1 to kg m-2 GA ts-1; kg = liter = mm
    !A( 23 , m ) = A( 23 , m ) * time%seconds_per_step / 1000. ! converted from g m-2 GA s-1 to kg m-2 GA ts-1; kg = liter = mm
    lambda  = 1000. * ( 2501.0 - 2.364 * temp )    ! latent heat of vapourisation, J kg-1
    totet   = 0.001 * time%seconds_per_step * transt(time%step) / lambda    !m3 m-2 t-1
    totevap = totevap + 1000. * totet    !mm m-2, or litre m-2
    ggppt   = time%seconds_per_step * gppt(time%step) * 1e-6 * 12.    !g C assimilated per time step
    call predictor( A , m , ggppt , time )    !DALEC - C allocation and turnover

  end subroutine timestep_calcs
  !
  !--------------------------------------------------------------------
  !
  subroutine set_leaf( clayer , frac , A , m , sun , predawn_lwp )

    ! sets Farquhar parameters and hydraulic parameters for each leaf layer

    use metab, only: ht, layer_capac, propresp, rn, rplant, rsoil, vcm, vjm , &
      rplant_canopy, rsoil_canopy
    use meteo, only: la, nit, psil ! Olli added psil
    use veg,   only: canopy_soil_resistance, capac, conductivity, gplant, kappac, kappaj, layerht, rbelow
    use enkf
    use scale_declarations


    implicit none

    integer, intent(in) :: clayer
    real, intent(in)    :: frac
    real,dimension(ndim,nrens),intent(inout) :: A
    integer,intent(in) :: m ! ensemble member subscript
    logical,intent(in) :: sun
    real, intent(in) :: predawn_lwp
    real :: illum_frac,layer_frac_soil_resistance

    ! local variables added by Olli
    real :: d, c, plant_conduct_max=2 !added by Olli. This is David's value1.221185571
    
    ! choose one of the following two lines for large or small trees, added by Olli
    d = 0.8 !1.1878 ! large trees, added by Olli. These are David's values
    !d = 2.4 ! small trees, added by Olli.

    ! choose one of the following two lines for large or small trees, added by Olli
    !c = 6.26082 ! large trees, added by Olli
    c = 1 !1.6686 ! small trees, added by Olli. These are David's values
    

    vcm = kappac * nit / la
    ! metabolic rates are umol/m2/s - so we assume here that nit is N/m2 (it's actually N/clayer)
    vjm = kappaj * nit / la
    rn = 0.105 * propresp  ! respiration constant umol CO2/g N at 10 deg C
    rn = rn / la           ! convert to resp per m2

    ht = layerht( clayer )

    ! plant hydraulic resistances are determined by the amount of leaf area in the sunlit or shaded fraction
    if ( conductivity .eq. 1 ) then
       !rplant = ht / ( gplant * la )  ! MPa s mmol-1 (per layer); gplant is A( 25 ); deactivated by Olli
       rplant = 1./( plant_conduct_max * exp( -( -predawn_lwp/d)**c )*la) ! added by Olli
    else
       !rplant = 1. / ( gplant * la )  ! conductance is constant with height; gplant is A( 25 ); deactivated by Olli
       rplant = 1./( plant_conduct_max * exp( -( -predawn_lwp/d)**c )*la) ! added by Olli
    endif  
    
    layer_capac = capac * la

    rsoil = canopy_soil_resistance( clayer )  ! soil+root resistance of each canopy layer

    ! now correct rsoil according to the fraction of total leaf area in this run
    rsoil = rsoil / frac
    if (below_frac) rsoil = (rplant * rbelow) / (1. - rbelow)

!    if (sun) then
!      rplant_canopy( clayer ) = 1./rplant
!      if (frac == 1) then
!        rplant_canopy( clayer ) = 1./rplant_canopy( clayer )
!      endif
!      rsoil_canopy( clayer )  = rsoil
!    else
!      if (frac == 1) then
!        rplant_canopy( clayer ) = 0. ; rsoil_canopy( clayer ) = 0.
!      endif
!      rplant_canopy( clayer ) = rplant_canopy( clayer ) + 1./rplant
!      rplant_canopy( clayer ) = 1./rplant_canopy( clayer )
!      rsoil_canopy( clayer )  = rsoil_canopy( clayer ) + rsoil
!    endif

  end subroutine set_leaf
  !
  !--------------------------------------------------------------------
  !
end module canopy
!
!--------------------------------------------------------------------
!
