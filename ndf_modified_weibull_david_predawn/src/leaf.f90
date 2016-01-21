! ------------------------------------------------------- !
!                                                         !
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module leaf

  !!           > module descriptor here <             !!

  use scale_declarations, only: nos_canopy_layers
  use math_tools

  implicit none

  !integer,parameter,private :: kmaxx  = 200    ! descriptors
  !integer,parameter,private :: MAXSTP = 10000  ! would be good!

  real :: gamma1, & ! gamma arrhensis curves for modification of limitation on photosynthesis
             gs2, & ! common block stomata conductance
              lt, & ! leaf temperature for photosynthesis (oC)
              kc, & ! half saturation concentration for carboxylation rate (umol.m-2)
              ko, & ! half saturation concentration for oxygenation rate (umol.m-2)
           vcmax, & ! Temperature modified vcmax rates
           vjmax    ! Temperature modifed vjmax rates

  real,dimension(nos_canopy_layers) :: & 
          numsum, & ! count of the number of times a given canopy layer has been counted
                    ! i.e. nos of timesteps passed through where photosynthesis has taken place
            suma, & ! Sum of photosynthesis per canopy layer but without shade/light area correction.
                    !  Is reset daily and is intended for output / use in daily NEE sum values
          sumaci    ! Sum of photosynthesis per canopy layer with correction for leaf area in
                    !  shade / light (umolC.m-2.s-1 canopy_layer-1)

contains
  !
  !-------------------------------------------------------------------
  !
  subroutine assimilate( clayer , time , etr , agr , res , gsm , A , m )

    ! Assimilation is determined through equalising assimilation !
    !  according to the Farquhar model and a diffusion model,    !
    !  through varying internal CO2 concentration (Ci).          !
    ! Rsoil calculations have been moved to soil_air.            !
    ! switch variables.. ded to specify                          !
    !    -conductance vs conductivity                            !
    !    -Farquhar parameters are N-linked or not                !

    use clim,               only: atmos_press
    use math_tools,         only: zbrent
    use meteo,              only: gbb, la, nit, psil, rad, rcon, temp, wdef
    use metab,              only: an, ci, et, layer_capac, resp, rn
    use scale_declarations, only: boltz, time_holder
    use spa_io,             only: handle_output
    use veg,                only: co2amb, flux, lai
    use enkf

    implicit none

    ! arguments
    integer,intent(in)                       :: clayer , m
    type(time_holder),intent(in)             :: time
    real,intent(out)                         :: agr, etr, gsm, res
    real,dimension(ndim,nrens),intent(inout) :: A

    ! local variables..
    real :: ad, asn, check1, check2, conv, darkresp, dpsildt, g1, g2, gs, &
            lambda, lt, maxg, netrad, prevpsil, xacc, to_mm_ts

    ! reset..
    et = 0.  ;  asn = 0.  ;  ci = 0.  ;  an = 0.  ;  ad = 0.
    xacc  = 0.0001
    gs    = 0.00001

    ! Brent's method to determine gs (m s-1)
    g1 = 0.00005
    g2 = 0.05

    ! check if conditions allow photosynthesis
    check1 = stomdiff( g1 )
    check2 = stomdiff( g2 )
    maxg = 0.002 
    if ( check1 * check2 .lt. 0 ) then 
      if ( maxg .gt. 0.00002 ) then 
        gs = zbrent( stomdiff , g1 , g2 , xacc )
      else        ! drought limitation - psil has fallen to minpsil (=minlwp)
        gs = maxg
      endif

      prevpsil = psil
      call deltapotential( psil )    ! change in leaf water potential
      dpsildt = psil - prevpsil
      lt = leaftemp( gs )
      netrad = rad - 4. * 0.96 * boltz * 1d-3 * ( temp + 273.15 )**3 * ( lt - temp )
      et  = evap( gs , lt , netrad , wdef , gbb ) * 1.0 / 18.0
      res = resp * la            ! determine total leaf respiration
      agr = ( an + resp ) * la   ! convert from m-2 to clayer total gross rates (resp is related to [N] not LAI)
    else    ! dark
       gs = 0.00004
       ci = co2amb
       lt = leaftemp( gs )
       netrad =rad - 4. * 0.96 * boltz * 1d-3 * ( temp + 273.15 )**3 * ( lt - temp )
       et = evap( gs , lt , netrad , wdef , gbb ) * 1.0 / 18.0    ! convert from g m-2 s-1 to mol m-2 s-1
       darkresp = rn * nit * exp( log(2.0) * ( lt - 10.0 ) / 10.0 )
       an  = -darkresp
       res = darkresp * la    ! convert from m-2 to clayer total for net rates
       agr = 0.0
       prevpsil = psil
       call deltapotential( psil )  ! change in leaf water potential
       dpsildt = psil - prevpsil
    endif
    lambda = 0.001 * ( 2501.0 - 2.364 * temp )    ! (kJ g-1)
    etr = et * 18. * 1000. * lambda * la          ! convert et (mol m-2 s-1) to latent heat (Wm-2)
    to_mm_ts = time%seconds_per_step * 18. * .001 * .001 ! convert sapflow from mmol m-2 GA s-1 to mm m-2 GA ts-1

    ! water flux at base of trunk..
    A( 22 , m ) = A( 22 , m ) + ( la * et * 1000. * to_mm_ts ) + ( ( layer_capac * dpsildt / ( time%seconds_per_step ) ) * to_mm_ts ) !(mm m-2 GA ts-1)
    !A( 23 , m ) = A( 23 , m ) + layer_capac * dpsildt / ( time%seconds_per_step ) * to_mm_ts ! capacitance component of sap flow (mm m-2 GA ts-1)

    !A( 22 , m ) = A( 22 , m ) + et * 18. * la ! transpiration component of sap flow (g m-2 s-1)
    !A( 23 , m ) = A( 23 , m ) + layer_capac * dpsildt * 18. / ( time%seconds_per_step * 1000. ) ! capacitance component of sap flow (g m-2 s-1)
   ! write(*,*) (et * 18. * la) + (layer_capac * dpsildt * 18. / ( time%seconds_per_step * 1000. ))


    !flux( time%step , clayer ) = flux( time%step , clayer ) + A( 22 , m ) + A( 23 , m ) ! only for output

    conv = 1000. * atmos_press / ( ( lt + 273.15 ) * Rcon )    !convert from m s-1 to mmol m-2 s-1
    gsm  = gs * conv

    if ( agr .gt. 0. ) then    ! If photosynthesis is occurring
      ! correct for fraction of total lai in sun or shaded part
      sumaci( clayer ) = sumaci( clayer ) + an * la / lai( clayer ) * ci 
        suma( clayer ) =   suma( clayer ) + an * la / lai( clayer )
      numsum( clayer ) = numsum( clayer ) + 1.
    endif

    ! Produce output..
    !call handle_output( 3 , time , output_data = (/ real(clayer) , &
    !                                 gsm, agr, res, etr, lt, conv /) )

  end subroutine assimilate
  !
  !-------------------------------------------------------------------
  !
  pure function arrhenious( a , b , t )

    implicit none

    ! arguments..
    real,intent(in) :: a , b , t
    real            :: arrhenious

    ! local variables..
    real :: answer, denominator, numerator

    numerator   = t - 25.0d0
    denominator = t + 273.15d0

    answer = a * exp( b * 1d0 * numerator / denominator )
    arrhenious = 1d0 * answer

  end function arrhenious
  !
  ! --------------------------------------------------------------------
  !
  real function cdiff( xmid )

    ! difference between metabolic assimilation !
    ! rate and diffusive assimilation rate      !

    use metab, only: et, resp 
    use meteo, only: gbb, par 

    implicit none

    ! arguments..
    real,intent(in) :: xmid

    ! local variables..
    real :: adx, anx

    anx = farquhar( vcmax , vjmax , kc , ko , gamma1 , resp , par , xmid )   
    adx = diffusion( gs2 , xmid , et , gbb , lt )

    if ( anx .lt. 0. ) anx = 0.
    cdiff = adx - anx

  end function cdiff
  !
  !--------------------------------------------------------------------
  !
  real function comp( ccc , g1 )

    ! determine light-limitation to photosynthesis !

    use metab !, only: opttemp, rn, vckurtosis, vcm, vjkurtosis, vjm
    use meteo

    implicit none

    ! arguments..
    real,intent(in) :: ccc, g1

    ! local variables..
    real            :: gamma1,gs,kc,ko,lt,resx,vcmax,vcmt,vjmax,vjmt

    !  {temperature modifications begin below: based on air - not leaf - temperature}
    gs     = g1
    lt     = leaftemp( gs )
    vcmt   = tempmet( 65.03 , opttemp , Vckurtosis , lt )    !  {temperature modifications begin below}
    vjmt   = tempmet( 57.05 , opttemp , Vjkurtosis , lt )
    vcmax  = vcmt * vcm
    vjmax  = vjmt * vjm
    kc     = arrhenious( 310.0 , 23.956 , lt )
    ko     = arrhenious( 155.0 , 14.509 , lt )
    gamma1 = arrhenious(  36.5 ,  9.46  , lt )
    resx   = rn * nit * exp( log(2.0) * ( lt - 10.0 ) / 10.0 )

    comp = farquhar( vcmax , vjmax , kc , ko , gamma1 , resx , par , ccc )   

  end function comp
  !
  !-------------------------------
  !
  subroutine deltapotential( leafwp )

    ! change in leaf water potential !

    use math_tools !,         only: odeint,rkqs,kmaxx
    use scale_declarations, only: nmax

    implicit none

    ! arguments..
    real,intent(inout)  :: leafwp

    ! local variables..
    integer,parameter   :: nvar = 1
    integer             :: kmax , layer !, kount, nbad, nok
    real                :: eps, h1, hmin, x1, x2, ystart(nvar)

    !common /path/ kmax, kount, dxsav, x, y

    eps  = 1.0e-4
    h1   = 0.01
    hmin = 0.0
    kmax = 100
    x1   = 1.
    x2   = 2.

    dxsav     = ( x2 - x1 ) / 20.0
    ystart(1) = leafwp    !initial conditions

    call odeint( ystart , nvar , x1 , x2 , eps , h1 , hmin , nok , nbad , lwp_diff_eqn )

    leafwp = ystart( 1 )

  end subroutine deltapotential
  !
  !--------------------------------------------------------------------
  !
  real function diffusion( gs , ci , e , gbbb , ttemp )

    ! diffusion limited assimilation rate !

    use clim,  only: atmos_press
    use meteo
    use veg

    implicit none

    ! arguments..
    real,intent(in) :: ci , gs , e , gbbb , ttemp

    ! local variables..
    real            :: ad , convert , gt , ca

    ad = 0.0
    ca = co2amb
    !gi=0.0025    now in declarations !has been set to 0.003, originally set to 0.0025, see 5 July 97
    ! see Jones Appendix 3 for conversion ms-1 to mol s-1)
    ! and Jones Appendix 2 for ratio CO2 diffusion to water diffusion (=1.65)
    convert = atmos_press / ( Rcon * ( 273.15 + ttemp ) )
    gt = convert / ( 1.65 / gs + 1.37 / gbbb + 1. / gi )! total leaf conductance (converts from ms-1 to mols-1)
    ! interaction between water flux out of and CO2 flux into the stomata
    ! determines diffusion limited assimilation rate
    ad = ( gt - 0.5 * e ) * ca - ( gt + 0.5 * e ) * ci
    diffusion = ad

  end function diffusion
  !
  !--------------------------------------------------------------------
  !
  real function evap( gs , tt , q , wdef , gbb )

    !  determine evapotranspiration rate (g m-2 s-1) !
    !  from q (kW m-2), tt (oC), wdef (g m-3) and    !
    !  gbb and gs (m s-1)                            !

    implicit none

    ! arguments..
    real,intent(in) :: gbb, gs, q, tt, wdef

    ! local variables..
    real :: eps, gcut, gleaf, lambda, psych, s, slope

    gcut   = 0.0005
    gleaf  = gcut + gs
    ! slope of saturation vapour pressure curve (t-dependent)
    s      = 6.1078 * 17.269 * 237.3 * exp( 17.269 * tt / ( 237.3 + tt ) )
    slope  = 0.1 * ( s / ( 237.3 + tt )**2 )          ! (kPa K-1)
    psych  = 0.1 * ( 0.646 * exp( 0.00097 * tt ) )    ! psych is temp-dependent (kPa K-1)
    eps    = slope / psych                            ! response of epsilon to temp
    lambda = 0.001 * ( 2501.0 -2.364 * tt )           ! latent heat of vapourisation (KJ g-1)
    evap   = ( eps * q / lambda + wdef * gbb ) / ( 1.0 + eps + gbb / gs ) ! (g m-2 s-1)

  end function evap
  !
  !--------------------------------------------------------------------
  !
  real function farquhar( vcmax , vjmax , kc , ko , gamma1 , resp , par , ci )

    !  metabolic assimilation rate !

    implicit none

    ! arguments..
    real,intent(in) :: ci, gamma1, kc, ko, par, resp, vcmax, vjmax

    ! local variables..
    real :: alphaj, an, bee, oi, theta, vc, vj, wc, wj

    oi     = 210.0
    alphaj = 0.385
    theta  = 0.7
    an     = 0.0
    ! determine Rubisco limited carboxylation rate...
    wc  = ( vcmax * ci ) / ( ci + kc * ( 1.0 + oi / ko ) )
    ! determine potential rate of RuBP regeneration...
    bee = alphaj * par + vjmax
    vj  = ( bee - sqrt( bee**2 - 4.0 * theta * alphaj * par * vjmax ) ) / ( 2.0 * theta )
    ! determine RuBP regeneration limited carboxylation rate...
    wj  = vj * ci / ( 4.5 * ci + 10.5 * gamma1 )
    ! determine limiting carboxylation rate...
    vc  = min( wc , wj )
    ! net photosynthetic rate... 
    an  = vc * ( 1.0 - gamma1 / ci ) - resp

    farquhar = an

  end function farquhar
  !
  !--------------------------------------------------------------------
  !
  real function leafcp(gs)

    ! leaf compensation point !

    use metab, only: an,ci
    use meteo, only: psil

    implicit none

    ! arguments..
    real,intent(in) :: gs

    ! local variables..
    real :: dummy  ! intent-out from stomata, but then we never use?!

    call stomata( gs, ci, dummy )

    leafcp = an

  end function leafcp
  !
  !--------------------------------------------------------------------
  !
  real function leaftemp( gs )

    ! determines leaf temperature. !

    use meteo,              only: cp, gbb, rad, temp, wdef
    use log_tools
    use scale_declarations, only: boltz

    implicit none

    ! arguments..
    real,intent(in) :: gs

    ! local variables..
    real :: de, denom, diff, emiss, ghr, gr, lambda, psych, Q, rho, rhr, rt, s, slope, ta, tt

    tt    = temp
    q     = rad
    emiss = 0.96
    ta    = tt + 273.15
    rho   = 353000.0 / ta   ! density of air g m-3 (t-dependent)

    ! slope of saturation vapour pressure curve (t-dependent)..
    s     = 6.1078 * 17.269 * 237.3 * exp( 17.269 * tt / ( 237.3 + tt ) )
    slope = 0.1 * ( s / ( 237.3 + tt )**2 )
    psych = 0.1 * ( 0.646 * exp( 0.00097 * tt ) )       ! psych is temp-dependent 
    lambda = 0.001 * ( 2501.0 - 2.364 * tt )            ! latent heat of vapourisation

    ! convert water deficit to vapour pressure deficit..
    de    = wdef * lambda * psych / ( rho * cp )  
    gr    = 4. * emiss * boltz * 1d-3 * ta**3 / ( rho * cp )   ! remove leaf area sensitivity from gr
    ghr   = gr + 2. * gbb / 0.93                        ! total thermal conductance
    rhr   = 1. / ghr

    rt    = 1.0 / gs + 1.0 / gbb                        ! combined leaf resistance
    denom = psych * rt + slope * rhr
    diff  = rhr * rt * psych * q / ( rho * cp * denom ) - rhr * de / denom  ! temperature difference

    leaftemp = diff + tt          ! leaf temp
    if ( abs(diff) .gt. 10. ) then
      call write_log( 'leaf problem' , msg_warning , __FILE__ , __LINE__ )
    endif

  end function leaftemp
  !
  !--------------------------------------------------------------------
  !
  subroutine lwp_diff_eqn( time, y, dydt )

    ! differential equation describing change in !
    ! leaf water potential given supply & demand !

    use metab,             only: et,ht,layer_capac,rplant,rsoil
    use meteo,             only: la,head,psis
    use scale_declarations,only: nmax,steps

    implicit none

    ! arguments..
    real,intent(in)    :: y(nmax)
    real,intent(in)    :: time ! dummy argument, provided for odeint
    real,intent(out)   :: dydt(nmax)

    ! local variables..
    real :: secs

    secs    = 24. / steps * 3600.    ! Number of seconds per timestep
    dydt(1) = secs * ( psis - ( head * ht ) - 1000. * la * et * &
         ( rplant + rsoil ) - y(1) ) / ( layer_capac * ( rplant + rsoil ) )

  end subroutine lwp_diff_eqn
  !
  !--------------------------------------------------------------------
  !
  subroutine minstom( lowg , maxg )

    ! minimum gs for net C fixation, !
    ! max gs in drought situation.   !

    use math_tools, only: zbrent
    use metab
    use meteo
    use veg

    implicit none

    ! arguments..
    real,intent(out) :: lowg, maxg

    ! local variables..
    real :: low, lwp, high, mings, xacc

    an   = 0.
    low  = 0.00002
    high = 0.05
    maxg = high
    xacc = 0.00000001
    lwp  = psil
    call stomata( low, ci, lwp )
    if ( (an .gt. 0.) .or. (temp .lt. 5.) ) then
       mings = 0.00005
    else
       mings = zbrent( leafcp, low, high, xacc )  !add delta
       mings = mings + 0.00004
    endif
    ! check lwp limit
    lwp = psil
    call stomata( mings, ci, lwp )
    if ( lwp .lt. minlwp ) then
    maxg = low
    endif

    lowg = mings

  end subroutine minstom
  !
  !--------------------------------------------------------------------
  !
  subroutine stomata( gs , ci_dummy , lwp )

    ! determines stable ci for given gs !

    use math_tools,         only: zbrent
    use metab,              only: an, ci, et, opttemp, resp, rn, vckurtosis, &
                                  vcm, vjkurtosis, vjm
    use meteo,              only: gbb, nit, par, rad, temp, wdef
    use scale_declarations, only: boltz
    use veg,                only: co2amb

    implicit none

    ! arguments..
    real,intent(in)  :: gs
    real,intent(out) :: ci_dummy, lwp  ! had to change ci to ci_dummy as ci is declared in metab

    ! local variables..
    real             :: ad, netrad, vcmt, vjmt, x1, x2, xacc

    gs2    = gs            ! swap argument gs into common block gs2

    lt     = leaftemp( gs )

    ! estimate net radiation from net isothermal and temperature difference (Jones, p.108)..
    netrad = rad - 4. * 0.96 * boltz * 1d-3 * ( temp + 273.15 )**3 * ( lt - temp )

    ! evaporation rate in g m-2 s-1, converted to clayer total mol m-2 s-1...
    et     = evap( gs , lt , netrad , wdef , gbb ) * 1.0 / 18.0 
    et     = max( 0. , et )

    vcmt  = tempmet( 65.03 , opttemp , Vckurtosis , lt )    !  {temperature modifications begin below}
    vjmt  = tempmet( 57.05 , opttemp , Vjkurtosis , lt )
    vcmax = vcmt * vcm
    vjmax = vjmt * vjm

    kc     = arrhenious( 310.0 , 23.956 , lt )
    ko     = arrhenious( 155.0 , 14.509 , lt )
    gamma1 = arrhenious(  36.5 ,  9.46 ,  lt )

    resp = rn * nit * exp( log( 2.0 ) * ( lt - 10.0 ) / 10.0  )

    x1   = co2amb
    x2   = 1.0
    xacc = 0.01
    ci   = zbrent( cdiff , x1 , x2 , xacc )

    ! NOW UPDATE THINGS OUTSIDE OF THE MODULE...
    an   = farquhar( vcmax , vjmax , kc , ko , gamma1 , resp , par , ci )   ! stable ci C uptake
    ad   = diffusion( gs , ci , et , gbb , lt)

    ! change in leaf water potential
    call deltapotential( lwp )

    ci_dummy = ci

  end subroutine stomata
  !
  !--------------------------------------------------------------------
  !
  real function stomdiff( gs )

    ! efficiency check and caviation check to determine maximum gs !

    use metab, only: an 
    use meteo, only: psil 
    use veg,   only: iota, minlwp 

    implicit none

    ! arguments..
    real,intent(in) :: gs

    ! local variables..
    real            :: an1, an2, ci1, ci2, delta, eff, lwp, minpsi

    delta = 0.00003
    lwp   = psil                       ! hold initial psil value
    call stomata( gs-delta, ci2, lwp )
    an2 = an                           ! uptake at lower gs
    lwp = psil                         ! reinitialise psil value
    call stomata( gs, ci1, lwp )
    an1 = an                           ! uptake at gs
    eff = an1 - an2 - ( iota - 1. )    ! efficiency check 

    minpsi = lwp - minlwp              ! cavitation check

    stomdiff = min( eff, minpsi )

  end function stomdiff
  !
  !-------------------------------------------------------------------
  !
  real function tempmet(max,opt,q,x)

    ! > function summary? < !

    implicit none

    ! arguments..
    real,intent(in) :: max,opt,q,x

    ! local variables..
    real            :: dum

    if (x.ge.max)then
       tempmet=0.0
    else
       dum=(max-x)/(max-opt)
       dum=exp(log(dum)*q*(max-opt))
       tempmet=dum*exp(q*(x-opt))
    endif

  end function tempmet
  !
  !-------------------------------------------------------------------
  !
end module leaf
!
!-------------------------------------------------------------------
!
