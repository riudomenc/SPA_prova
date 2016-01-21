! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module light

  !! SPA DALEC hourly version                                             !!
  !! correction mar 03 to sun and shade routines                          !!
  !! version 4 allows input of diffuse fraction data                      !!
  !! corrected 9-10 April 1999 (sun and shade component of ATTENUATE_PAR) !!
  !! modified SOLAR, to partition absorption between                      !!
  !! sunlit and shaded fractions in each layer                            !!
  !! new frac_diffuse_rad - hourly variation in diffuse fraction          !!
  !! version 4, 3 Dec 1997                                                !!
  !! NIRABS is now truly just near-infra-red, not total short-wave as in  !!
  !! previous versions                                                    !!
  !! correction 15 Oct 1998 - set intbm and intdf =0 at start of          !!
  !!  attentuate_long and attenute_nir routines (par was correct already) !!

  use scale_declarations, only: nos_canopy_layers

  implicit none

  real :: fdiff,    & ! ?
          skyabsp,  & ! ?
          soilabsn, & ! ?
          soilabsp, & ! ?
          soilabsl    ! ?

  real,dimension(nos_canopy_layers) :: &
         abspar_sun,   & ! ?
         abspar_shade, & ! ?
         checkpar,     & ! Check PAR radiation balances
         leaf_sun,     & ! Fraction of leaf in a given canopy layer exposed to sun
         longem,       & ! long wave emission from each canopy layer
         nirabs,       & ! NIR abosrbed by each canopy layer
         parabs_shade, & ! PAR absorbed by each canopy layer in sun exposed leaves
         parabs_sun      ! PAR absorbed by each canopy layer in shade covered leaves

  real,parameter :: nirrefl   = 0.43  ! ?
  real,parameter :: nirtrans  = 0.26  ! ?
  real,parameter :: parrefl   = 0.11  ! Baldocchi
  real,parameter :: partrans  = 0.16  ! ?
  real,parameter :: soilpar   = 0.033 ! ?
  real,parameter :: soilnir   = 0.023 ! ?
  real,parameter :: spherical = 2.    ! spherical leaf angle distribution has a value 2.

contains
  !
  !-------------------------------------------------------------------
  !
  subroutine attenuate_NIR(ff,lai,bm,df,refl,trans,sum,reflsoil,absorb,updf, &
       downdf,soilrad,skyrad,kbm,clump)
    use scale_declarations
    implicit none

    ! arguments..
    real,intent(in)    :: bm,clump,df,ff,kbm,lai(nos_canopy_layers),refl,reflsoil,trans
    real,intent(inout) :: absorb(10),downdf(90:101),skyrad,soilrad,updf(90:101)
    real,intent(out)   :: sum

    ! local variables..
    integer :: i
    real    :: beam(90:101),decay,intbm,intdf,kdf,leafrad

    ! calculations..
    intbm=0.;intdf=0.
    beam=0.
    beam(101)=bm
    downdf(101)=df
    ! diffuse radiation approximated by beta=30 degrees
    kdf=0.5
    do i=100,91,-1
       ! determine unintercepted radiation
       beam(i)=beam(i+1)*exp(-kbm*clump*lai(101-i))
       ! and thus intercepted radiation
       intbm=beam(i+1)-beam(i)
       ! now for diffuse
       decay=exp(-kdf*clump*lai(101-i))
       downdf(i)=downdf(i)+downdf(i+1)*decay
       intdf=downdf(i+1)*(1.-decay)
       ! correct for transmittance (trans) and reflectance (refl)
       absorb(i-90)=absorb(i-90)+intbm*(1.-trans-refl)/ff
       absorb(i-90)=absorb(i-90)+intdf*(1.-trans-refl)/ff
       ! add transmitted beam & diffuse to downward diffuse
       downdf(i)=downdf(i)+trans*(intbm+intdf)
       ! add reflected beam & diffuse to upward diffuse
       updf(i)=updf(i)+refl*(intbm+intdf)
    enddo
    ! reflectance by soil of radiation that penetrates all vegetation
    updf(90)=reflsoil*(beam(91)+downdf(91))
    soilrad=soilrad+(1.-reflsoil)*(beam(91)+downdf(91))
    !     write(31,'(I4,14(",",F7.2))')t,(beam(i)+downdf(i),
    !    1      i=100,91,-1),bm,df
    do i=91,101
       downdf(i)=0.
    enddo
    do i=91,100
       ! now return upwards through the canopy, dealing with refelcted radiation
       ! unintercepted
       decay=exp(-kdf*clump*lai(101-i))
       updf(i)=updf(i)+updf(i-1)*decay
       ! intercepted
       intdf=updf(i-1)*(1.-decay)
       ! absorbed determined from transmittance/reflectance
       absorb(i-90)=absorb(i-90)+intdf*(1.-trans-refl)/ff
       ! add reflected beam & diffuse to downward diffuse
       downdf(i)=downdf(i)+refl*intdf
       ! add transmitted beam & diffuse to upward diffuse
       updf(i)=updf(i)+trans*intdf
       updf(i-1)=0.
    enddo
    skyrad=skyrad+updf(100)
    updf(100)=0.
    leafrad=0.
    do i=1,10
       leafrad=leafrad+absorb(i)
    enddo
    sum=soilrad+skyrad+leafrad

  end subroutine attenuate_NIR
  !
  !-------------------------------------------------------------------
  !
  subroutine attenuate_PAR( ff , lai , bm , df , refl , trans , sum_rad ,     &
                            reflsoil , absorb_sun, absorb_shade , updf ,      &
                            downdf , soilrad , skyrad , kbm , sunfrac , clump )

    use scale_declarations, only : nos_canopy_layers

    implicit none

    ! arguments..
    real,intent(in)    :: clump, bm, df, ff, kbm, lai(nos_canopy_layers), refl, reflsoil, trans
    real,intent(inout) :: absorb_sun(nos_canopy_layers), absorb_shade(nos_canopy_layers), skyrad, &
                          soilrad, updf(90:101)
    real,intent(out)   :: downdf(90:101), sum_rad, sunfrac(nos_canopy_layers)

    ! local variables..
    integer :: i
    real    :: beam(90:101), cumlai, decay, intbm, intdf, kdf, leafrad, suncum, sunla(nos_canopy_layers), sunprev

    ! calculations..
    beam = 0.  ;  intbm = 0.  ;  intdf = 0.  ;  sunla = 0.  ;  sunfrac = 0.
    beam(101)   = bm
    downdf(101) = df
    cumlai = 0.  ;  sunprev = 0.
    ! diffuse radiation approximated by beta=30 degrees
    kdf = 0.5
    do i = 100 , 91 , -1    ! determine unintercepted radiation
      if ( ( kbm .gt. 0. ) .and. ( bm .gt. 0. ) ) then ! determine sun-lit leaf area 
        cumlai = cumlai + lai( 101 - i )               ! first calculate cumulative leaf area
        suncum = ( 1 - EXP( -kbm * cumlai ) ) / kbm    ! then total sunlit for cumlai. kbm = 1/(2sinBeta)
        sunla(101-i) = suncum - sunprev                ! sunlit area in this layer is cumulative sunlit minus previous
        sunprev = suncum                               ! save previous
        if ( sunla( 101 - i ) .gt. 0. ) then           ! determine sunlight fraction
          sunfrac( 101 - i ) = sunla( 101 - i ) / lai( 101 - i )
        else
          sunfrac( 101 - i ) = 0.
        endif
        beam( i ) = bm    ! Beam radiation on sunlit leaves is constant 
        ! Intercepted radiation dI = Io.k.Ls
        ! (dI=intercepted rad, Io=downwelling beam, Ls =sunlit leaf area, k=extinction)
        intbm = bm * kbm * sunla( 101 - i )
      else
        intbm = 0.
      endif
      ! now for diffuse
      decay     = exp( -kdf * clump * lai(101-i) )            !attenuation factor
      downdf(i) = downdf(i) + downdf(i+1) * decay        !diffuse radiation that passes through layer
      intdf     = downdf(i+1) * ( 1. - decay )                !intercepted diffuse radiation in each layer
      ! Absorption of direct beam (correct interception for transmittance and reflectance)..
      absorb_sun( i - 90 )   = absorb_sun( i - 90 )   + intbm * ( 1. - trans - refl ) / ff
      ! Absorption of diffuse..
      absorb_shade( i - 90 ) = absorb_shade( i - 90 ) + intdf * ( 1. - trans - refl ) / ff

      ! add transmitted beam & diffuse to downward diffuse
      downdf( i ) = downdf( i ) + trans * ( intbm + intdf )
      ! add reflected beam & diffuse to upward diffuse
      updf( i )   = updf( i ) + refl * ( intbm + intdf )
    enddo
    if ( bm .gt. 0. ) then
      ! Direct beam radiation reaching soil surface is directly calculated by this eqn...
      beam( 91 ) = ( 1. / kbm - suncum ) * bm * kbm
    endif
    ! reflectance by soil of radiation that penetrates all vegetation
    updf(90) = reflsoil * ( beam(91) + downdf(91) )
    soilrad  = soilrad + ( 1. - reflsoil) * (beam(91) + downdf(91) )
    do i = 91 , 101
      downdf(i) = 0.
    enddo
    do i = 91 , 100
      ! now return upwards through the canopy, dealing with reflected radiation
      ! unintercepted
      decay = exp( -kdf * clump * lai(101-i) )
      updf(i) = updf(i) + updf(i-1) * decay
      ! intercepted
      intdf = updf( i - 1 ) * ( 1. - decay )
      ! absorbed determined from transmittance/reflectance
      absorb_shade(i-90) = absorb_shade(i-90) + intdf * ( 1. - trans - refl ) / ff
      ! add reflected beam & diffuse to downward diffuse
      downdf(i) = downdf(i) + refl * intdf
      ! add transmitted beam & diffuse to upward diffuse
      updf(i)   = updf(i) + trans * intdf
      updf(i-1) = 0.
    enddo
    skyrad    = skyrad + updf(100)
    updf(100) = 0.
    leafrad   = 0.
    do i = 1 , 10
      leafrad = leafrad + absorb_sun(i) + absorb_shade(i)
    enddo
    sum_rad = soilrad + skyrad + leafrad

  end subroutine attenuate_PAR
  !
  !-------------------------------------------------------------------
  !
  real function frac_diffuse_rad( dec , lat , parsteps , rtime )
    
    ! Determines the ratio of actual-to-potential radiation (ksko)   !
    ! for the day, from the total potential daily PAR, and then uses !
    ! a relationship from Erbs et al (1982) to estimate fraction of  !
    ! incoming radiation that is diffuse.                            !

    use scale_declarations, only: pi, steps

    implicit none

    ! arguments..
    real,intent(in) :: dec, lat, parsteps, rtime

    ! local variables..
    real :: ff, hourangle, ko, ks, ksko, So

    ! calculations..
    So        = 1360.                ! solar constant (Wm-2)
    ks        = parsteps / 2.3       ! (MJ m-2 d-1)
    hourangle = pi / 180. * 15. * ( rtime - 0.5 * real(steps) ) * 24. / real(steps)
    ! extraterrestrial radiation (MJ m-2)
    ko        = So * ( sin( dec ) * sin( lat ) + cos( lat ) * cos( dec ) * cos( hourangle ) )
    ksko      = ks / ko  ! hourly ksko
    ! Diffuse fraction
    if ( ksko .lt. 0.22 ) then
      ff = 1.0 - 0.09 * ksko
    else if ( ksko .lt. 0.8 ) then
      ff = 0.9511                  &
          - 0.1604 * ksko          &
           + 4.388 * ksko**2       &
            - 16.638 * ksko**3     &
             + 12.336 * ksko**4
    else 
      ff = 0.165
    endif
    frac_diffuse_rad = ff 

  end function frac_diffuse_rad
  !
  !-------------------------------------------------------------------
  !
  subroutine longwave( ff, lai, df, sumrad, absorb, updf, downdf, soilrad, skyrad, count, & 
                       Ta, Ts, totemit, soilem, totab, clump)

    use scale_declarations, only: boltz, nos_canopy_layers

    implicit none

    ! arguments..
    integer,intent(in) :: count
    real,intent(in)    :: clump, df, ff, lai(nos_canopy_layers), Ta(nos_canopy_layers), Ts
    real,intent(inout) :: absorb(10), downdf(90:101), skyrad, soilrad, totemit, updf(90:101)
    real,intent(out)   :: soilem, sumrad, totab

    ! local variables..
    integer :: i
    real    :: decay, emiss, eup, edown, intdf, kdf, leafrad


    intdf       = 0.
    emiss       = 0.96
    downdf(101) = df

    ! diffuse radiation approximated by beta=30 degrees
    kdf=0.5
    do i = 100 , 91 , -1
       decay = exp( -kdf * clump * lai(101-i) )
       if ( count .eq. 0 ) then
         ! longwave radiative heat loss from top side of leaf (KW m-2)
         eup = 0.001 * emiss * boltz * ( Ta(101-i) + 273.15 )**4 * ( 1. - decay )
            !'*(1-decay)' corrects emissions into the 1-D, vertical
         totemit = totemit + 2. * eup * lai(101-i)
       else
         eup = 0.
       endif

       ! and bottom side..
       edown     = eup
       downdf(i) = downdf(i) + downdf(i+1) * decay + edown
       intdf     = downdf(i+1) * ( 1. - decay )

       ! correct for transmittance (trans) and reflectance (refl) & PAI/LAI ratio
       absorb(i-90) = absorb(i-90) + intdf * emiss / ff - eup - edown

       ! add transmitted diffuse to downward diffuse
       downdf(i) = downdf(i) + 0.5 * ( 1. - emiss ) * intdf

       ! add reflected diffuse to upward diffuse
       updf(i) = updf(i) + 0.5 * ( 1. - emiss ) * intdf + eup
       totab = totab + intdf * emiss
    enddo

    ! reflectance by soil of radiation that penetrates all vegetation
    soilem = emiss * boltz * ( Ts + 273.15 )**4 * 0.001
    if ( count .eq. 0 ) then
       updf(90) = 0.04 * downdf(91) + soilem
    else
       updf(90) = 0.04 * downdf(91)
    endif
    soilrad = soilrad + 0.96 * downdf(91)
    do i = 91 , 101
       downdf(i) = 0.
    enddo
    do i = 91 , 100
       ! now return upwards through the canopy, dealing with refelcted radiation
       ! unintercepted
       decay   = exp( -kdf * clump * lai(101-i) )
       updf(i) = updf(i) + updf(i-1) * decay

       ! intercepted
       intdf = updf(i-1) * ( 1. - decay )

       ! absorbed determined from transmittance/reflectance
       absorb(i-90) = absorb(i-90) + intdf * emiss / ff

       ! add reflected beam & diffuse to downward diffuse
       downdf(i) = downdf(i) + 0.02 * intdf

       ! add transmitted beam & diffuse to upward diffuse
       updf(i)   = updf(i) + 0.02 * intdf
       updf(i-1) = 0.
       totab     = totab + intdf * 0.96
    enddo
    skyrad = skyrad + updf(100)
    updf(100) = 0.
    leafrad = 0.
    do i = 1 , 10
       leafrad = leafrad + absorb(i)
    enddo
    sumrad = soilrad + skyrad + leafrad

  end subroutine longwave
  !
  !----------------------------------------------------------------------------
  !
  subroutine solar( time )

    use clim,                only: par_ratio, par_top, sw_rad, temp_bot, temp_top 
    use irradiance_sunshade
    use scale_declarations,  only: boltz, nos_canopy_layers, pi, steps, time_holder
    use spa_io,              only: handle_output
    use veg,                 only: lai, lat, modrnet  

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time

    ! local variables..
    integer :: counter, i
    real    :: ab, beam, clump,dayd,dec,diff,em,estfdiff,ff,kbm,lnet,long, &
         nirbeam,nirdiff,period,PIb,PId,radnet,rtime,skyabsl,skyabsn,soem,   &
         soilt,suml,sumn,sump,sun,SWR,totlong,totnir,totpar,xhour
    real,dimension(90:101) :: downlong,downnir,downpar,uplong,upnir,uppar
    real                   :: temp(10)
    real,dimension(nos_canopy_layers)   :: abslong,absnir,sunfrac,totparabs

    ! calculations..
    period    = real(steps)
    clump     = 1.    ! clumping factor - =1 means unclumped, random distribution
    ff        = 1.0   ! correction factor for stems etc on attenuation, currently unused and untested, set at 1.0.
    dayd      = real( time%day )
    totparabs = 0.
    rtime     = real( time%step )

    ! solar characteristics..
    dec = -23.4 * cos( ( 360. * ( dayd + 10. ) / 365. ) * pi / 180. ) * pi / 180.
    ! diffuse fraction of radiation (Amthor,1994)

    ! Estimate diffuse fraction of radiation from measured vs maximum radiation..
    estfdiff = frac_diffuse_rad( dec , lat , par_top , rtime )
    ! day is divided into 48 1/2hour-long periods
    fdiff = estfdiff

    ! reset diffuse arrays
    do i = 90 , 101
      downnir(i)  = 0.
      upnir(i)    = 0.
      downpar(i)  = 0.
      uppar(i)    = 0.
      downlong(i) = 0.
      uplong(i)   = 0.
    enddo
    do i = 1 , 10
      absnir(i)       = 0.
      abspar_sun(i)   = 0.
      abspar_shade(i) = 0.
      abslong(i)      = 0.
      temp( i ) = temp_top - 0.1111 * (i-1) * ( temp_top - temp_bot )
    enddo

    ! detemine light extinction coefficent from sun angle..
    xhour = pi / 180. * 15. * ( rtime - 0.5 * period ) * 24. / period    ! use period = real(steps)
    ! sin(beta) - solar geometry affects attenuation
    sun = cos(lat) * cos(xhour) * cos(dec) + sin(lat) * sin(dec)

    ! determine extinction coefficient for direct beam radiation..
    if ( sun .lt. 0.06 ) then
      ! At low sun angles extinction coefficient is zero..
      ! (only occurs near sunrise and sunset - without this correction we get v unrealistic estimates)
      sun = 0.
      kbm = 0.0
    else
      kbm = 1. / ( spherical * sun )
    endif

    if ( sun .gt. 0. ) then
      ! If sun has risen then partition par between beam and diffuse
      beam = ( 1. - fdiff ) * par_top    ! light attentuation(PPFD) initialise layer above canopy
      diff =    fdiff       * par_top
    else 
      beam = 0.
      diff = par_top
    endif

    ! you can estimate short wave radiation from PAR data - 4.6 umol J-1, 50% of SWR is PAR
    !    new 19/Mar/97, i.e.  =par_top(t)/(0.5*par_ratio)    
    SWR     = sw_rad
    PIb     = 0.48 * ( 1. - fdiff ) * SWR
    PId     = 0.5 * SWR - PIb
    nirbeam = ( 1. - fdiff ) * SWR - PIb
    nirdiff = fdiff * SWR - PId

    ! determine longwave incident from sky (Jones p.27) kW m-2 long=.350..
    long = 0.001 * boltz * ( temp_top + 273.15 - 20. )**4

    totnir  = nirbeam + nirdiff
    totpar  = beam + diff
    totlong = long

    ! set values to zero
    soilabsn = 0.  ;  soilabsp = 0.  ;  skyabsp = 0. ;  skyabsn = 0.
    soilabsl = 0.  ;  skyabsl  = 0.  ;  counter = 0  ;  em = 0.  ;  ab = 0.

    ! start the attenuation routines..
    do while ( counter .lt. 3 )
      ! multiple passes through the canopy - 3 is enough to account for every photon in general
      if ( totpar .gt. 1. ) then    ! if the sun is up, then call the routines....

        ! firstly calculate PAR attenuation - only do sunlit and shaded calculations for PAR..
        call attenuate_PAR( ff , lai , beam , diff , parrefl , partrans , sump ,    &
                            soilpar , abspar_sun , abspar_shade , uppar , downpar , &
                            soilabsp , skyabsp , kbm , sunfrac , clump              )
        if ( counter .eq. 0. ) then ! only set leaf_sun after first pass, when beam radiation is incident
          do i = 1 , 10
            leaf_sun( i ) = sunfrac( i )
          enddo
        endif

        ! next do Near Infra Red..
        call attenuate_NIR( ff , lai , nirbeam , nirdiff , nirrefl , nirtrans , sumn , soilnir, &
                            absnir , upnir , downnir , soilabsn , skyabsn , kbm , clump         )
      else
        abspar_sun = 0.  ;  abspar_shade = 0.
        absnir     = 0.  ;  sunfrac      = 0.
      endif

      soilt = temp_top    ! use air temp from current timestep

      call longwave( ff , lai , long , suml , abslong , uplong , downlong , &        ! longwave
                     soilabsl , skyabsl , counter , temp , soilt , em , soem , ab , clump )
      beam = 0.  ;  diff   = 0.  ;  nirbeam = 0.  ;  nirdiff = 0.
      long = 0.  ;  radnet = 0.  ;  lnet    = 0.
      counter = counter + 1 
    enddo

    ! calculate output for use in other routines
    do i = 1 , 10
      parabs_sun(i)   = abspar_sun(11-i)      ! beam par absorbed in each layer
      parabs_shade(i) = abspar_shade(11-i)    ! diffuse par absorbed in each layer
      ! total shortwave radiation for heating = par + nir , per m2 ground area per layer
      nirabs(i)    = absnir(11-i)
      longem(i)    = abslong(11-i) * 1000.
      lnet         = lnet + longem(i)
      radnet       = radnet + nirabs(i) 
      totparabs(i) = abspar_sun(11-i) + abspar_shade(11-i)
    enddo

    ! net radiation equals shortwave + longwave radiation
    soilnet = soilabsn + 1000. * soilabsl + soilabsp / par_ratio

    check = sum(abspar_sun) + sum(abspar_shade) + soilabsp + skyabsp
    if ( par_top .gt. 0 ) then
      !call handle_output( 4 , output_data =  (/ soilnet, soilabsn, soilabsl,  &
      !                  soilabsp, par_top, fdiff, check, skyabsp, abspar_sun, &
      !                  parabs_sun, abspar_shade, parabs_shade, leaf_sun /)   )
    endif
    modRnet  = SWR - skyabsn - skyabsp / par_ratio - 1e3 * ( skyabsl - totlong )
    checkpar = 0.

  end subroutine solar
  !
  !-----------------------------------------------------------------
  !
end module light
!
!-----------------------------------------------------------------
!
