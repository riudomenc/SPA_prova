! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module allocate_carbon

  !! Carbon water model for data assimilation !!
  !! part of the DALEC model.                 !!
  !!                                          !!
  !! 'N' - N-limitation added to v3           !!
  !! Cfmax - max foliar carbon (LAI)          !!
  !! temperature controlled phenology         !!

  implicit none

  logical :: phen_sim = .false. ! phenology simulated or prescribed?

contains

  subroutine predictor( A , m , gppsum , time )    ! model
    !
    ! Predictor adjusts calculations depending
    ! on whether evergreen or deciduous run.

    use clim,               only: avtemp, daytempsum, gdd, max_fol, mint, multtf, multtl, temp_top
    use scale_declarations, only: steps, time_holder
    use soil_structure,     only: resprate, rootbiomass 
    use spa_io,             only: handle_output,user_opts
    use veg,                only: theta , totlai , LMA , root_leaf_ratio !, A
    use enkf

    implicit none

    ! arguments
    real,dimension(ndim,nrens),intent(inout) :: A
    integer,intent(in)                       :: m ! ensemble member subscript
    real,intent(in)                          :: gppsum
    type(time_holder),intent(in)             :: time ! hour of day

    ! local variables..
    real    :: Aar,npp,ralabfrom,ralabto,rate,ts_length

    ! Evergreen => 17 SVs = Ra,Af,Aw,Ar,Lf,Lw,Lr,Cf,Cw,Cr,Rh1,Rh2,D,Clit,Csom,G,Caresp
    ! Deciduous => 21 SVs = Ra,Af,Aw,Ar,Lf,Lw,Lr,Cf,Cw,Cr,Rh1,Rh2,D,Clit,Csom,G,
    !                                                   ET,Clab,Atolab,Afromlab,Caresp

    ts_length = 24. / real(steps)    ! length of time step in hours

    ! Carbon fluxes
    A(16,m) = gppsum
    rate  = 0.5 * exp( resprate * temp_top )

    if ( user_opts%veg_is_deciduous ) then

       if ( time%day .le. 100 ) then
          gdd     = 0.    ! growing degree day summation starts after day 100
          max_fol = 1.    ! switch to determine if max foliage has been attained
       endif
       if ( time%step .eq. 1 ) then   ! reset at start of each day
          daytempsum = 0.
          mint       = 100
       endif
       ! time switches
       mint = min( temp_top , mint )
       daytempsum = daytempsum + temp_top
       if ( time%step .eq. steps ) then          ! determine average temp
          avtemp = daytempsum/real(steps)
          gdd    = gdd+avtemp        ! growing degree day heat sum 
          !call handle_output( 2 , time ) ! write this to output
       endif
       if (gdd.lt.theta(10)) then    ! Winter
          multtf = 1.                !  turnover of foliage on
          multtl = 0.                !  turnover of labile C off
       endif
       if ( gdd .ge. theta(10) ) then
          if ( max_fol .eq. 1 ) then ! Spring
             multtl = 1.             !  turnover of foliage on
             multtf = 0.             !  turnover of labile C off
          else                       ! Summer
             multtl = 0.             !  turnover of foliage off
             multtf = 0.             !  turnover of labile C off
          endif
       endif
       if ( (A(8,m) .ge. theta(15) ) .or. ( time%day .ge. 200 ) ) then  ! theta15 = max Cf
          max_fol = 0.
          multtl  = 0.
       endif
       if ( ( time%day .ge. 200 ) .and. ( mint .lt. theta(11) ) ) then ! drop leaves - N hemisphere
          multtf = 1.
       endif

       ralabfrom = theta(13) * A(18,m) * theta(14) * multtl * rate * ts_length
       ralabto   = ( 1. - theta(12) ) * theta(5) * A(8,m) * theta(14) * multtf * rate * ts_length

    endif ! ends decid/evergreen if block.

    !-- Now update all the A's and the pools --

    if (user_opts%veg_is_deciduous) then
       ! autotrophic respiration..
       A(1,m) = A(21,m) * theta(16) * ts_length
       ! allocation to Car..
       Aar  = A(16,m) * theta(2) + ralabfrom + ralabto
    else
       ! autotrophic respiration..
       A(1,m) = A(17,m) * theta(10) * ts_length
       ! allocation to Car..
       Aar  = A(16,m) * theta(2)
    endif

    npp   = ( 1 - theta(2) ) * A(16,m)                           ! npp remaining

    if ( user_opts%veg_is_deciduous ) then
       A(2,m) = min( theta(15) - A(8,m), npp * theta(3) ) * multtl + A(20,m)
                                                               ! ...Af-deciduous
       npp  = npp - min( theta(15) - A(8,m) , npp * theta(3) ) * multtl
                                                               ! ...NPP remaining after leaf growth
    else
       A(2,m) = npp * theta(3)                                   ! Af-evergreen
       npp  = npp - A(2,m)                                       ! NPP remaining after leaf growth
    endif

    A(4,m)  = npp * theta(4)                                     ! Ar
    A(3,m)  = ( 1 - theta(4) ) * npp                             ! Aw
    if (user_opts%veg_is_deciduous) then
       if (phen_sim) A(5,m) = theta(5) * A(8,m) * ts_length * theta(12) * multtf ! Lf
    else
       if (phen_sim) A(5,m) = theta(5) * A(8,m) * ts_length      ! Lf
    endif
    A(6,m)  = theta(6) * A(9,m)  * ts_length                       ! Lw
    if (phen_sim) A(7,m)  = theta(7) * A(10,m) * ts_length         ! Lr
    A(11,m) = theta(8) * A(14,m) * ts_length * rate                ! Rh1
    A(12,m) = theta(9) * A(15,m) * ts_length * rate                ! Rh2
    A(13,m) = theta(1) * A(14,m) * ts_length * rate                ! D
    if (user_opts%veg_is_deciduous) then
       ! Atolab: 
       A(19,m) = ( 1. - theta(12) ) * theta(5) * A(8,m) * ( 1. - theta(14) ) &
                    * multtf * rate * ts_length
       ! Afromlab
       A(20,m) = theta(13) * A(18,m) * (1. - theta(14)) * multtl * rate * ts_length
    endif

    ! Calculate the pools
    if ( user_opts%veg_is_deciduous ) then
       A(8,m)  = A(8,m)  +  A(2,m) -  A(5,m) - A(19,m) - ralabto         ! Cf
       A(18,m) = A(18,m) + A(19,m) - A(20,m) - ralabfrom               ! Labile Pool
       A(21,m) = A(21,m) + Aar   - A(1,m)                            ! Car
    else
       if (phen_sim) then
         A(8,m) = A(8,m)  + A(2,m) - A(5,m)             ! Cf
       else
         A(8,m) = totlai * LMA
       endif
       A(17,m) = A(17,m) + Aar  - A(1,m)                             ! Car
    endif
    A(9,m)  = A(9,m)  +  A(3,m) -  A(6,m)                              ! Cw
    if (phen_sim) then
      A(10,m) = A(10,m) +  A(4,m) -  A(7,m)                    ! Cr
    else
      A(10,m) = A(8,m) * root_leaf_ratio
    endif
    A(14,m) = A(14,m) +  A(5,m) +  A(7,m) - A(11,m) - A(13,m)              ! Clit
    A(15,m) = A(15,m) + A(13,m) - A(12,m) +  A(6,m)                      ! Csom

    rootbiomass = A(10,m) * 2.     ! convert from gC to g biomass

  end subroutine predictor
  !
  !----------------------------------------------
  !
  subroutine roots ( A , m )

    ! determines root distribution given rootbiomass !

    use math_tools,         only: zbrent
    use scale_declarations, only: pi, core
    use soil_structure,     only: layer_depth, max_depth, rootbiomass, rootdens, rootl, rootlength, &
                                  rootmass, rootrad, rootreach, root_k, surfbiomass, thickness,     &
                                  root_distr_param, model_rootreach
    use veg,                only: !A
    use enkf

    implicit none

    ! arguments
    integer,intent(in) :: m ! ensemble member subscript
    real,dimension(ndim,nrens),intent(inout) :: A

    ! local variables..
    integer        :: i
    real           :: cumdepth, depth, mult, preb, rootdepth, rootxsecarea, slpa, x1, x2, xx(1), prev, curr
    real,parameter :: xacc = 0.0001
    real :: dummyR
    integer :: dummyI

    depth = 0. ; preb  = 0. ; prev = 0. ; curr = 0. ; rootmass = 0.

    ! convert from gC to g biomass
    rootbiomass = A(10,m) * 2.

    ! always provide a minimum root biomass
    rootbiomass = max( 5. , rootbiomass )

    rootxsecarea = pi * rootrad**2    ! root X-sectional area (m2)

    rootdepth = max_depth * rootbiomass / ( root_k + rootbiomass )    ! rmass=fine root C

    i = size( shape( layer_depth ) )  
    xx = minloc( layer_depth , MASK=layer_depth > rootdepth )
    if (model_rootreach == 1.) then
        rootl = int( xx(1) )
        rootreach = layer_depth( rootl )    ! to what soil layer do the roots pentrate?
    else
    rootl = int( rootreach * 10.)
    endif

    ! ensure 50% of root mass is in top 25% of rooted layers 
    ! see pipochronosequence.xls for derivation of this relationship
    !root_distr_param: SPA default = -0.006
    mult = min( 1. / thickness(1) , max( 2.0 , 11. * exp( -root_distr_param * rootbiomass ) ) ) ! original exponent: -0.006

    ! assume surface root density is related to total root biomass by a scalar dependent on root biomass
    surfbiomass = rootbiomass * mult

    if ( rootl .gt. 1 ) then
      x1 = 0.1
      x2 = 10.
      ! determine slope of root distribution given rooting depth 
      ! and ratio of root mass to surface root density
      slpa = zbrent( rootdist , x1 , x2 , xacc )
      prev = 1. / slpa
      cumdepth = 0.
      do i = 1 , rootl
        cumdepth      = cumdepth + thickness( i )
        curr          = 1 / slpa * exp( -slpa * cumdepth )
        ! rootmass is g biomass, i.e. twice the C content
        rootmass(i)   = ( ( prev - curr ) * surfbiomass ) / thickness( i ) !g root m-3 soil
        rootlength(i) = rootmass( i ) / ( rootdens * rootxsecarea )    ! (m root m-3 soil)
        prev          = curr
      enddo
    else
       rootmass( 1 ) = rootbiomass / thickness(i)
    endif

  end subroutine roots
  !
  !------------------------------
  !
  real function rootdist(xin)    ! parameters

    use soil_structure

    implicit none

    ! arguments..
    real, intent(in) :: xin       ! slope parameter

    ! local variables..
    real    :: one,two

    one=(1.-exp(-xin*rootl*thickness(1)))/xin
    two=rootbiomass/surfbiomass

    ! see pipo chronosequence.xls for this relationship
    rootdist=(1.-exp(-xin*rootreach))/xin-rootbiomass/surfbiomass

  end function rootdist
  !
  !------------------------------
end module allocate_carbon
