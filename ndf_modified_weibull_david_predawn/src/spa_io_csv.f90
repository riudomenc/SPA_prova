! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2011 SPA contributors.                   !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module spa_io_csv

  !! contains routines specific to reading/writing csv files. !!

  use scale_declarations, only: nos_canopy_layers, steps
  use enkf

  implicit none

  ! to keep track of whether write-routine already called this timestep..
  real,private :: last_call(nos_canopy_layers) = 0.

  ! Used in the write_csv_output s/r..
  real,dimension( ndim ),private :: sum_A
  real,dimension(steps),private :: lemod, mmmod, neemod, timeflux, wetle


  ! procedure private to this module..
  private :: open_file

contains
  !
  !-------------------------------------------------------------------
  !
  subroutine read_met_csv( filename , met )

    ! open the met driver (csv) file, and load the data from it !

    use log_tools
    use scale_declarations, only: met_drivers, user_opts

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    type(met_drivers)           :: met

    ! local variables..
    character(len=200) :: line
    integer            :: i, num_lines, status
    real               :: dummy

    ! open file
    call open_file( filename , 23 , readonly=.true. )

    ! get header out of the way..
    read(unit=23,fmt=*) line

    ! count the number of remaining lines in the file..
    status = 0    ;     num_lines = 0
    do
      read(unit=23,fmt=*,iostat=status) line
      if ( status .ne. 0. ) exit
      num_lines = num_lines + 1
    enddo

    ! make the variables the correct size..
    allocate( met%ambient_co2( num_lines ) , &
                      met%par( num_lines ) , &
                   met%precip( num_lines ) , &
                   met%sw_rad( num_lines ) , &
             met%sfc_pressure( num_lines ) , &
                     met%temp( num_lines ) , &
                      met%vpd( num_lines ) , &
                 met%wind_spd( num_lines )   )

    ! Go back to start of file..
    rewind(23)

    ! get header out of the way again..
    read(unit=23,fmt=*) line

    ! and now load all data..
    do i = 1 , num_lines
      read(23,fmt=*)dummy,met%temp(i),met%ambient_co2(i),met%wind_spd(i),&
                    met%sw_rad(i),met%vpd(i),met%par(i),met%precip(i)
    enddo

    ! Check sanity of temperature..
    if ( ( maxval(met%temp) .gt. 150. ) .or. ( minval(met%temp) .lt. -100 ) ) then
      write(message,*) "Temperature appears to be in units other than Celcius or Kelvin."
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )
    endif

    ! Surface pressure not provided..
    met%sfc_pressure = met%const_sfc_pressure

    ! In case user doesn't want to use co2 in the file..
    if ( .not. user_opts%co2_in_met_file ) met%ambient_co2 = met%const_ambient_co2

    ! In case user doesn't want to use the par in the file..
    if ( .not. user_opts%par_in_met_file )  met%par = 2.3 * met%sw_rad 

    ! Adjust zero windspeeds to ensure always some (small) turbulence..
    where ( met%wind_spd .lt. 0.2 ) met%wind_spd = 0.2

  end subroutine read_met_csv
  !
  !-------------------------------------------------------------------
  !
  subroutine read_soils_csv( filename , decid_flag )

    ! open and read the contents of the soils.csv file. !
 
    use hydrol,            only: claypc, sandpc
    use scale_declarations,only: core 
    use snow_info,         only: snowheight, snowweight
    use soil_structure,    only: iceprop, layer_depth, mineralfrac, organicfrac, resprate, &
                                 rootrad, thickness, waterfrac, soiltemp

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    logical,intent(in)          :: decid_flag

    ! local variables..
    character(len=200) :: header
    integer            :: i

    ! Open the soils parameters file..
    call open_file( filename , 27 , readonly=.true. )

    ! read soil properties
    read(unit=27,fmt=*)header
    read(unit=27,fmt=*)header,(thickness(i),  i=1,core)
    read(unit=27,fmt=*)header,(layer_depth(i),i=1,core)
    read(unit=27,fmt=*)header,(organicfrac(i),i=1,core)
    read(unit=27,fmt=*)header,(mineralfrac(i),i=1,core)
    read(unit=27,fmt=*)header,(waterfrac(i),  i=1,core)
    read(unit=27,fmt=*)header,(soiltemp(i),   i=1,core)
    read(unit=27,fmt=*)header,(iceprop(i),    i=1,core)
    if (decid_flag) then
      read(unit=27,fmt=*)header,sandpc(1)    ! convert to array if sand content varies in all layers
      sandpc=sandpc(1)                       ! set parameter constant for all soil layers
      read(unit=27,fmt=*)header,claypc(1)    ! convert to array if clay content varies in all layers
      claypc=claypc(1)                       ! set parameter constant for all soil layers
    else
      read(unit=27,fmt=*)header,(sandpc(i),i=1,core)
      read(unit=27,fmt=*)header,(claypc(i),i=1,core)
    endif
    read(unit=27,fmt=*)header,rootrad
    read(unit=27,fmt=*)header,snowweight
    read(unit=27,fmt=*)header,snowheight
    read(unit=27,fmt=*)header,resprate

  end subroutine read_soils_csv
  !
  !-------------------------------------------------------------------
  !
  subroutine read_veg_csv(filename,decid_flag)

    ! open and read the contents of the veg.csv file. !

    use log_tools
    use scale_declarations,only: nos_canopy_layers, pi
    use soil_structure,    only: max_depth,max_storage,root_k,rootresist,through_fall,root_distr_param,rootreach,model_rootreach
    use veg                ! 28 out of 41 vars in veg are read or given values in here...
                           ! only: A,avN,capac,canht,conductivity,dimen,gplant,iota,kappac,kappaj,&
                           !       lafrac,lat,layerht,LMA,minlwp,nfrac,towerht,theta

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename
    logical,intent(in)          :: decid_flag

    ! local variables..
    character(len=200) :: header,dummy
    integer            :: i

    ! Open the plant physiological parameters
    call open_file( filename , 25 , readonly=.true. )

    !read vegetation parameters
    read(unit=25,fmt=*)header,(lafrac(i),i=1,nos_canopy_layers)  ! LA fraction
    read(unit=25,fmt=*)header,(nfrac(i),i=1,nos_canopy_layers)   ! nitrogen fraction in each layer
    read(unit=25,fmt=*)header,(layerht(i),i=1,nos_canopy_layers) ! layer heights
    read(unit=25,fmt=*)header,avN                   ! average foliar N, gN m-2 leaf area
    read(unit=25,fmt=*)header,gplant                ! conductivity
    read(unit=25,fmt=*)header,minlwp                ! critical LWP
    read(unit=25,fmt=*)header,iota                  ! stomatal efficiency
    read(unit=25,fmt=*)header,capac                 ! leaf capacitance    
    read(unit=25,fmt=*)header,lat                   ! latitude
    read(unit=25,fmt=*)header,dummy
    read(unit=25,fmt=*)header,dimen                 ! leaf size
    read(unit=25,fmt=*)header,rootresist            ! root resistivity
    read(unit=25,fmt=*)header,towerht               ! height of measurement tower
    read(unit=25,fmt=*)header,conductivity          ! Does conductance vary with stem length? 0=NO, 1=YES
    read(unit=25,fmt=*)header,kappac
    read(unit=25,fmt=*)header,kappaj
    read(unit=25,fmt=*)header,LMA                   ! leaf mass per area
    read(unit=25,fmt=*)header,root_leaf_ratio
    read(unit=25,fmt=*)header,max_depth             ! max rooting depth
    read(unit=25,fmt=*)header,root_k                ! root mass for reaching 50% max depth
    read(unit=25,fmt=*)header,theta(1)              ! dc        decomposition rate
    read(unit=25,fmt=*)header,theta(2)              ! fa        fraction of gpp respired
    read(unit=25,fmt=*)header,theta(3)              ! nf        npp allocated to foliage
    read(unit=25,fmt=*)header,theta(4)              ! nrr       remaining npp allocated to fine roots
    read(unit=25,fmt=*)header,theta(5)              ! tf        turnover rate of foliage
    read(unit=25,fmt=*)header,theta(6)              ! tw        turnover rate of wood
    read(unit=25,fmt=*)header,theta(7)              ! tr        turnover rate of fine roots
    read(unit=25,fmt=*)header,theta(8)              ! ml        mineralisation rate of litter
    read(unit=25,fmt=*)header,theta(9)              ! ms        mineralisation rate of SOM/CWD
    if (decid_flag) then
      read(unit=25,fmt=*)header,theta(10)           ! gddthresh GDD threshhold
      read(unit=25,fmt=*)header,theta(11)           ! lfal      Minimum temperature threshold
      read(unit=25,fmt=*)header,theta(12)           ! fl        Fraction of leaf loss to litter
      read(unit=25,fmt=*)header,theta(13)           ! tl        turnover rate of labile pool
      read(unit=25,fmt=*)header,theta(14)           ! cr        respiratory cost of labile transfers
      read(unit=25,fmt=*)header,theta(15)           ! cfmax     max foliar carbon stock
      read(unit=25,fmt=*)header,theta(16)           ! tar       turnover rate of autotrophic respiration pool
    else
      read(unit=25,fmt=*)header,theta(10)           ! tar       turnover rate of autotrophic respiration pool
    endif
    read(unit=25,fmt=*)header, dummy !A(8)                  ! Cf
    read(unit=25,fmt=*)header, dummy !A(9)                  ! Cw
    read(unit=25,fmt=*)header, dummy !A(10)                 ! Cr
    read(unit=25,fmt=*)header, dummy !A(14)                 ! Clit
    read(unit=25,fmt=*)header, dummy !A(15)                 ! Csom
    if (decid_flag) then
      read(unit=25,fmt=*)header, dummy !A(18)               ! Clab
      read(unit=25,fmt=*)header, dummy !A(21)               ! Caresp
    else
      read(unit=25,fmt=*)header, dummy !A(17)               ! Caresp
    endif
    read(unit=25,fmt=*)header,through_fall
    read(unit=25,fmt=*)header,max_storage
    read(unit=25,fmt=*)header,root_distr_param
    read(unit=25,fmt=*)header,rootreach
    read(unit=25,fmt=*)header,model_rootreach

    canht = layerht(1)
    lat   = lat * pi / 180.    ! convert degrees to radians

  end subroutine read_veg_csv
  !
  !-------------------------------------------------------------------

  subroutine read_enkf_csv(filename,initmean,initerr)

  ! open and read the contents of the EnKF_setup.csv file. !
    use log_tools
    use enkf
    use veg

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename

    ! local variables
    real,intent(out) :: initmean(ndim), initerr(ndim)			! initial mean values of ensemble state vector
    real :: perturbation
    character(len=10) :: dummy
    integer :: i

    ! Open the EnKF setup file
    call open_file( filename , 26 , readonly=.true. )

    read(26,*)dummy,(initmean(i),i=1,ndim)
    read(26,*)dummy,(modvar(i),i=1,ndim)
    read(26,*)dummy,(initerr(i),i=1,ndim)
    read(26,*)dummy,(obsvar(i),i=1,maxobs)
    read(26,*)dummy,(abs_errors(i),i=1,ndim)
    read(26,*)dummy,(analyse(i),i=1,ndim)
    read(26,*)dummy,(lo(i),i=nsv+1,ndim)
    read(26,*)dummy,(hi(i),i=nsv+1,ndim)
    read(26,*)dummy,(SV_names(i),i=1,ndim)

  ! (OPTIONAL) scale input initerr from percentual to absolute for variables in SV
  ! initerr((nsv+1):ndim) = ( (hi - lo) * initerr((nsv+1):ndim) ) / b((nsv+1):ndim)

  end subroutine read_enkf_csv

  !-------------------------------------------------------------------
  subroutine read_obs_csv(filename)

  ! open and read the contents of the EnKF_setup.csv file. !
    use log_tools
    use enkf
    use obsdrivers,only : all_baseobs

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename

    ! local variables..
    character(len=200) :: line
    integer            :: i, j , num_lines, status
    real               :: dummy

    ! open file
    call open_file( filename , 28 , readonly=.true. )

    ! get header out of the way..
    read(unit=28,fmt=*) line

    ! count the number of remaining lines in the file..
    status = 0    ;     num_lines = 0
    do
      read(unit=28,fmt=*,iostat=status) line
      if ( status .ne. 0. ) exit
      num_lines = num_lines + 1
    enddo

    ! make the observation array the correct size
    if (Prades) then
        allocate(all_baseobs( 3 , num_lines ))
    else
        allocate(all_baseobs( maxobs , num_lines ))
    endif

    ! go back to start of file
    rewind(28)

    ! get header out of the way again..
    read(unit=28,fmt=*) line

    ! and now load all data..
    do i = 1 , num_lines
      if (Prades) then
        read(28,*)dummy,(all_baseobs(j,i),j=1,3)
      else
        read(28,*)dummy,(all_baseobs(j,i),j=1,maxobs)
      endif
    enddo

  end subroutine read_obs_csv

  !-------------------------------------------------------------------

  subroutine read_phen_csv(filename)

  ! open and read the contents of the EnKF_setup.csv file. !
    use log_tools
    use scale_declarations

    implicit none

    ! arguments..
    character(len=*),intent(in) :: filename

    ! local variables..
    character(len=200) :: line
    integer            :: i, j , num_lines, status
    real               :: dummy

    ! open file
    call open_file( filename , 29 , readonly=.true. )

    ! get header out of the way..
    read(unit=29,fmt=*) line

    ! count the number of remaining lines in the file..
    status = 0    ;     num_lines = 0
    do
      read(unit=29,fmt=*,iostat=status) line
      if ( status .ne. 0. ) exit
      num_lines = num_lines + 1
    enddo

    ! make the phenology array the correct size
    allocate(all_phenology( num_lines ))

    ! go back to start of file
    rewind(29)

    ! get header out of the way again..
    read(unit=29,fmt=*) line

    ! and now load all data..
    do i = 1 , num_lines
      read(29,*)dummy,all_phenology(i)
    enddo

  end subroutine read_phen_csv

  !-------------------------------------------------------------------

  subroutine open_output_csv( outdir , decid_flag )

    ! This opens the standard output files !
    !  and writes a header to each.        !

    use scale_declarations, only: fname_length, nos_canopy_layers
    use veg,                only: lafrac

    implicit none

    ! arguments..
    character(len=*),intent(in) :: outdir
    logical,intent(in)          :: decid_flag

    ! local variables..
    character(fname_length) :: filename
    integer                 :: i

    ! Open all output files..
    call open_file( trim(outdir)//'daily.csv', 30 , header = &
         'daymm,modle,modess,modwet,SWP,rplant,CSRes,lsc,totalflux' )

    if (decid_flag) then
      call open_file( trim(outdir)//'gdd.csv', 31 , header='time,gdd,mint' )
    endif
    call open_file( trim(outdir)//'budget.csv', 32 , header='delC,neesum' )
    call open_file( trim(outdir)//'hourly.csv', 35 , header = &
                      'gpp,ra,rh1,rh2,le,trans,ess,wetle,gpp,nee' )
    call open_file( trim(outdir)//'energy.csv', 40 , header = &
                     'Qh,Qe,Qn,Qc,airt,surfacet,soilt2,drythick' )
    do i = 1 , nos_canopy_layers
      if ( lafrac(i) .ne. 0. ) then 
        write(filename,"(A,i0,A)")"layer_",i,".csv"
        call open_file( trim(outdir)//trim(filename), 40+i , header = &
                 'LWP,stom_conduct,rsoil,rplant' )
               ! 'gs,agr,res,psil,ci,etr,tempdf,wdef,par,rad,cca,la,cistar,an' )
      endif
    enddo
    call open_file( trim(outdir)//'solar_part1.csv', 56 , header=&
                   'soilnet,soilabsn,soilabls,soilabsp,par_top' )
    call open_file( trim(outdir)//'solar_part2.csv', 57 , header=&
       'sum(abspar_sun),(1.-fdiff)*par_top,(parabs_sun(i),i=1,10)' )
    call open_file( trim(outdir)//'solar_part3.csv', 58 , header=&
     'sum(abspar_shade),fdiff*par_top,(parabs_shade(i),i=1,10),skyabsp' )
    call open_file( trim(outdir)//'solar_part4.csv', 59 , header=&
                           'fdiff,(leaf_sun(i),i=1,10),check' )
    call open_file( trim(outdir)//'soilwater.csv',   60 , header=&
     'w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w_swp' )
    call open_file( trim(outdir)//'upfrac.csv',      62 , header=&
     'up1,up2,up3,up4,up5,up6,up7,up8,up9,up10,up11,up12,up13,up14,up15' )
    call open_file( trim(outdir)//'iceprop.csv',     68 , header=&
            'i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15' )
    if (decid_flag) then
      ! Deciduous (21 SVs): Ra, Af, Aw, Ar, Lf, Lw, Lr, Cf, Cw, Cr, Rh1, Rh2, D, Clit, Csom, G, ET, Clab, Atolab, Afromlab, Caresp
      call open_file( trim(outdir)//'stocks.csv',    70 , header=&
        'A(8)=Cf,A(9)=Cw,A(10)=Cr,A(14)=Clit,A(15)=Csom,A(18)=Clab,A(21)=Caresp' )
      call open_file( trim(outdir)//'fluxes.csv',    71 , header=&
        'Ra,Af,Aw,Ar,Lf,Lw,Lr,Rh1,Rh2,D,G,Atolab,Afromlab,neesum,sumresp,daymm' )
    else
      ! Evergreen (17 SVs): Ra, Af, Aw, Ar, Lf, Lw, Lr, Cf, Cw, Cr, Rh1, Rh2, D, Clit, Csom, G, Caresp
      call open_file( trim(outdir)//'stocks.csv',    70 , header=&
        'A(8)=Cf,A(9)=Cw,A(10)=Cr,A(14)=Clit,A(15)=Csom,A(17)=Caresp')
      call open_file( trim(outdir)//'fluxes.csv',    71 , header=&
        'Ra,Af,Aw,Ar,Lf,Lw,Lr,Rh1,Rh2,D,G,neesum,sumresp,daymm')
    endif
    call open_file( trim(outdir)//'soiltemp.csv',    83 , header=      &
       'soil-temp(layer 1),soil-temp(layer 2),soil-temp(layer 3),'     &
     //'soil-temp(layer 4),soil-temp(layer 5),soil-temp(layer 6),'     &
     //'soil-temp(layer 7),soil-temp(layer 8),soil-temp(layer 9),'     &
     //'soil-temp(layer 10),soil-temp(layer 11),soil-temp(layer 12),'  &
     //'soil-temp(layer 13),soil-temp(layer 14),soil-temp(layer 15)'   )

    call open_file( trim(outdir)//'waterfluxes.csv', 94 , header=  &
                   'modess,runoff,ppt,trans,delw,disc,cwtr,check,' &
                 //'diff,modwet,unint,canst,snow_watermm'          )

    call open_file( trim(outdir)//'soilstatus.csv', 101 , header=      &
                   'of,totet,uf,sin,wtr,wch,flux,chk,w+1,w-1,delw,'    &
                 //'chkdff,ppt+,w+,w-,w+c,w-c,snow weight,snow height' )

    ! EnKF output files
    call open_file( trim(outdir)//'B_sd.csv', 113 )
    call open_file( trim(outdir)//'spurious.csv', 114 , header=       &
                   'mean_stddev,inflate')
    call open_file( trim(outdir)//'inflate.csv', 115 )
    call open_file( trim(outdir)//'parcor.csv', 116 , header =        &
                   'iota_gpl,iota_rootr,iota_lwp,iota_root_leaf_ratio,iota_rbelow,' &
                  // 'gpl_rootr,gpl_lwp,gpl_root_leaf_ratio,gpl_rbelow,'       &
                  // 'rootr_lwp,rootr_root_leaf_ratio,rootr_rbelow,'           &
                  // 'lwp_root_leaf_ratio,lwp_rbelow,root_leaf_ratio_rbelow')

    call open_file( trim(outdir)//'SWC.csv', 117)
    call open_file( trim(outdir)//'transparcor.csv', 118 , header =        &
                   'iota,gplant,rootresist,minlwp,root_leaf_ratio,rbelow')
    call open_file( trim(outdir)//'capacitanceparcor.csv', 119 , header =        &
                   'iota,gplant,rootresist,minlwp,lappac,rbelow')
    call open_file( trim(outdir)//'soilR.csv', 120 , header =        &
                   'soilR,soilR1,soilR2,CSR(1),CSR(2),CSR(3),CSR(4),CSR(5)')
    if (decid_flag) then
      call open_file( trim(outdir)//'filter.csv', 111 , header=          &
                   'Ra,Af,Aw,Ar,Lf,Lw,Lr,Cf,Cw,Cr,Rh1,Rh2,D,Clit,'  &
                 //'Csom,G,ET,Clab,Atolab,Afromlab,Caresp,trans,capacitance,iota,gplant,capacitance,root_leaf_ratio,,rbelow,NEE' )
      call open_file( trim(outdir)//'filter_sd.csv', 112 , header=          &
                   'Ra,Af,Aw,Ar,Lf,Lw,Lr,Cf,Cw,Cr,Rh1,Rh2,D,Clit,'  &
                 //'Csom,G,ET,Clab,Atolab,Afromlab,Caresp,trans,capacitance,iota,gplant,capacitance,root_leaf_ratio,,rbelow,NEE' )
    else
      call open_file( trim(outdir)//'filter.csv', 111 , header=          &
                   'Ra,Af,Aw,Ar,Lf,Lw,Lr,Cf,Cw,Cr,Rh1,Rh2,D,Clit,'  &
                 //'Csom,G,Caresp,,,,,trans,capacitance,iota,gplant,capacitance,root_leaf_ratio,,rbelow,NEE' )
      call open_file( trim(outdir)//'filter_sd.csv', 112 , header=          &
                   'Ra,Af,Aw,Ar,Lf,Lw,Lr,Cf,Cw,Cr,Rh1,Rh2,D,Clit,'  &
                 //'Csom,G,Caresp,,,,,trans,capacitance,iota,gplant,capacitance,root_leaf_ratio,,rbelow,NEE' )
    endif

    ! open dummy output file for writing variables to be analysed
    call open_file( trim(outdir)//'dummy.csv', 999 , header='V1,V2,V3,V4,V5,V6,V7,V8,V9,V10' )

  end subroutine open_output_csv
  !
  !-------------------------------------------------------------------
  !
  subroutine write_output_csv( time , prevwater , decid_flag , A )

    ! write to output files. HARD-CODED. !

    use clim,               only: daypar, dayppt, temp_bot, temp_top, wetev
    use hourscale,          only: canopy_store, discharge, runoff, unintercepted
    use irradiance_sunshade,only: check
    use scale_declarations, only: met, steps, time_holder
    use snow_info,          only: snow_watermm
    use soil_structure,     only: watericemm, weighted_swp
    use veg,                only: canopy_soil_resistance, conductivity, ess, flux, gplant, &
                                  gppt, lai, layerht, prevC, soiletmm, totevap, transt !, A
    use enkf

    implicit none

    ! arguments..
    type(time_holder),intent(in) :: time
    real,intent(inout)           :: prevwater
    logical,intent(in)           :: decid_flag
    real,dimension(ndim,nrens)   :: A

    ! local variables..

    integer               :: i, j
    real                  :: check2, countle, countnee, currentC, currentwater, daymm, dayrad, &
                             daytrans, delC, delwater, gppsum, lambda, lsc, modess, modess2,   &
                             modess3, modle, modwet, modwetle, nightflux, rplant, sumresp,     &
                             timestep, totalflux, totle, totmodle, totmodnee, totnee, neesum
    real,dimension(steps) :: posess
    real,dimension(ndim)  :: Amean

    Amean = sum(A, dim=2)/real(nrens)

    ! calculations...
    ! Deciduous (21 SVs): Ra, Af, Aw, Ar, Lf, Lw, Lr, Cf, Cw, Cr, Rh1, Rh2, D, Clit, Csom, G, ET, Clab, Atolab, Afromlab, Caresp
    ! Evergreen (17 SVs): Ra, Af, Aw, Ar, Lf, Lw, Lr, Cf, Cw, Cr, Rh1, Rh2, D, Clit, Csom, G, Caresp
    if ( time%step .eq. 1 ) then
      ! reset sums
      countle  = 0.  ;  countnee  = 0.   ;  lemod    = 0.  ;  mmmod    = 0.
      neemod   = 0.  ;  nightflux = 3.0  ;  sum_A    = 0.  ;  timeflux = 0.
      totle    = 0.  ;  totmodnee =  0.  ;  totmodle = 0.  ;  totnee=0.
      wetle    = 0.
    endif

    lambda              = 1000. * ( 2501.0 - 2.364 * temp_top )                  ! latent heat of vapourisation (J kg-1)
    wetle( time%step )  = lambda * wetev( time%step ) / time%seconds_per_step    ! wet leaves evap convert (mm t-1 to Wm-2)
    neemod( time%step ) = (Amean(1) + Amean(11) + Amean(12) - Amean(16)) / time%seconds_per_step * 1e6 / 12. ! estimate NEE (umol m-2 s-1)
    ! sum modelled canopy & modelled soil LE, and wet canopy evap..
    lemod( time%step )  = transt( time%step ) + ess( time%step ) + wetle( time%step )
    ! lemod x length of timestep..
    mmmod( time%step )  = lemod( time%step) * ( time%seconds_per_step / lambda ) 

    ! sum flux from each tree layer..
    do j = 1 , 10
       timeflux( time%step ) = timeflux( time%step ) + flux( time%step , j )   ! (mmol m-2 GA s-1)
    enddo

    call write_to_file( 35 , (/ Amean(16), Amean(1), Amean(11), Amean(12), lemod(time%step), transt(time%step), &
                        ess(time%step), wetle(time%step), gppt(time%step), neemod(time%step) /) )

    sum_A = sum_A + Amean    !sum SV

    if ( time%step .eq. steps ) then
      gppsum    = sum(gppt) * time%seconds_per_step * 1e-6 * 12. ! total GPP (gC m-2 d-1)
      modle     = sum(lemod) * time%seconds_per_step * 1e-6      ! modelled ecosystem water loss to atmosphere (MJ m-2 d-1)
      daytrans  = sum(transt) * time%seconds_per_step * 1e-6     ! modelled canopy transpiration (MJ m-2 d-1)
      modess    = sum(soiletmm)                                  ! modelled soil evaporation (mm d-1)
      modess2   = sum(ess) * time%seconds_per_step * 1e-6        ! modelled soil evaporation (MJ m-2 d-1)
      posess    = max( 0. , ess )                                ! remove dew component
      modess3   = sum(posess) * time%seconds_per_step * 1e-6     ! modelled soil evaporation (MJ m-2 d-1)
      dayrad    = daypar * 1e-6 / 2.3                            ! measured irradiance (MJ m-2 d-1)
      modwet    = sum(wetev)                                     ! modelled evap from wetted canopy (mm d-1)
      modwetle  = sum(wetle) * time%seconds_per_step * 1e-6      ! modelled evap from wet canopy (MJ m-2 d-1)
      daymm     = sum(mmmod)                                     ! total ET (mm d-1)
      totalflux = sum(timeflux) * time%seconds_per_step * 0.001 * 18 * 0.001  ! total flux  (mm d-1 or kg m-2 ga d-1)
      sumresp   = sum_A(1) + sum_A(11) + sum_A(12)
      neesum    = sumresp - sum_A(16)

      ! determine leaf specific conductance for 2nd canopy layer
      if ( lai(2) .gt. 0 ) then
        if ( conductivity .eq. 1 ) then
          rplant = layerht(2) / ( Amean(25) * lai(2) ) ! MPa s mmol-1 layer
        else
          rplant = 1. / ( Amean(25) * lai(2) )         ! Conductance is constant with height
        endif
        lsc = ( 1. / ( rplant + canopy_soil_resistance(2) ) ) / lai(2)
      else
        rplant = -999.
        lsc    = -999.
      endif

      call write_to_file( 30, (/daymm,modle,modess,modwet,weighted_SWP,rplant,&
                                 canopy_soil_resistance(2),lsc,totalflux/)    )
      if (decid_flag) then
        ! deciduous stocks..
        call write_to_file( 70, (/Amean(8:10),Amean(14:15),Amean(18),Amean(21)/) )
        ! fluxes, daily sums..
        call write_to_file( 71, (/sum_A(1:7),sum_A(11:13),sum_A(16), &
                                  sum_A(19:20),neesum,sumresp,daymm/) )
        currentC     = Amean(8) + Amean(9) + Amean(10) + Amean(14) + Amean(15) + Amean(18) + Amean(21)
      else
        ! evergreen stocks..
        call write_to_file( 70, (/Amean(8:10),Amean(14:15),Amean(17)/) )
        ! fluxes, daily sums..
        call write_to_file( 71, (/sum_A(1:7),sum_A(11:13),sum_A(16),neesum,sumresp,daymm/) )                     
        currentC     = Amean(8) + Amean(9) + Amean(10) + Amean(14) + Amean(15) + Amean(17)
      endif

      currentwater = sum( watericemm )
      delwater     = currentwater - prevwater    ! change in soil water storage (mm)
      check2       = -1. * ( modess + runoff * 1e3 + totevap + discharge - unintercepted ) ! total water fluxes
      !check and delwater should be the same magnitude

      call write_to_file( 94, (/ modess, runoff*1e3, dayppt, totevap,   &
            delwater, discharge, currentwater, check2, check-delwater , &
            modwet, unintercepted, canopy_store, snow_watermm /)        )
            ! all output in mm d-1 (canopys_store=mm)

      delC      = currentC - prevC
      prevwater = currentwater
      prevC     = currentC

      call write_to_file( 32, (/delC, neesum/) )

    endif

  end subroutine write_output_csv

  !-------------------------------------------------------------------

  subroutine write_enkf_output( A , B , SWC , LWPout , LWPout_sd , runoff_out, &
    BGtot, BG1, BG2 , CSR , rsoilout , rplantout )

    use enkf
    use math_tools
    use scale_declarations !,only: time
    use veg, only : lai , stom_conduct

    implicit none

    ! arguments
    real,dimension(ndim , nrens),intent(in)      :: A , B
    real,dimension(core),intent(in)              :: SWC, BGtot, BG1, BG2
    real,dimension(nos_canopy_layers),intent(in) :: LWPout,LWPout_sd,CSR,rsoilout,rplantout
    real,intent(in) :: runoff_out

    ! local variables
    real      :: Aplus(ndim + 1 , nrens) ! Ensemble matrix with extra space for NEE
    real      :: AT(nrens , ndim + 1)    ! Transpose of Aplus, for stats
    real,save :: AT_cum(nrens , ndim + 1)  ! matrix holding cumulative A values
    real      :: BT(nrens , ndim)
    real      :: spurious_out(2)
    real      :: stat(15,ndim+1)			! Matrix holding stats on analysis
    integer   :: i , k , j , tm

    tm = time%steps_count - (( time%year -1 ) * user_opts%days_in_year * user_opts%timesteps_per_day)
    j = tm

    Aplus = 0.

    do i = 1 , ndim
      do k = 1 , nrens
        if ( i .le. nsv ) then
          Aplus( i , k ) = A( i , k )
        else !transform parameters
          Aplus( i , k ) = lo( i ) + ( hi( i ) - lo( i ) ) * ( sin( A( i , k ) ) ) ** 2
        endif
      enddo
    enddo

    do k = 1 , nrens
      Aplus( ndim + 1 , k ) = A( 1 , k ) + A( 11 , k ) + A( 12 , k ) - A( 16 , k ) ! NEE
    enddo

    ! transpose A for stats
    AT = transpose( Aplus )
    BT = transpose( B     )

    ! build cumulative A
    AT_cum = AT_cum + AT

    ! do stats on A
    do i = 1 , ndim + 1 ! loop over dimensions of A
      analysismean( j , i ) = mean( AT( : , i ) )
      stat( 1 , i ) = analysismean( j , i )
      stat( 2 , i ) = stddev( AT( : , i ) )
      analysismean_cum( j , i ) = mean( AT_cum( : , i ) )
      stat( 3 , i ) = analysismean_cum( j , i )
      stat( 4 , i ) = stddev( AT_cum( : , i ) )

    enddo

    ! do stats on B
    do i = 1 , ndim
      stat( 5 , i ) = mean( BT( : , i ) )
      stat( 6 , i ) = stddev( BT( : , i ) )
    enddo
    spurious_out(1) = spurious_mean; spurious_out(2) = spurious

    call write_to_file( 111, stat( 1 , : ) )
    call write_to_file( 112, stat( 2 , : ) )
    call write_to_file( 113, stat( 6 , : ) )
    call write_to_file( 114, (/spurious_out, runoff_out/) )
    call write_to_file( 115, inflate_exe)
    call write_to_file( 116, parcor )
    call write_to_file( 117, SWC )
    call write_to_file( 118, transparcor )
    call write_to_file( 119, capacitanceparcor )
    call write_to_file( 120, (/sum(BGtot),sum(0.2/BG1),sum(0.2/BG2),CSR(1:5)/) )
    !call write_to_file( 121, rplantout(1:5) )
    do i = 1 , nos_canopy_layers
      if (lai(i).gt.0) then
        call write_to_file( 40+i , (/LWPout(i),stom_conduct(i),rsoilout(i),rplantout(i)/) )
      endif
    enddo

  end subroutine write_enkf_output
  !
  !-------------------------------------------------------------------
  !
  subroutine close_output_csv

     use scale_declarations, only: nos_canopy_layers
     use veg,                only: lafrac

     implicit none

     ! local variables..
     integer :: i
 
     close(unit=30,status='keep')
     close(unit=31,status='keep')
     close(unit=32,status='keep')
     close(unit=35,status='keep')
     close(unit=40,status='keep')
     do i=1,nos_canopy_layers
       if (lafrac(i) .ne. 0.) close(unit=40+i,status='keep')
     enddo
     close(unit=56,status='keep')
     close(unit=57,status='keep')
     close(unit=58,status='keep')
     close(unit=59,status='keep')
     close(unit=60,status='keep')
     close(unit=62,status='keep')
     close(unit=68,status='keep')
     close(unit=70,status='keep')
     close(unit=71,status='keep')
     close(unit=83,status='keep')
     close(unit=94,status='keep')
     close(unit=101,status='keep')
     close(unit=111,status='keep')
     close(unit=112,status='keep')
     close(unit=113,status='keep')

  end subroutine close_output_csv
  !
  !-------------------------------------------------------------------
  !
  subroutine write_assimilate_output_csv( time, output_data )

    ! moved writes from assimilate (leaf.f90) to here !

    use clim,               only: coa, gdd, mint
    use log_tools
    use metab,              only: an, ci
    use meteo,              only: gbb, la, par, psil, rad, temp, wdef
    use scale_declarations, only: time_holder

    implicit none

    ! arguments..
    type(time_holder),intent(in)    :: time
    real,             intent(inout) :: output_data(:)

    ! local variables..
    integer                 :: output_layer

    output_layer = int( output_data(1) )   ! convert back from real to integer
    output_layer = max( output_layer , 1 ) ! make sure index is at least 1

    ! Check we have not been called already this time step!
    if ( last_call(output_layer) .lt. time%daytime ) then

      call write_to_file( 40+output_layer , (/ &
!      write(unit=40+output_layer,fmt='(f0.2,14(",",f0.7))')         &
                         output_data(2:4), psil, ci,  &
                         output_data(5), output_data(6)-temp, wdef, &
                         par, rad, coa, la, an, gbb*output_data(7) /) )
 
      ! reset the carrier..
      output_data = 0.

      last_call(output_layer) = time%daytime

    else

      write(message,*)'write_assimilate_output_csv has already been called this timestep: ignoring additional call'
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )

    endif

  end subroutine write_assimilate_output_csv
  !
  !-------------------------------------------------------------------
  !
  subroutine write_soilday_output_csv( flag , output_data )

    ! moved writes from soilday (soil_functions.f90) to here !

    use clim,               only: gdd, mint, temp_bot, temp_top
    use hourscale,          only: freeze, hourtemp, overflow, Qc, Qe, Qh, Qn, surface_watermm, &
                                  totet, underflow
    use log_tools
    use scale_declarations, only: core
    use snow_info,          only: snowheight, snowweight
    use soil_structure,     only: drythick, fraction_uptake, iceprop, pptgain, soiltemp, waterfrac, &
                                  watergain, watericemm, waterloss, weighted_SWP

    implicit none

    ! arguments..
    integer,intent(in)    :: flag
    real,   intent(inout) :: output_data(:)

    ! local variables..
    integer :: i

    select case (flag)
    case (1)
      call write_to_file(40,(/ Qh, Qe, Qn, Qc, hourtemp-freeze, &
           output_data(1)-freeze, soiltemp(2)-freeze, drythick*1e3/) )

      call write_to_file(101, (/ 1e3*overflow, 1e3*totet,                      &
           1e3*underflow, surface_watermm, sum(watericemm), output_data(2:4),  &
           1e3*watergain(1), 1e3*waterloss(1), output_data(5:6), sum(pptgain), &
           sum(watergain), sum(waterloss) , watergain(core), waterloss(core),  &
           snowweight, snowheight /) )
      output_data = 0.

      call write_to_file( 60, (/ waterfrac(1:15), weighted_swp /) )
      call write_to_file( 68, iceprop(1:15) )
      call write_to_file( 83, (/soiltemp(1:15)-273.15, temp_bot/) )
    case (2)
      call write_to_file( 62, fraction_uptake(1:15) )
    case default
      write(message,*)"flag supplied to handle_output:",flag," was not recognised"
      call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
    end select

  end subroutine write_soilday_output_csv
  !
  !-------------------------------------------------------------------
  !
  subroutine write_solar_output_csv( output )

    ! moved writes from solar (light.f90) to here !

    use clim,               only: par_top

    implicit none

    ! arguments..
    real,          intent(inout) :: output(:)

    ! local variables..
    real,dimension(nos_canopy_layers) :: abspar_shade_out, abspar_sun_out, &
                               leaf_sun_out, parabs_shade_out, parabs_sun_out
    integer :: i

    ! moved solar output here..

    ! output_data1 =  (/ soilnet, soilabsn, soilabsl, soilabsp, par_top, fdiff, check, skyabsp,  &
    !                      abspar_sun(nos_canopy_layers),  &
    !                      parabs_sun(nos_canopy_layers),  &
    !                    abspar_shade(nos_canopy_layers),  &
    !                    parabs_shade(nos_canopy_layers),  &
    !                        leaf_sun(nos_canopy_layers)   /)

    abspar_sun_out   = output( 9                     : 9+nos_canopy_layers-1   )
    parabs_sun_out   = output( 9+nos_canopy_layers   : 9+2*nos_canopy_layers-1 )
    abspar_shade_out = output( 9+2*nos_canopy_layers : 9+3*nos_canopy_layers-1 )
    parabs_shade_out = output( 9+3*nos_canopy_layers : 9+4*nos_canopy_layers-1 )
    leaf_sun_out     = output( 9+4*nos_canopy_layers : 9+5*nos_canopy_layers-1 )

    ! WRITE: soilnet, soilabsn, soilabsl, soilabsp, par_top
    call write_to_file( 56 , output(1:5) )

    ! WRITE: sum(abspar_sun), (1.-fdiff)*par_top, (parabs_sun(i),i=1,10)
    call write_to_file( 57 , (/ sum(abspar_sun_out), (1.-output(6))*output(5), &
                                   (parabs_sun_out(i),i=1,nos_canopy_layers) /) )

    ! WRITE: sum(abspar_shade), fdiff*par_top, (parabs_shade(i),i=1,10), skyabsp
    call write_to_file( 58 , (/ sum(abspar_shade_out), output(6)*output(5) , &
                    (parabs_shade_out(i),i=1,nos_canopy_layers) , output(8) /) )

    ! WRITE: fdiff, (leaf_sun(i),i=1,10), check
    call write_to_file( 59 , (/ output(6), (leaf_sun_out(i), &
                           i=1,nos_canopy_layers), output(7) /) )

  end subroutine write_solar_output_csv
  !
  !-------------------------------------------------------------------
  !
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
  !-------------------------------------------------------------------
  !
  subroutine write_to_file( unit_nos , out_data )

    use log_tools
    use scale_declarations, only: time

    implicit none

    ! arguments..
    integer,intent(in) :: unit_nos
    real,   intent(in) :: out_data(:)

    ! local variables..
    integer           :: i, n
    logical           :: file_open
    character(len=7)  :: write_permit
    character(len=20) :: write_fmt

    ! check file is open and can be written to..
    inquire( unit=unit_nos, opened=file_open, write=write_permit )

    if ( .not. file_open ) then

      ! Unit-nos is not for a file that is open!!
      write(message,*)"Unit_nos:",unit_nos," is not connected to an open file!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )

    else if ( trim(write_permit) .ne. 'YES' ) then

      ! Unit-nos is for a file in read-only setting!!
      write(message,*)"Unit number:",unit_nos," is not open in write-mode!"
      call write_log( message , msg_fatal , __FILE__ , __LINE__ )

    else

      ! Okay to do the write..
      n = size( out_data )
      write(write_fmt,'(a,i0,a)') "(f0.2," , n ,'(",",f0.7))'
      write( unit=unit_nos , fmt=trim(write_fmt) ) &
                        time%daytime , ( out_data(i),i=1,n )
    endif

  end subroutine write_to_file
  !
  !-------------------------------------------------------------------
  !
end module spa_io_csv
!
!-------------------------------------------------------------------
!
