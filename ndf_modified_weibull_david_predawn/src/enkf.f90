module enkf

  implicit none

  logical            :: Vallcebre = .false., Prades = .true.
  logical            :: defoliated = .false., nondefoliated = .true.
  logical            :: measured_sapflow = .false. ! set to true if ONLY measured sapflow is affecting soil water balance
  integer, parameter :: ndim=29			! Dimension of model state
  integer, parameter :: nsv=23				! number of non-parameter state variables
  integer, parameter :: nrens=1  ! Number of ensemble members
  integer, parameter :: maxobs=1	        ! Maximum number of observations
  integer, parameter :: maxerr=1	        ! Maximum number of obseration errors
  integer, parameter :: totaltime=96*365   ! EC potential half hours, as of 22/8/07
  integer            :: switch						! 0=fitting parameters, 1=final filter
  integer            :: smooth = 0						! 0=no smooth, 1=smoothing
  integer            :: outon							! 0=low output, 1=high output
  logical            :: forward       = .false.	! 0=assimilation, 1=forward run
  logical            :: obserr_abs    = .true.  ! errors on obs absolute or percentual?
  logical            :: spurious_flag = .false. ! account for spurious correlation effects?
  logical            :: modvarswitch  = .true. ! add model uncertainty?
  logical            :: below_frac    = .false. ! belowground resistance a (fixed) fraction of total (i.e. below- + aboveground) resistance?
  logical            :: inflate_flag  = .false. ! inflate model parameters?
  logical            :: parallel_flag = .false. ! code parallelization activated?
  logical            :: read_initial_distribution = .false. ! read initial parameter distribution from file?
  real               :: forecastmean(totaltime,ndim+1)
  real               :: analysismean(totaltime,ndim+1),analysismean_cum(totaltime,ndim+1)
  real               :: hi(nsv+1:ndim), lo(nsv+1:ndim) , spurious_all(ndim) , spurious_mean
  real               :: obsvar(maxobs), modvar(ndim), abs_errors(ndim), analyse(ndim),countobs(maxobs), avobs(maxobs), cinit(nrens)
  real               :: inflate ! ensemble spread inflation factor
  real               :: spurious ! correction factor for spurious correlations
  real               :: parcor((ndim-nsv-1)*(ndim-nsv)/2) ! array to store parameter correlations
  real               :: transparcor(6),capacitanceparcor(6) ! array to store correlations between parameters and trans/capacitance
  real,dimension(ndim-nsv) :: minerr
  real,dimension(3)  :: inflate_exe = 0.
  real               :: sap_df, sap_ndf, sap_oak
  character(len=10)  :: SV_names(ndim)

contains
!------------------------------------------

subroutine enkf_calcs ( A , B , time , hbase )
! executes all EnKF related routines after model forward calculations (i.e. after the forecast)

  use obsdrivers, only : baseobs, sapobs
  use scale_declarations,only: time_holder
  use math_tools
  use clim

  implicit none

  ! arguments
  real,dimension(ndim,nrens),intent(inout) :: A , B
  real,dimension(ndim,maxobs),intent(in)   :: hbase
  type(time_holder),intent(inout) :: time

  ! local variables
  real                        :: sigma , rho , exist , ss , Aold(ndim , nrens ) , av_old( ndim ), astd
  real                        :: X5( nrens , nrens )  ! X5
  real                        :: covmat(ndim,ndim), cormat(ndim,ndim)	!covariance and correlation matrices
  real,dimension(ndim)        :: sigma2 , av
  real                        :: func				! sum-of-squares errors
  !real,dimension(maxobs)      :: baseobs              ! raw input data, including missing values
  real,dimension(maxobs)      :: measvar
  real,dimension(nrens)       :: R
  !real,dimension(ndim,maxobs) :: Hbase	! baseline observation operator
  real,parameter              :: dt = 1.0
  real,parameter              :: alpha = 0.	!time correlation of model errors
  real,allocatable,dimension (:,:) :: D 	! Matrix holding innovations
  real,allocatable,dimension (:,:) :: H		!observation operator
  real,allocatable,dimension (:,:) :: HA	! Matrix holding HA
  real,allocatable,dimension (:,:) :: S 	! Matrix holding HA transposed (?)
  real,allocatable,dimension (:,:) :: E		! Matrix holding observation perturbations
  real,allocatable,dimension (:)   :: avHA,avD,avE,avS,actobs
  integer                           :: i , k , n , nrobs , l , analysis_AB
  integer,dimension(4),save         :: iseed
  integer,allocatable,dimension (:) :: obloc
  integer :: n_analysis = 1
  character(len=20) :: write_fmt
  real,save :: modvarflag = 0.

  if (modvarswitch) modvarflag = 1.

  rho = sqrt( 1. / dt * ( 1. - alpha ) ** 2 / ( 1. / dt - 2. * alpha - 1. / dt * alpha ** 2 + 2 * alpha ** ( 1. / dt + 1. ) ) )
  ss = 0 ; n = 0

  av = sum( A , dim = 2 ) / real( nrens )

  !call errors( sigma2 , measvar ,av )

  if ( time%steps_count .eq. 1) iseed = (/213,876,1962,3755/)

  !Read the ensemble forecast into A (matrix holding ensemble members)
  do i = 1 , ndim
    !call slarnv( 3, iseed, nrens, R ) !random numbers from a normal distribution
    do k = 1 , nrens
      sigma = sqrt( sigma2( i ) )
      A( i , k ) = A( i , k ) + sqrt( dt ) * sigma * rho * R( k ) * modvarflag
    enddo
    R = R - mean( R ) ! remove non-zero mean from R, and scale stddev to 1 for creating B matrix
    R = R / stddev( R )
    do k = 1 , nrens
      B( i , k ) = R( k ) ! B is a matrix with zero mean and unit variance in each row,
        ! used for analysing spurious correlations and estimating size of inflation parameter
    enddo
  enddo

  if (Prades) then
    sap_df = sapobs(1); sap_ndf = sapobs(2); sap_oak = sapobs(3)
    if (defoliated)    baseobs = sapobs(1) ! select df  entry if df
    if (nondefoliated) baseobs = sapobs(2) ! select ndf entry if ndf
  endif

  exist = maxval( baseobs )  !are there any obs?
  if (exist.le.0.) exist = -999.
  if (forward) exist = -999.
  if (time%day.le.60) exist = -999.
  !if (exist < 0.005) exist = -999. ! don't assimilate potentially biased synthetic data
  !if (par_top .lt. 1. ) exist = -999. ! no nighttime assimilation

! If there are observations for this timestep....
  if ( exist .gt. -999. ) then

    !call varobs( maxobs , nrobs ) ! determine number of observations this timestep

    if (spurious_flag) n_analysis = 2

    do analysis_AB = 1 , n_analysis ! cycle over analysis with A and B

    if (analysis_AB .EQ. 2) A = B

    allocate (D(nrobs,nrens))		! Matrix holding innovations
    allocate (HA (nrobs,nrens))		! Matrix holding HA
    allocate (S (nrobs,nrens))		! Matrix holding HA�
    allocate (E(nrobs,nrens))		! Matrix holding observation perturbations
    allocate (H(nrobs,ndim))		! observation operator
    allocate (AvHa(nrobs),avD(nrobs),avE(nrobs),avS(nrobs))	!averages
    allocate (actobs(nrobs))		! the actual observations (missing data stripped out)
    allocate (obloc(nrobs))			! the locations of the actual observations in the full obs vector

    !call masks( nrobs , actobs , H , Hbase , obloc )

    !Compute the matrix HA, where H = measurement operator relating true model state to the observations d.
    HA = matmul( H , A )					! conversion from A to HA
    avHA = sum( HA , dim = 2 ) / real( nrens )	! determine expected value of HA

    ! Compute the measurement perturbations E
    do i = 1 , nrobs
      !call slarnv( 3, iseed, nrens, R ) !random numbers from a normal distribution
      R = R - mean( R ) ! subtract nonzero mean so that random observation errors have mean equal 0 and thus are not biased
      do k = 1, nrens
        E(i,k) = R(k) * sqrt( measvar( obloc( i ) ) )
      enddo
    enddo
    avE = sum( E , dim = 2 ) / real( nrens )

    ! Compute the innovations, D', where D' = D - HA.
    do i = 1 , nrens
      do k = 1 , nrobs
        D( k , i ) = E ( k , i ) + actobs( k ) - HA( k , i )
        ss = ss + ( D( k , i ) / ( countobs( obloc( k ) ) * avobs( obloc( k ) ) ) ) ** 2
        n = n + 1
      enddo
    enddo
    avD = sum( D , dim = 2 ) / real( nrens )

    ! Subtract expected values of HA from HA to get HA' (requires H to be linear)
    do i = 1 , nrens
      do k = 1 , nrobs
        S( k , i ) = HA( k , i ) - avHA( k )
      enddo
    enddo
    avS = sum( S , dim = 2 ) / real( nrens )
    av = sum( A , dim = 2 ) / real( nrens )
    av_old = av
    if (analysis_AB .EQ. 1) Aold = A

    !call analysis( A , D , E , S , nrobs , x5 )  ! Ensemble Kalman Filter

    ! restore ensemble values of non-analysis state variables
    if (analysis_AB .EQ. 1) then
      do i = 1 , ndim
        if (analyse(i).EQ.0) A( i , : ) = Aold( i , : ) ! restore pre-analysis value if state variable should not be updated by analysis
      enddo
      Aold = A
    endif

    if ( smooth .eq. 1 ) then
      !call smoother( X5 , timeid - 1 )
    endif
    deallocate (D, HA, S, E, H, AvHa , avD , avE , avS , actobs , obloc )

    !if (analysis_AB == 1) call covariancemat( A , covmat , cormat )

    func = ss / n * 1e6

    if (analysis_AB .EQ. 2) then
      B = A
      A = Aold
    endif

    enddo ! end loop over A and B

  endif ! ends if condition for whether observation exists

do i = 1, ndim
  spurious_all( i ) = stddev( B( i , : ) )
enddo
spurious_mean = mean( spurious_all(1:ndim) )
spurious = 1. / spurious_mean
av = sum( A , dim = 2 ) / real( nrens )
if (spurious_flag) then
do i = nsv+1,ndim
  do k = 1, nrens
    A( i , k ) = spurious * ( A( i , k ) - av(i) ) + av(i)
  enddo
enddo
endif

end subroutine enkf_calcs

!------------------------------------------

subroutine analysis( A , D , E , S , nrobs , x5 )
! Computes the analysed ensemble in the EnKF
! Written by G. Evensen ( Geir.Evensen@nersc.no )
! This routine uses subroutines from BLAS and EISPACK
! modified by M. Williams (addition of smoother calculations (X5))
! calls the additional multiplication routine multa

implicit none

! arguments
integer, intent(in) :: nrobs				! Number of observations
real, intent(inout) :: A(ndim, nrens)		! Ensemble matrix
real, intent(in)    :: D(nrobs,nrens)			! Matrix holding innovations
real, intent(in)    :: S (nrobs , nrens)		! Matrix holding HA�
real, intent(in)    :: E(nrobs , nrens)		! Matrix holding observation perturbations
real, intent(out)   :: X5(nrens , nrens)		! Matrix holding X5

! Local variables
real, allocatable, dimension (:)   :: sig, work
real, allocatable, dimension (:,:) :: X1 , X2 , U, X4 , Reps, Reps2
real ES (nrobs, nrens) , X3 (nrobs, nrens) , V(nrens, nrens)
real sigsum , sigsum1, Av(20)
integer nrsigma , i , j , nrmin , iblkmax, IRANK, lwork, ierr

!minimum of nrobs and nrens
nrmin=min(nrobs,nrens)  !why nrobs+1?

!compute HA'+E
ES=S+E

!compute SVD of HA'+E -> U and sig
allocate (U(nrobs,nrmin))
allocate (sig(nrmin))
lwork = 2 * max ( 3 * nrens + nrobs, 5 * nrens )
allocate (work(lwork))
sig=0.0

!call sgesvd ('S', 'N', nrobs, nrens, ES, nrobs, sig, U, nrobs, V, nrens, work, lwork, ierr) ! LAPACK routine

deallocate (work)
!
if (ierr.NE.0) then
	print*,'ierr from call sgesvd =', ierr
	stop
endif

!convert to eigenvalues
do i=1,nrmin
	sig(i)=sig(i)**2
enddo

!compute number of significant eigenvalues
sigsum=sum(sig(1:nrmin))
sigsum1=0.0
nrsigma=0
do i=1,nrmin
    if (sigsum1/sigsum<0.999)then
        nrsigma=nrsigma+1
        sigsum1=sigsum1+sig(i)
        sig(i)=1.0/sig(i)
    else
        sig(i:nrmin)=0.0
        exit
    endif
enddo

!compute X1
allocate(X1(nrmin,nrobs))
do j=1,nrobs
    do i=1,nrmin
        X1(i,j)=sig(i)*U(j,i)
    enddo
enddo
deallocate(sig)

!compute X2=X1*D
allocate(X2(nrmin,nrens))
!
!call sgemm('n' , 'n' , nrmin, nrens, nrobs, 1.0, X1, nrmin, D, nrobs, 0.0, X2, nrmin)	! BLAS routine

deallocate(X1)

!compute X3=U*X2
!call sgemm('n', 'n', nrobs, nrens, nrmin, 1.0, U, nrobs, X2, nrmin, 0.0, X3, nrobs)

deallocate(U)
deallocate(X2)

!compute final analysis

!compute X4=(HA')^T*X3
allocate(X4(nrens,nrens))
!call sgemm('t', 'n', nrens, nrens, nrobs, 1.0 , S, nrobs, X3, nrobs, 0.0, X4, nrens)

!compute X5=X4+I (stored in X4)
do i=1,nrens
	X4(i,i)=X4(i,i)+1.0
enddo
X5=X4

if(ndim>2000)then

!if(2*ndim*nrobs>nrens*(nrobs+ndim))then
!case with nrobs 'large'
!compute A=A*X5
	iblkmax=min(ndim,200)
	!call multa(A,X4,ndim,nrens,iblkmax)
	deallocate(X4)

else
!case with nrobs 'small'
!compute representers Reps=A'*S^T
	allocate (Reps(ndim,nrobs))
	!call sgemm('n', 't', ndim, nrobs, nrens, 1.0, A, ndim, S, nrobs, 0.0, Reps, ndim)
! Compute A = A + Reps * X3
	!call sgemm('n', 'n', ndim, nrens, nrobs, 1.0, Reps, ndim, X3, nrobs, 1.0, A, ndim)
	deallocate(Reps)

endif

end subroutine analysis

!---------------------------------------------------------------------------------------

subroutine varobs (maxobs,numobs) !determine number of observations this timestep

use obsdrivers, only : baseobs

implicit none

! arguments
integer, intent(in)   :: maxobs	! max number of observations
integer, intent (out) :: numobs	! actual number of observations

! local variables
integer :: j

numobs=0
do j=1,maxobs	!for each of the possible observations
	if (baseobs(j) .ge. -900.) then	!is there an observation?  very negative value means no
		numobs=numobs+1
	endif
enddo

end subroutine varobs

!------------------------------------------

subroutine masks( numobs , y , z , zbase , obloc )

! designed to construct the observation matrix z (aka H),
! according to the number of data available for each particular timestep.
! The full observation operator is defined by zbase
! M. Williams, 13 May 2003
! observation matrix referenced by z and H

use obsdrivers , only : baseobs

implicit none

! arguments
integer, intent(in)  :: numobs		! actual number of observations
real,dimension(numobs),intent(out)    :: y		! truncated observation file, missing data removed
real, intent(out)    :: Z(numobs,ndim)	! =H, the observation operator
real, intent(in)     :: ZBASE(ndim,maxobs)		!baseline observation operator, if maxobs data are available
integer, intent(out) :: obloc(numobs)	! location in the full obs matrix of the actual available obs

!local variables
integer :: j, i, x
real                           :: row(ndim)
logical,dimension(ndim,maxobs) :: maska , maskb
real,dimension(ndim,numobs) :: znew

maska = .false.  !initially fill the mask with FALSE values
obloc = 0

x = 1
znew = 0.		!initialise the final observation matrix
maskb = maska

DO j=1,maxobs	!for each of the possible observations
	IF(baseobs(j).ge.-900.)THEN	!is there an observation?  very negative value means no
		y(x)=baseobs(j)		!if so then load it into the y array
		obloc(x)=j		!record location of data point
		!first construct the observation matrix
		DO i=1,ndim		!define extraction row from zbase
			maskb(i,j)=.TRUE.	!fill the relevant row of the mask with TRUE values
		ENDDO
		row=PACK(zbase,maskb)	!extract the relevant row from the full obs. matrix
		maskb=maska		!reset the mask
		DO i=1,ndim		!define insertion row into znew
			maskb(i,x)=.TRUE.
		ENDDO
		znew=UNPACK(row,maskb,znew)	!insert the extracted row into the next empty row of the new observation matrix.
		maskb=maska		!reset the mask
		x=x+1
	ENDIF
ENDDO
z=transpose(znew)		!transpose the final matrix to give z, which is returned to the main programme

end subroutine masks

!------------------------------------------

subroutine smoother (X5,ct)		!see Evensen (2003) eqn 103

implicit none

! arguments
integer, intent(in) :: ct
real, intent(in)    :: X5(nrens, nrens)		! Matrix holding X5

!local variables
real    :: As(ndim, nrens)		! Smoothed Ensemble matrix
integer :: j,start,lag

lag=30
start=max(1,ct-lag)

do j=start,ct
  !read(unit=16,rec=j)As					! read smoothed ensemble
  As=matmul(As,X5)					!eqn 103
  !write(unit=16,rec=j)As					! write re-smoothed ensemble
enddo

end subroutine smoother

!------------------------------------------

subroutine smoothout()  ! smoother output

  implicit none

  !local variables
  real       :: As(ndim, nrens)		! Smoothed Ensemble matrix
  real       :: Aplus(ndim+1, nrens)	! Ensemble matrix with extra space for NEE calculations (deleted Rtot, changed 2 to 1.)
  real       :: Ast(nrens,ndim+1)		! Transpose of Aplus, for stats
  real       :: stat(2,ndim+2) , statn(2,1)			! Matrix holding stats on analysis
  real, save :: NEEs(nrens)
  integer    :: i , k , j , NRMISS

  NEEs = 0.

  do j = 1 , totaltime
    read(16,rec=j)As					! read smoothed ensemble
    !determine modelled NEE and Rtot for each member of the ensemble
    Aplus=0.
    do i = 1 , ndim
      do k = 1 , nrens
        Aplus( i , k ) = As( i , k )
      enddo
    enddo

    do k = 1 , nrens
!		Aplus(ndim+1,k)=As(1,k)+As(11,k)+As(12,k)	!Rtot
      Aplus( ndim + 1 , k ) = Aplus( ndim + 1 , k ) - As( 1 , k ) - As( 2 , k )	!NEE
    enddo
    !determine accumulated NEE

    do k = 1 , nrens
      NEEs( k ) = NEEs( k ) + Aplus ( ndim + 1 , k )	!NEE
    enddo

    Ast = transpose( Aplus )

    ! do stats on Ast
    do i = 1 , ndim + 1
      stat( 1 , i ) = SUM( AsT( : , i ) ) / nrens !  stat(1,i) ! compute mean for each column of AT
      stat( 2 , i ) = SQRT( SUM( ( AsT( : , i ) - stat( 1 , i ) ) ** 2 ) / nrens ) ! compute standard deviation for each column of AT
    enddo

    !write(17,'(I6,",",50(F9.3,","))')j,(stat(1,k),k=1,ndim+2)
    !write(20,'(I6,",",50(F9.3,","))')j,(stat(2,k),k=1,ndim+2)	!st.devs.
    !write(48,'(I6,",",50(F9.3,","))')j,(stat(1,k),k=1,ndim+2),(stat(2,k),k=1,ndim+2)
  enddo

  statn( 1 , 1 ) = SUM( NEEs ) / nrens !  stat(1,i) ! compute mean for NEEs
  statn( 2 , 1 ) = SQRT( SUM( ( NEEs - statn( 1 , 1 ) ) ** 2 ) / nrens ) ! compute standard deviation for NEEs

  !write(46,'(I6,",",20(F9.3,","))')j,statn(1,1),statn(2,1)

end subroutine smoothout

!------------------------------------------

subroutine covariancemat(A,covmat,cormat)	!calculate covariance matrix (see Evensen 2004, p.540)

use clim

implicit none

! arguments
real, intent(in)  :: A(ndim, nrens)	! Ensemble matrix
real, intent(out) :: covmat(ndim,ndim), cormat(ndim,ndim)	!covariance and correlation matrices

! local variables
real    :: pertmat(ndim,nrens)	!perturbation matrix
real    :: oneN(nrens,nrens),abar(ndim,nrens)
integer :: i,j

oneN=1./nrens
Abar=matmul(a,oneN)
pertmat=A-Abar
covmat=matmul(pertmat,transpose(pertmat))/(nrens-1.)
Do j=1,ndim
	Do i=1,ndim
		cormat(i,j)=covmat(i,j)/(covmat(i,i)**0.5*covmat(j,j)**0.5)
	enddo
enddo

DO i=1,ndim
    !write(164,'(A10,",",50(F12.5,","))')TRIM(SV_names(i)),(covmat(i,j),j=1,ndim)
    !write(165,'(A10,",",50(F12.5,","))')TRIM(SV_names(i)),(cormat(i,j),j=1,ndim)
ENDDO
parcor(1)  = cormat(nsv+1, nsv+2) ! iota and gplant
parcor(2)  = cormat(nsv+1, nsv+3) ! iota and rootresist
parcor(3)  = cormat(nsv+1, nsv+4) ! iota and minLWP
parcor(4)  = cormat(nsv+1, nsv+5) ! iota and root_leaf_ratio
parcor(5)  = cormat(nsv+1, nsv+6) ! iota and rbelow
parcor(6)  = cormat(nsv+2, nsv+3) ! gplant and rootresist
parcor(7)  = cormat(nsv+2, nsv+4) ! gplant and minLWP
parcor(8)  = cormat(nsv+2, nsv+5) ! gplant and root_leaf_ratio
parcor(9)  = cormat(nsv+2, nsv+6) ! gplant and rbelow
parcor(10) = cormat(nsv+3, nsv+4) ! rootresist and minLWP
parcor(11) = cormat(nsv+3, nsv+5) ! rootresist and root_leaf_ratio
parcor(12) = cormat(nsv+3, nsv+6) ! rootresist and rbelow
parcor(13) = cormat(nsv+4, nsv+5) ! minLWP and root_leaf_ratio
parcor(14) = cormat(nsv+4, nsv+6) ! minLWP and rbelow
parcor(15) = cormat(nsv+5, nsv+6) ! root_leaf_ratio and rbelow

transparcor(1) = cormat(nsv-1, nsv+1) ! trans and iota
transparcor(2) = cormat(nsv-1, nsv+2) ! trans and gplant
transparcor(3) = cormat(nsv-1, nsv+3) ! trans and rootresist
transparcor(4) = cormat(nsv-1, nsv+4) ! trans and minLWP
transparcor(5) = cormat(nsv-1, nsv+5) ! trans and root_leaf_ratio
transparcor(6) = cormat(nsv-1, nsv+6) ! trans and rbelow

capacitanceparcor(1) = cormat(nsv, nsv+1) ! capacitance and iota
capacitanceparcor(2) = cormat(nsv, nsv+2) ! capacitance and gplant
capacitanceparcor(3) = cormat(nsv, nsv+3) ! capacitance and rootresist
capacitanceparcor(4) = cormat(nsv, nsv+4) ! capacitance and minLWP
capacitanceparcor(5) = cormat(nsv, nsv+5) ! capacitance and root_leaf_ratio
capacitanceparcor(6) = cormat(nsv, nsv+6) ! capacitance and rbelow

!write(164,'(A50,I3)')"-----------------------------------------------",yearday
!write(165,'(A50,I3)')"-----------------------------------------------",yearday

end subroutine covariancemat

!------------------------------------------

subroutine obsop( Hbase )	!observation operator

implicit none

! arguments
real, intent(out) :: Hbase(ndim,maxobs)	! baseline observation operator

Hbase = RESHAPE ((/&
0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,0.,0.,0.,0.,0.,0.,0. &
/), (/ndim, maxobs/))

end subroutine obsop

!------------------------------------------

subroutine errors( sigma2 , measvar , av )	!errors on model and observations

use obsdrivers, only : baseobs

implicit none

real, intent(out) :: sigma2(ndim)			!model error variance
real, intent(out) :: measvar(maxobs)		!measurement error variance
real, intent(in)  :: av(ndim)
!real, intent(in) :: baseobs(maxobs)
integer i,j

! error variance in the model
DO i=1,ndim
	IF(abs_errors(i).eq.0)THEN
		sigma2(i)=(modvar(i)*abs(av(i)))**2		!square the RHS
	ELSE
		sigma2(i)=modvar(i)**2
	ENDIF
ENDDO

! error in the observation equation.
DO i=1,maxobs
    IF (obserr_abs) THEN
      measvar( i ) = obsvar( i ) ** 2
    ELSE
      measvar( i ) = ( obsvar( i ) * abs( baseobs( i ) ) ) **2 !stand. error
    ENDIF
ENDDO

!write(57,'(20(F9.3,","))')(measvar(i)**.5,i=1,2)

end subroutine errors

!------------------------------------------

  subroutine biastest(mtest)

    use biasmod

    implicit none

    ! arguments
    real, intent(out) :: mtest
    ! local variables
    integer           :: i
    real              :: ro, totaln , meaninn
    real              :: inn( totaltime )

    ro = 0. ; inn = 0. ; totaln = 0.

    !calculate NEE innovation
    do i = 1 , totaltime
      if ( neeobs( i ) .gt. -900 ) then
        inn( i ) = forecastmean( i , ndim + 1 ) - neeobs( i )
        totaln = totaln + 1.
      endif
    enddo
    meaninn = sum( inn ) / totaln	!mean innovation

    do i = 1 , totaltime
      if ( neeobs( i ) .gt. -900 ) then
        ro = ro + ( inn( i ) - meaninn ) ** 2
      endif
    enddo
    ro = ro / totaln
    Mtest = meaninn * totaln ** 0.5 / ( ro ** 0.5 )

  end subroutine biastest

!------------------------------------------

  subroutine update_ensemble_variables ( i , load , frstrn , SWC , LWPout , &
    LWPout_sd , update , runoff_out , BGtot , BG1 , BG2 , CSR , rsoilout , rplantout)

    use veg
    use soil_structure
    use snow_info
    use clim
    use hourscale
    use irradiance_sunshade
    use metab
    use meteo
    use scale_declarations !, only : nos_canopy_layers , steps , user_opts , nos_soil_layers , core
    use math_tools

    implicit none

    ! arguments
    integer,intent(in) :: i ! ensemble member subscript
    logical,intent(in) :: load ! switch determining whether ensemble value is loaded into variable
    logical,intent(in) :: frstrn ! switch determining whether model has looped over first timestep or not
    logical,intent(in) :: update ! switch determining whether update is performed
    real,dimension(core),intent(out)   :: SWC,BGtot,BG1,BG2
    real,dimension(nos_canopy_layers),intent(out) :: LWPout,LWPout_sd,CSR,rsoilout,rplantout
    real,intent(out) :: runoff_out

    ! local variables
    integer :: j

    ! module math_tools variables
    real,dimension(nmax,kmaxx,nrens),save :: yp_ens
    real,dimension(kmaxx,nrens),save      :: xp_ens
    real,dimension(nrens),save            :: dxsav_ens
    ! module veg variables
    real,dimension(nrens),save ::  co2amb_ens , &
                                    prevC_ens , &
                                  totevap_ens , &
                                     totn_ens
    real,dimension(nos_canopy_layers , nrens),save :: &
                            LWPstore_ens = 0 , &
              canopy_soil_resistance_ens , &
                                 lai_ens , &
                                 Nla_ens
    real,dimension(steps , nrens),save :: ess_ens , &
                                gppt_ens , &
                               respt_ens , &
                            soiletmm_ens , &
                              transt_ens
    real,dimension(steps , nos_canopy_layers , nrens ),save :: flux_ens

    ! module soil_structure variables
    real,dimension(nrens),save :: drythick_ens , &
                               rootbiomass_ens , &
                                 rootreach_ens , &
                               surfbiomass_ens , &
                              weighted_SWP_ens , &
                                   thermal_ens
    real,dimension(nos_soil_layers , nrens),save :: watericemm_ens
    real,dimension(core , nrens),save :: conduc_ens , &
                                          cond1_ens , &
                                          cond2_ens , &
                                          cond3_ens , &
                                        iceprop_ens , &
                                           potA_ens , &
                                           potB_ens , &
                                        pptgain_ens , &
                                     rootlength_ens , &
                                       rootmass_ens , &
                                       soiltemp_ens , &
                                soiltemp_nplus1_ens , &
                                            SWP_ens , &
                                          soilR_ens , &
                                         soilR1_ens , &
                                         soilR2_ens , &
                                      waterfrac_ens , &
                                      waterloss_ens , &
                                      watergain_ens
    real,dimension(10 , nrens),save :: wettingtop_ens = 0. , wettingbot_ens = 0.

    ! module snow_info variables
    real,dimension(nrens),save :: snow_watermm_ens

    ! module clim variables
    real,dimension(nrens),save :: daytempsum_ens , &
                                 daypar_ens , &
                                 dayppt_ens , &
                                    gdd_ens , &
                                max_fol_ens , &
                                 multtf_ens , &
                                 multtl_ens
    real,dimension(steps , nrens),save :: wetev_ens

    ! module hourscale variables
    real,dimension(nrens),save :: canopy_store_ens , &
                                       dayevap_ens , &
                                     discharge_ens , &
                                    evap_store_ens , &
                                      overflow_ens , &
                                        runoff_ens , &
                               surface_watermm_ens , &
                                         totet_ens , &
                                     underflow_ens , &
                                 unintercepted_ens

    ! module irradiance_sunshade variables
    real,dimension(nrens),save :: soilnet_ens

    ! module metab variables
    real,dimension(nrens),save :: resp_ens , &
                                 rsoil_ens
    real,dimension(nos_canopy_layers,nrens) :: rsoil_canopy_ens , &
                                               rplant_canopy_ens

    ! module meteo variables
    real,dimension(nrens),save :: la_ens , &
                                 nit_ens , &
                                 par_ens , &
                                psil_ens , &
                                psis_ens , &
                                 rad_ens

    if ( load .eqv. .true.) then !load the variable states from ensemble into current variable

    ! On first hour of day, set counters to zero..
    if (( time%step .eq. 1 ).and.( i .eq. 1 )) then
      daypar_ens  = 0. ; dayppt_ens  = 0. ; discharge_ens = 0. ; flux_ens   = 0.
      gppt_ens    = 0. ; respt_ens   = 0. ; runoff_ens    = 0.
      totevap_ens = 0. ; transt_ens    = 0. ; wetev_ens  = 0.
    endif

    if (update) then
        !======================================= math tools variables
        yp = yp_ens( : , : , i ) ; xp = xp_ens( : , i ) ; dxsav = dxsav_ens( i )
        !======================================= veg variables
        co2amb = co2amb_ens( i ) ; prevC = prevC_ens( i ) ; totevap = totevap_ens( i )
        totn   = totn_ens  ( i )
        LWPstore = LWPstore_ens( : , i ) ; canopy_soil_resistance = canopy_soil_resistance_ens( : , i )
        lai      = lai_ens     ( : , i ) ; Nla                    = Nla_ens                   ( : , i )
        ess      = ess_ens     ( : , i ) ; gppt                   = gppt_ens                  ( : , i )
        respt    = respt_ens   ( : , i ) ; soiletmm               = soiletmm_ens              ( : , i )
        transt   = transt_ens  ( : , i )
        flux     = flux_ens ( : , : , i )

        !======================================= soil_structure variables
        drythick    = drythick_ens   ( i ) ; rootbiomass  = rootbiomass_ens ( i )
        surfbiomass = surfbiomass_ens( i ) ; weighted_SWP = weighted_SWP_ens( i )
        rootreach   = rootreach_ens  ( i ) ; thermal      = thermal_ens     ( i )
        watericemm  = watericemm_ens ( : , i ) ; conduc     = conduc_ens     ( : , i )
        cond1       = cond1_ens      ( : , i ) ; cond2      = cond2_ens      ( : , i )
        cond3       = cond3_ens      ( : , i ) ; iceprop    = iceprop_ens    ( : , i )
        potA        = potA_ens       ( : , i ) ; potB       = potB_ens       ( : , i )
        pptgain     = pptgain_ens    ( : , i ) ; rootlength = rootlength_ens ( : , i )
        rootmass    = rootmass_ens   ( : , i ) ; soiltemp   = soiltemp_ens   ( : , i )
        SWP         = SWP_ens        ( : , i ) ; waterfrac  = waterfrac_ens  ( : , i )
        waterloss   = waterloss_ens  ( : , i ) ; watergain  = watergain_ens  ( : , i )
        wettingtop  = wettingtop_ens ( : , i ) ; wettingbot = wettingbot_ens ( : , i )
        soiltemp_nplus1 = soiltemp_nplus1_ens ( : , i )

        !======================================= snow_info variables
        snow_watermm = snow_watermm_ens( i )

        !======================================= clim variables
        daytempsum = daytempsum_ens( i ) ; daypar  = daypar_ens ( i ) ; dayppt = dayppt_ens( i )
        gdd        = gdd_ens       ( i ) ; max_fol = max_fol_ens( i ) ; multtf = multtf_ens( i )
        multtl     = multtl_ens    ( i )
        wetev      = wetev_ens     ( : , i )

        !======================================= hourscale variables
        canopy_store    = canopy_store_ens   ( i ) ; dayevap       = dayevap_ens      ( i )
        discharge       = discharge_ens      ( i ) ; evap_store    = evap_store_ens   ( i )
        overflow        = overflow_ens       ( i ) ; runoff        = runoff_ens       ( i )
        surface_watermm = surface_watermm_ens( i ) ; totet         = totet_ens        ( i )
        underflow       = underflow_ens      ( i ) ; unintercepted = unintercepted_ens( i )

        !======================================= irradiance_sunshade variables
        soilnet = soilnet_ens( i )

        !======================================= metab variables
        resp = resp_ens( i ) ; rsoil = rsoil_ens( i )
        rsoil_canopy = rsoil_canopy_ens( : , i )
        rplant_canopy = rplant_canopy_ens( : , i )

        !======================================= meteo variables
        la   = la_ens  ( i ) ; nit = nit_ens( i ) ; par = par_ens( i ) ; psil = psil_ens( i )
        psis = psis_ens( i ) ; rad = rad_ens( i )

      endif ! ends if condition for whether update should be performed

    else !save the variable states into corresponding ensemble member

        !======================================= math tools variables
        yp_ens( : , : , i ) = yp ; xp_ens( : , i ) = xp ; dxsav_ens( i ) = dxsav
        !======================================= veg variables
        co2amb_ens  ( i )     = co2amb ; prevC_ens  ( i ) = prevC
        totn_ens    ( i )     = totn   ; totevap_ens( i ) = totevap
        LWPstore_ens( : , i ) = LWPstore ; canopy_soil_resistance_ens( : , i ) = canopy_soil_resistance
        lai_ens     ( : , i ) = lai      ; Nla_ens                   ( : , i ) = Nla
        ess_ens     ( : , i ) = ess      ; gppt_ens                  ( : , i ) = gppt
        respt_ens   ( : , i ) = respt    ; soiletmm_ens              ( : , i ) = soiletmm
        transt_ens  ( : , i ) = transt
        flux_ens    ( : , : , i ) = flux

        !======================================= soil_structure variables
        drythick_ens   ( i ) = drythick    ; rootbiomass_ens ( i ) = rootbiomass
        surfbiomass_ens( i ) = surfbiomass ; weighted_SWP_ens( i ) = weighted_SWP
        rootreach_ens  ( i ) = rootreach   ; thermal_ens     ( i ) = thermal
        watericemm_ens     ( : , i ) = watericemm ; conduc_ens    ( : , i ) = conduc
        cond1_ens          ( : , i ) = cond1      ; cond2_ens     ( : , i ) = cond2
        cond3_ens          ( : , i ) = cond3      ; iceprop_ens   ( : , i ) = iceprop
        potA_ens           ( : , i ) = potA       ; potB_ens      ( : , i ) = potB
        pptgain_ens        ( : , i ) = pptgain    ; rootlength_ens( : , i ) = rootlength
        rootmass_ens       ( : , i ) = rootmass   ; soiltemp_ens  ( : , i ) = soiltemp
        SWP_ens            ( : , i ) = SWP        ; waterfrac_ens ( : , i ) = waterfrac
        waterloss_ens      ( : , i ) = waterloss  ; watergain_ens ( : , i ) = watergain
        wettingtop_ens     ( : , i ) = wettingtop ; wettingbot_ens( : , i ) = wettingbot
        soiltemp_nplus1_ens( : , i ) = soiltemp_nplus1
        soilR1_ens         ( : , i ) = soilR1     ; soilR2_ens    ( : , i ) = soilR2
        soilR_ens          ( : , i ) = soilR

        !======================================= snow_info variables
        snow_watermm_ens( i ) = snow_watermm

        !======================================= clim variables
        daytempsum_ens( i ) = daytempsum ; daypar_ens ( i ) = daypar  ; dayppt_ens( i ) = dayppt
        gdd_ens       ( i ) = gdd        ; max_fol_ens( i ) = max_fol ; multtf_ens( i ) = multtf
        multtl_ens    ( i ) = multtl
        wetev_ens ( : , i ) = wetev

        !======================================= hourscale variables
        canopy_store_ens   ( i ) = canopy_store    ; dayevap_ens      ( i ) = dayevap
        discharge_ens      ( i ) = discharge       ; evap_store_ens   ( i ) = evap_store
        overflow_ens       ( i ) = overflow        ; runoff_ens       ( i ) = runoff
        surface_watermm_ens( i ) = surface_watermm ; totet_ens        ( i ) = totet
        underflow_ens      ( i ) = underflow       ; unintercepted_ens( i ) = unintercepted

        !======================================= irradiance_sunshade variables
        soilnet_ens( i ) = soilnet

        !======================================= metab variables
        resp_ens( i ) = resp ; rsoil_ens( i ) = rsoil
        rsoil_canopy_ens( : , i ) = rsoil_canopy
        rplant_canopy_ens( : , i ) = rplant_canopy

        !======================================= meteo variables
        la_ens  ( i ) = la   ; nit_ens( i ) = nit ; par_ens( i ) = par ; psil = psil_ens( i )
        psis_ens( i ) = psis ; rad_ens( i ) = rad

        if ( i .EQ. nrens ) then
          do j = 1 , core
            SWC( j )   = mean( waterfrac_ens( j , : ) )
            BGtot( j ) = mean( soilR_ens( j , : ) )
            BG1( j )   = mean( soilR1_ens( j , : ) )
            BG2( j )   = mean( soilR2_ens( j , : ) )
          enddo
          ! write LWP for each layer to output files, averaged over ensemble
          do j = 1 , nos_canopy_layers
            LWPout( j ) = mean(LWPstore_ens( j , : ))
            LWPout_sd( j ) = stddev(LWPstore_ens( j , : ))
            CSR( j ) = mean( canopy_soil_resistance_ens( j , : ))
            rsoilout( j ) = mean(rsoil_canopy_ens( j , : ))
            rplantout( j ) = mean(rplant_canopy_ens( j , : ))
          enddo
          runoff_out = mean( runoff_ens( : ) )
        endif

    endif ! ends if condition for loading to or from ensemble variables

  end subroutine update_ensemble_variables

!------------------------------------------

  subroutine transform_parameters( A , trans , time , memb )

    use math_tools
    use scale_declarations,only: time_holder

    implicit none

    ! arguments
    real, intent(inout) :: A(ndim, nrens)		! Ensemble matrix
    logical             :: trans
    type(time_holder),intent(inout) :: time
    integer,intent(in) :: memb

    ! local variables
    integer                     :: j , k , i
    real                        :: reduce = 0.47  ! ensemble spread reduction factor, see value in Chen08, p 319
    real                        :: av , astd , minA , maxA
    real,dimension(nrens)       :: R
    integer,dimension(4),save   :: iseed
    real,dimension(ndim-nsv)    :: minerrr
    logical :: delete

    if ( time%steps_count .eq. 1) iseed = (/10,300,500,701/)

    i = 0

    do j = nsv + 1 , ndim

      i = i + 1

      !call slarnv( 3, iseed, nrens, R ) !random numbers from a normal distribution
      R = R - mean( R )
      R = R / stddev( R )

      !do k = 1 , nrens
        ! transform FROM original parameter values (i.e. after roots and timestep_calcs, trans is .true.)
        if ( trans ) then
          A( j , memb ) = asin( sqrt( ( A( j , memb ) - lo( j ) ) / ( hi( j ) - lo( j ) ) ) )
        ! transform TO original parameter values (i.e. before roots and timestep_calcs, trans is .false.)
        else
          A( j , memb ) = lo( j ) + ( hi( j ) - lo( j ) ) * ( sin( A( j , memb ) ) ) ** 2
        endif
      !enddo

      if (trans .eqv. .false.) then
        !inflate_exe( i ) = 0.
        minerrr = minerr
        av = mean( A( j , : ) ) ; astd = stddev( A( j , : ) )
        if (inflate_flag) then
          if ( astd .lt. minerr( i ) ) then
            !inflate_exe( i ) = 1.
            !inflate = 1.01 ! when parameter stddev lower than minerr, inflation factor is greater than 1 (1.01 as in EvensenBook, p 239)
            inflate = minerr( i ) / astd ! inflation factor set to inflate standard error back to minimum value, see Aksoy06
            do k = 1 , nrens
              A( j , k ) = inflate * ( A( j , k ) - av ) + av
            enddo
          else
            !inflate = sqrt( 1. - ( reduce ** 2. ) )
          endif
          !inflate = 0.
        endif ! closes if condition on inflate flag
        av = mean( A( j , : ) ) ; astd = stddev( A( j , : ) )
        !do k = 1 , nrens
          !A( j , k ) =  ( reduce * A( j , k ) ) + ( ( 1. - reduce ) * av ) + ( A( j , k ) * R( k ) * inflate * astd )
          if ( A( j , memb ) .lt. lo( j ) ) A( j , memb ) = lo( j )  ! parameters to be within preset limits (lo/hi)
          if ( A( j , memb ) .gt. hi( j ) ) A( j , memb ) = hi( j )
        !enddo
        minA = minval( A( j , : ) )
        !if ( minA .lt. lo( j ) ) A( j , : ) = A( j , : ) + lo( j ) - minA
        maxA = maxval( A( j , : ) )
        !if ( maxA .gt. hi( j ) ) A( j , : ) = A( j , : ) + hi( j ) - maxA
      endif

    enddo

  end subroutine transform_parameters

end module enkf
