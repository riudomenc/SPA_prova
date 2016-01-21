module m_multa

  contains

  subroutine multa(A,X,ndim,nrens,iblkmax)

    implicit none

    integer, intent (in) :: ndim
    integer, intent (in) :: nrens
    integer, intent (in) :: iblkmax
    real, intent (in) :: X(nrens,nrens)
    real, intent (inout)::A(ndim,nrens)
    real v (iblkmax,nrens)	!automatic work array
    integer ia, ib

    do ia = 1,ndim,iblkmax
      ib=min(ia+iblkmax-1,ndim)
      v(1:ib-ia+1,1:nrens)=A(ia:ib,1:nrens)
      !call sgemm( 'N' , 'N' , ib-ia+1 , nrens , nrens , 1.0 , v(1,1) , iblkmax , X(1,1) , nrens , 0.0 , A(ia,1) , ndim )
    enddo

  end subroutine multa

end module m_multa

!-----------------------------------------------------------------------------------------------------------------------------------

module obsdrivers

  implicit none

  real                            :: dectime,obssd(2)
  real,allocatable,dimension(:,:) :: all_baseobs
  real,allocatable,dimension(:)   :: baseobs
  real,dimension(3)               :: sapobs
  integer                         :: timeid,newday,hr

end module obsdrivers

!-----------------------------------------------------------------------------------------------------------------------------------

module biasmod
IMPLICIT NONE
real :: neeobs(1096),neetot,neeerr
end module biasmod

!-----------------------------------------------------------------------------------------------------------------------------------
