module parameters
implicit none

! Pi and imaginary unit
real*8, parameter     :: pi=dacos(-1.d0)
complex*16, parameter :: Iu=(0.d0,1.d0)

! Time grid
integer :: Nts, NTimePoints
real*8  :: Time, Eps

! System parameters
real*8  :: Mass, Width, qIn, pIn

! DVR grid
integer :: NR, NXpm
real*8  :: Rmin, Rmax
real*8  :: dR
real*8, allocatable     :: gridR(:)
real*8, allocatable     :: gridxP(:), gridxM(:)
real*8, allocatable     :: gridY(:), gridZ(:)
integer, allocatable    :: Pidx(:), Midx(:)

! Filinov parameters
real*8  :: cFil, dFil
end module

subroutine input
use parameters
implicit none

integer         :: i, j, k
character*75    :: infostr

open(555,file='input',status='old')

read(555,'(a75)') infostr
read(555,*)       Mass

read(555,'(a75)') infostr
read(555,*)       Rmin, Rmax, NR

read(555,'(a75)') infostr
read(555,*)       Time, Eps

NTimePoints = int(Time/Eps)

read(555,'(a75)') infostr
read(555,*)       qIn, pIn, Width

read(555,'(a75)') infostr
read(555,*)       cFil, dFil

close(555)

! Dimension of xPM grid
NXpm =  NR*NR

! Spacing in x grid
dR = (Rmax-Rmin)/dble(NR-1)
allocate(gridR(NR))

allocate(gridxP(NXpm),gridxM(NXpm))
allocate(Pidx(NXpm),Midx(NXpm))
allocate(gridY(NXpm),gridZ(NXpm))

!Build x grid
do i = 0, NR-1
  gridR(i+1) = Rmin + dble(i)*dR
enddo

!Build x+ and x- grid and indices
do i = 1, NR
  do j = 1, NR
     k = (i-1)*NR+j
     GridxP(k) = GridR(i)
     GridxM(k) = GridR(j)
     Pidx(k)   = i
     Midx(k)   = j
     GridY(k)  = (GridxP(k)+GridxM(k))/2.d0
     GridZ(k)  = GridxP(k)-GridxM(k)
   end do
end do
 

end subroutine
