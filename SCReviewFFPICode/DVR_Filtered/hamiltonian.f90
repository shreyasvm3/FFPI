! This program diagonalizes the Hamiltonian
! It outputs the eigenvectors in H.dat
!   and the eigenvalues in Energy.dat
!

program energymatrix
use parameters

integer                 :: i, j, k
real*8                  :: VR, dV

! Hamiltonian matrix and eigenvalues
real*8, allocatable     :: H(:,:), Energy(:)

! Diagonalization variables
integer                 :: lwork, liwork, info
real*8, allocatable     :: work(:)
integer, allocatable    :: iwork(:)

call cpu_time(tstart)

call input

liwork = 5*NR + 3
lwork  = 2*NR**2 + 6*NR + 1

allocate(H(NR,NR),Energy(NR),work(lwork),iwork(liwork))

H = 0.d0

! Diagonal elements of Hamiltonian
do i = 1, NR
  call VofR(gridR(i),VR,dV)
  H(i,i) = (pi)**2/(6.d0*Mass*dR**2) + VR
enddo
! Upper off-diagonal elements of Hamiltonian
do i = 1, NR-1
  do j = i+1, NR
    H(i,j) = (-1)**(i-j)/(Mass*(dR*(i-j))**2)
  enddo
enddo
! Lower off-diagonal elements of Hamiltonian
do i = 1, NR-1
   do j = i+1, NR
      H(j,i)=H(i,j)
   enddo
enddo

! Diagonalization     
call dsyevd('V','U',NR,H,NR,Energy,work,lwork,iwork,liwork,info)

if (info.ne.0) write(*,*) 'Problems with diagonalization'

deallocate(work, iwork)

call cpu_time(tfinish)

write(*,*) 'Diagonalization time (mins): ', (tfinish-tstart)/60.d0

open(800,file='Energy.dat',status='unknown')
open(900,file='H.dat',status='unknown')

do i = 1, NR
   write(800,*) Energy(i)
   do j = 1, NR
      write(900,*) H(i,j)
   enddo
enddo

close(800)
close(900)

end program
