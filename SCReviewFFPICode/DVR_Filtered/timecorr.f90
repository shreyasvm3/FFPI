! This program uses the eigensystem generated from 
! hamiltonian.f90 and computes the time-dependent
! position correlation function
!

program dynamics
use parameters
implicit none

include 'mpif.h'

real*8, allocatable     :: H(:,:), Energy(:)
real*8, allocatable     :: Amat(:), Bmat(:)
complex*16, allocatable :: FpropXdt(:,:), BpropXdt(:,:), PropPMdt(:,:)
complex*16, allocatable :: Rt(:), AverR(:), Prop(:,:), Filpropdt(:,:)
real*8                  :: cf0
integer                 :: i, j
integer                 :: ierr, myrank, nprocs, istart, iend

call input
! Setting up the parallelization
call mpi_init(ierr) !<--- Create child processes
call mpi_comm_size(mpi_comm_world,nprocs,ierr) !<--- Find the total number of processes
call mpi_comm_rank(mpi_comm_world,myrank,ierr) !<--- Find the rank of each process
call para_range(1,NTimePoints,nprocs,myrank,istart,iend) !<--- Divide configurations between processes
call mpi_barrier(mpi_comm_world,ierr)

allocate(Amat(NXpm),Bmat(NXpm))
allocate(FpropXdt(NR,NR),BpropXdt(NR,NR))
allocate(Prop(NXpm,NXpm))
allocate(PropPMdt(NXpm,NXpm),FilPropdt(NXpm,NXpm))
allocate(H(NR,NR),Energy(NR))

! Read Hamiltonian matrix and energies
open (800,file='Energy.dat',status='old')
open (600,file='H.dat',status='old')
   do i = 1, NR
      read(800,*) Energy(i)
      do j = 1, NR
         read(600,*) H(i,j)
     enddo
   enddo
close(800)
close(600)
! The vector A in the PM basis
call Amatrix(Amat)
! The matrix B in the PM basis
call Bmatrix(Bmat)
allocate(Rt(0:NTimePoints),AverR(0:NTimepoints))

cf0       = 0.d0
Rt        = 0.d0
AverR     = 0.d0
FpropXdt  = 0.d0
BpropXdt  = 0.d0
PropPMdt  = 0.d0
FilPropdt = 0.d0
Prop      = 0.d0


! Calculate the zero-time cf
if (myrank.eq.0) then
   call cf0time(Amat,Bmat,Rt(0))
end if
! Calculate th F&B props in X basis (for eps/2 propagation)
call STPropXbasis(H,Energy,FpropXdt,BpropXdt)
! Calculate the short time (eps/2) Prop in PM basis
call STPropPMbasis(FpropXdt,BpropXdt,PropPMdt) 
! Calculate the filtered propagator
call FilteredProp(PropPMdt,FilPropdt)
! Calculate the Fil prop for t = istart
call Propst(istart,FilPropdt,Prop)
!Get the c.f. for t= istart
call corrf(Amat,Bmat,Prop,Rt(istart))
do j = istart+1, iend
   Nts = j
   ! Calculate the full Prop by adding eps prop
   call Fullprop(FilPropdt,Prop)
   ! Calculate the corr. function 
   call corrf(Amat,Bmat,Prop,Rt(j))
enddo
call mpi_barrier(mpi_comm_world,ierr)
! Collecting data from different processes
call mpi_reduce(Rt,AverR,NTimePoints+1,mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)


if (myrank.eq.0) then
   open(555,file='R.out',status='unknown')
      cf0 = real(AverR(0))
      do j = 0, NTimePoints
         write(555,*) j*Eps, real(AverR(j))/cf0, aimag(AverR(j))/cf0
      enddo
   close(555)
endif

call mpi_finalize(ierr)

end program
