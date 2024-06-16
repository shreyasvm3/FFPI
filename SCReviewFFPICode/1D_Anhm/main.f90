!###############################################################
!       This file contains the main program that calculates the
!       correlation function using all the other modules.
!###############################################################

program YZTraj_TCF
use parameters
use MD
use Estimator
use stat
implicit none

include 'mpif.h'

real*8,allocatable      :: pY(:), yTraj(:,:), zTraj(:,:), pZ(:)
integer                 :: MDBreakC, MDTrajC, MDBreakT, MDTrajT
integer                 :: ierr, myrank, nprocs, MDFlag
integer                 :: istart, iend, i, j, k
complex*16,allocatable  :: AvgCorrf(:,:), NormInt(:,:)
complex*16,allocatable  :: FAvgCorrf(:,:), FNormInt(:,:)


call input()                                                    ! Read input file on root                      
allocate(yTraj(Ndof,0:Nts),zTraj(Ndof,0:Nts))                   ! Allocate arrays                
allocate(pY(Ndof),pZ(Ndof))
allocate(AvgCorrf(3,0:Nts),NormInt(3,0:Nts))
allocate(FAvgCorrf(3,0:Nts),FNormInt(3,0:Nts))


call mpi_init(ierr)                                             ! Create child processes
call mpi_comm_size(mpi_comm_world,nprocs,ierr)                  ! Find the total number of processes
call mpi_comm_rank(mpi_comm_world,myrank,ierr)                  ! Find the rank of each process
call init_random_seed
call para_range(1,NTraj,nprocs,myrank,istart,iend) 

AvgCorrf        = 0.d0
FAvgCorrf       = 0.d0
NormInt         = 0.d0
FNormInt        = 0.d0
yTraj           = 0.d0
zTraj           = 0.d0
pY              = 0.d0
pZ              = 0.d0

!{Loop over 3 separate calculations to get error bars
do i = 1, 3
   MDFlag          = 0
   MDBreakC        = 0
   MDTrajC         = 0
   MDTrajT         = 0
   MDBreakT        = 0
   call mpi_barrier(mpi_comm_world,ierr)
   do j = istart, iend
   101   continue
      call sampleIC(yTraj(:,0),pY,zTraj(:,0),pZ)
      call PropYZ(yTraj(:,0),pY,zTraj(:,0),pZ,yTraj,zTraj,MDBreakC,MDTrajC,MDFlag)
      if (MDflag.eq.1) goto 101
      call Est(yTraj,zTraj,AvgCorrf(i,:),NormInt(i,:))
   end do
   call mpi_barrier(mpi_comm_world,ierr)
end do
call mpi_reduce(AvgCorrf,FAvgCorrf,3*(Nts+1),mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(NormInt,FNormInt,3*(Nts+1),mpi_double_complex,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(MDBreakC,MDBreakT,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)
call mpi_reduce(MDTrajC,MDTrajT,1,mpi_integer,mpi_sum,0,mpi_comm_world,ierr)

!Final calculations only on the root.
if (myrank.eq.0) then
  open(1,file='TCF.out',status='replace')
  call corrfdata(FAvgCorrf,FNormInt,MDBreakT,MDTrajT)
  close(1)
end if
call mpi_finalize(ierr)

end program YZTraj_TCF        


