!###############################################################################
!###############################################################################

!      This module contains subroutines that generate inital coordinates 
!      and momenta, and propogate them as classical trajectories 

!      Also calculates energy conservation: (EInt-Efin)/EInt
!      Calls a traj "broken" if EngError > 1d-5
!###############################################################################

! Last updated: Shreyas

module MD
implicit none
contains
!!###############################################################################
!       Generate Initial y & z, py 
!!###############################################################################
subroutine sampleIC(y0,py,z0,pz)
use parameters, only : Ndof, qIn, pIn, Width0, AlphaY, AlphaZ
implicit none

real*8, intent(out)     :: y0(Ndof), pY(Ndof), z0(Ndof), pZ(Ndof)
integer                 :: i, j

y0      = 0.d0
z0      = 0.d0
pY      = 0.d0
pZ      = 0.d0

do i = 1, Ndof
   call gauss(dsqrt(0.5d0/Width0(i)),y0(i))
   y0(i) = y0(i) + qIn(i) 
   call gauss(dsqrt(2.d0*Width0(i)*AlphaY(i)),pY(i))
   pY(i) = pY(i) + pIn(i)
   call gauss(dsqrt(2.d0/Width0(i)),z0(i))
   call gauss(dsqrt(0.5d0*Width0(i)*AlphaZ(i)),pZ(i))
end do

end subroutine sampleIC
!###############################################################################
!       Propagate YZ trajectories
!###############################################################################
subroutine PropYZ(y0,pY0,z0,pZ0,yTraj,zTraj,MDBreakC,MDTrajC,MDFlag)
use parameters, only : Ndof, Mass, Eps, Nts
implicit none

real*8, intent(in)      :: y0(Ndof), pY0(Ndof), z0(Ndof), pZ0(Ndof)
real*8,dimension(Ndof)  :: y, py, z, pZ, mY, mZ
real*8                  :: Vp, Vm ,dVp(Ndof), dVm(Ndof), Ei, Ef, dE
integer                 :: i, j
real*8,intent(out)      :: yTraj(Ndof,0:Nts), zTraj(Ndof,0:Ndof)
integer, intent(inout)  :: MDBreakC, MDTrajC
integer, intent(out)    :: MDFlag

MDFlag  = 0
y       = y0
py      = 2.d0*pY0           !The momentum sampled is m.dy/dt 
z       = z0
pz      = pZ0/2.d0           !The momentum sampled is m.dz/dt
mY      = 2.d0*Mass
mZ      = Mass/2.d0
      
call VdV(y+z/2.d0,Vp,dVp)
call VdV(y-z/2.d0,Vm,dVm)

Ei = Vp + Vm
do i = 1, Ndof
   Ei = Ei + 0.5d0*(py(i)**2/mY(i)+pz(i)**2/mZ(i))
end do

do j = 1, Nts     ! Loop over timesteps
   ! Propagate for one timestep
   call Integrator(y,py,z,pz)
   ! Check energy conservation
   call VdV(y+z/2.d0,Vp,dVp)
   call VdV(y-z/2.d0,Vm,dVm)
   Ef = Vp + Vm
   do i = 1, Ndof
      Ef = Ef + 0.5d0*(py(i)**2/mY(i)+pz(i)**2/mZ(i))
   enddo
  ! Save phase space coordinates at each timestep
   yTraj(:,j) = y
   zTraj(:,j) = z
   if (dabs(1.d0-Ef/Ei).ge.1.d-5) then
     MDflag = 1
     MDBreakC = MDBreakC + 1
     goto 112
   endif
enddo

112 continue
MDTrajC = MDTrajC + 1

end subroutine PropYZ
!######################################################!
!               Symplectic Integrator                  !       
!######################################################!
subroutine integrator(y,py,z,pz)
use parameters, only : Mass, Ndof, Eps
implicit none

real*8, intent(inout) :: py(Ndof), y(Ndof), z(Ndof), pz(Ndof)
integer               :: j, i, k
real*8                :: Vp, Vm, const, mY(Ndof), mZ(Ndof)
real*8                :: dVp(Ndof), dVm(Ndof), a(4), b(4)
real*8                :: tempY(Ndof), tempZ(Ndof), dt

dt    = Eps
mY    = 2.d0*Mass
mZ    = Mass/2.d0
const = 1.d0/dsqrt(3.d0)
a(1)  = 0.5d0*(1.d0-const)*dt
a(2)  =  const*dt
a(3)  = -const*dt
a(4)  = 0.5d0*(1.d0+const)*dt
b(1)  = 0.d0
b(2)  = 0.5d0*(0.5d0+const)*dt
b(3)  = 0.5d0*dt
b(4)  = 0.5d0*(0.5d0-const)*dt

do j = 1, 4
   if (j.gt.1) then
      call VdV(y+z/2.d0,Vp,dVp)
      call VdV(y-z/2.d0,Vm,dVm)
      py  = py - b(j)*(dVp+dVm)
      pz  = pz - b(j)*(dVp-dVm)/2.d0
   end if
   do i = 1, Ndof
      tempY(i) = pY(i)/mY(i)
      tempZ(i) = pZ(i)/mZ(i)
   enddo
   y   = y + a(j)*tempY
   z   = z + a(j)*tempZ
end do

end subroutine integrator
!############################################################
!       Inititialize Random Seed (taken from Matt's code)
!############################################################
subroutine init_random_seed()
implicit none 
       
integer                 :: i, n, clock
integer, allocatable    :: seed(:)

call random_seed(size=n)
allocate(seed(n))
        
call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /) 

call random_seed(put=seed)

deallocate(seed)

end subroutine init_random_seed
!#########################################################

!#########################################################
!      Gaussian Random number generator
!      Taken from Matt's code.
!      Based on the polar form of the Box-Muller algorithm.
!#########################################################

subroutine gauss(sigma,g)

implicit none
real*8, intent(in)      :: sigma
real*8                  :: w, v1, v2, l
real*8                  :: s1, s2
real*8,intent(out)      :: g
w = 2.d0

do
call random_number(s1)
call random_number(s2)
v1 = 2.d0*s1 - 1.d0
v2 = 2.d0*s2 - 1.d0
w = v1*v1 + v2*v2
if (w.lt.1.d0) exit
end do

l       = v1*sqrt(-2.d0*log(w)/(w))
g       = sigma*l
end subroutine gauss

!######################################################
!######################################################!
!             Compute MPI iteration range              !       
!              This is Matt's subroutine
!######################################################!
subroutine para_range(n1,n2,nprocs,irank,ista,iend)
implicit none
integer, intent(in)  :: n1, n2
integer,intent(in)   :: nprocs, irank
integer, intent(out) :: ista, iend
integer              :: iwork1, iwork2

iwork1 = (n2 - n1 + 1)/nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista   = irank*iwork1 + n1 + min(irank, iwork2)
iend   = ista + iwork1 -1

if (iwork2 .gt. irank) then
        iend = iend + 1
endif

return

end subroutine para_range
!######################################################!

end module MD
