!###############################################################

!       This module contains subroutines that are required to 
!       calculate the estimator at each samping point.

!       Check is the action is defined correctly
!       Check how to deal with avg over parallelization

!###############################################################

module Estimator
implicit none
contains

!###############################################################
!       Calculate the estimator
!###############################################################
subroutine Est(yCoord,zCoord,AvgCorrf,NormInt)
use parameters, only : II, Ndof, Nts
implicit none

real*8,intent(in)               :: yCoord(Ndof,0:Nts), zCoord(Ndof,0:Nts)
complex*16                      :: Corrf(0:Nts), NormC(0:Nts)
complex*16,intent(inout)        :: AvgCorrf(0:Nts), NormInt(0:Nts)  
real*8                          :: i, AMat, BMat, Sfb(0:Nts), Fexp(0:Nts)
integer                         :: t

call Act(yCoord,zCoord,Sfb)
call FilExp(yCoord,zCoord,FExp)
call A_Matrix(AMat)

do t = 0, Nts
   call B_Matrix(yCoord(1,t),BMat)
   Corrf(t)     = AMat*BMat*Fexp(t)*cdexp(Sfb(t)*II)
   NormC(t)     = Fexp(t)*cdexp(Sfb(t)*II)
   AvgCorrf(t)  = AvgCorrf(t) + Corrf(t)
   NormInt(t)   = NormInt(t)  + NormC(t)
end do

end subroutine

!###############################################################
!       A matrix
!###############################################################
subroutine A_Matrix(AMat)
implicit none

real*8,intent(out) :: AMat

AMat = 1.d0

end subroutine

!###############################################################
!       B Matrix
!###############################################################
subroutine B_Matrix(YIn,YOut)
implicit none
real*8,intent(in)       :: yIn
real*8,intent(out)      :: YOut

YOut = YIn

end subroutine

!###############################################################
!       Action
!###############################################################
subroutine Act(yCoord,zCoord,Sfb)
use parameters, only : Ndof, Nts, Mass, Eps

real*8,intent(in)       :: yCoord(Ndof,0:Nts), zCoord(Ndof,0:Nts)
real*8,intent(out)      :: Sfb(0:Nts)
integer                 :: i, j, t
real*8                  :: V, dV(Ndof), VPlus, VMinus

Sfb = 0.d0

! The (y1-y0)z0 term that needs to be included for all t.
! Second square bracket in Eq.(5)
do i = 1, Ndof
   call VdV(yCoord(:,0)+zCoord(:,0)/2.d0,VPlus,dV)
   call VdV(yCoord(:,0)-zCoord(:,0)/2.d0,VMinus,dV)
   Sfb(1) = Sfb(1) - Mass(i)/Eps*(yCoord(i,1)-yCoord(i,0))*zCoord(i,0)
end do
Sfb(1) = Sfb(1) - Eps*(VPlus-VMinus)/2.d0

!First square bracket in Eq. (5)        
do j = 1, Nts-1
   call VdV(yCoord(:,j)+zCoord(:,j)/2.d0,VPlus,dV)
   call VdV(yCoord(:,j)-zCoord(:,j)/2.d0,VMinus,dV)
   Sfb(j+1) = Sfb(j)
   do i = 1,Ndof
      Sfb(j+1) = Sfb(j+1) - Mass(i)/Eps*zCoord(i,j)*(yCoord(i,j+1)-2.d0*yCoord(i,j)+yCoord(i,j-1)) 
   end do
   Sfb(j+1) = Sfb(j+1) - Eps*(VPlus-VMinus)
end do
!Last square bracket in eq.(5) should be zero because z(t) = 0.d0

end subroutine
!###############################################################
!###############################################################
!       Filinov exponent
!###############################################################
subroutine FilExp(yTraj,zTraj,FExp)
use parameters, only : Ndof, Nts, cFil, dFil
implicit none

real*8,intent(in)       :: yTraj(Ndof,0:Nts), zTraj(Ndof,0:Nts)
integer                 :: i, j
real*8                  :: fj(Ndof)
real*8,intent(out)      :: Fexp(0:Nts)

Fexp = 1.d0

do j = 1, Nts-1
   call fjCalc(yTraj(:,j+1),yTraj(:,j),yTraj(:,j-1),fj)
   Fexp(j+1) = Fexp(j)
   do i = 1, Ndof
      Fexp(j+1) = Fexp(j+1)*dexp(-0.5d0*(cFil(i)*zTraj(i,j)**2+dFil(i)*fj(i)**2))
   end do
end do

end subroutine FilExp
!###############################################################################
!       Subroutine that calculates fj, given y.
!###############################################################################
subroutine fjCalc(yjp1,yj,yjm1,fj)
use parameters, only : Ndof, Mass, Eps

real*8, intent(in)      :: yjp1(Ndof), yj(Ndof), yjm1(Ndof)
real*8                  :: V, dV(Ndof)
integer                 :: i
real*8, intent(out)     :: fj(Ndof)

call VdV(yj,V,dV)
do i = 1, Ndof
   fj(i) = Mass(i)/Eps*(yjp1(i)-2.d0*yj(i)+yjm1(i))+Eps*dV(i)
end do

end subroutine fjCalc
!###############################################################################

end module Estimator
