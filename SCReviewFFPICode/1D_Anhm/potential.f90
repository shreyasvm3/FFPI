!##############################################
!##############################################

!       This subroutine computes the potential V,
!       the gradient dV and the hessian d2V

!##############################################


! Last updated: Shreyas

subroutine VdV(x,V,dV)
use parameters, only : Mass, Ndof
implicit none

real*8, intent(in)      :: x(Ndof)
integer*8               :: i
real*8                  :: w1
real*8, intent(out)     :: V, dV(Ndof)

V       =       0.d0
dV      =       0.d0
w1      =       1.414

V       =       0.5d0*Mass(1)*w1**2*x(1)**2 - 0.1d0*x(1)**3 + 0.1d0*x(1)**4
dV(1)   =       Mass(1)*w1**2*x(1)          - 0.3d0*x(1)**2 + 0.4d0*x(1)**3

end subroutine VdV

!##############################################
 
