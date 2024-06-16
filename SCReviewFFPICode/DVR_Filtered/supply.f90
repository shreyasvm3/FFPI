!########################################
!       Potential and its gradient
!########################################
subroutine VofR(R,VR,dV)
use parameters, only : Mass
implicit none

real*8, intent(in)      :: R
real*8, intent(out)     :: VR, dV
real*8                  :: k

k   = 2.d0
VR  = 0.5d0*k*R**2-0.1d0*R**3+0.1d0*R**4
dV  = k*R-0.3d0*R**2+0.4d0*R**3

end subroutine
!########################################
!       Split loop over procs
!########################################
subroutine para_range(n1,n2,nprocs,irank,ista,iend)
implicit none
integer, intent(in) :: n1,n2,nprocs,irank
integer, intent(out) :: ista,iend
integer :: iwork1,iwork2
iwork1 = (n2 - n1 + 1)/nprocs
iwork2 = mod(n2 - n1 + 1, nprocs)
ista = irank*iwork1 + n1 + min(irank, iwork2)
iend = ista + iwork1 -1
if (iwork2 .gt. irank) then
        iend = iend + 1
endif
return
end subroutine para_range
!########################################
!       A matrix in yz basis
!########################################
subroutine Amatrix(Amat)
use parameters, only : NxPM, GridxP, GridxM, Width, qIn, NR

real*8,intent(out)      :: Amat(NXpm)
integer                 :: i, j

Amat = 0.d0

do i = 1, NXpm
      Amat(i) = dexp(-Width*((GridxP(i)-qIn)**2+(GridxM(i)-qIn)**2)/2.d0)
end do

!Amat = Amat/dsqrt(dble(NR))

end subroutine Amatrix
!########################################
!       B matrix in yz basis
!########################################
subroutine Bmatrix(Bmat)
use parameters, only : NXpm, GridxP, GridxM, NR, GridR

real*8,intent(out)      :: Bmat(NXpm)
integer                 :: i

Bmat = 0.d0

do i = 1, NXpm
   if(GridxP(i).eq.GridxM(i)) then
      Bmat(i) = GridxP(i)
   end if
end do

!Bmat = Bmat/dsqrt(dble(NR))

end subroutine Bmatrix
!##############################################################
!       Calculates the fwd and bck prop for time Delta(t)
!###############################################################
subroutine STPropXbasis(H,Energy,FpropX,BpropX)
use parameters, only : NR, Iu, Eps
implicit none

integer                 :: i
real*8,intent(in)       :: H(NR,NR), Energy(NR)
complex*16              :: temp(NR,NR)
real*8                  :: Ht(NR,NR)
complex*16              :: FpropX(NR,NR), BpropX(NR,NR)

FpropX = 0.d0
BpropX = 0.d0
temp   = 0.d0 

do i = 1, NR
   FpropX(i,i) = cdexp(-Iu*Energy(i)*Eps/2.d0)
end do

Ht = transpose(H)

temp      = matmul(H,FpropX)
FpropX    = matmul(temp,Ht)

do i = 1, NR
   BpropX(i,i) = cdexp(Iu*Energy(i)*Eps/2.d0)
end do

temp      = matmul(H,BpropX)
BpropX    = matmul(temp,Ht)

end subroutine STPropXbasis
!##############################################################
!       Calculates the fwd and bck prop for time Delta(t)
!###############################################################
subroutine STPropPMbasis(FpropX,BpropX,PropPM)
use parameters, only : NR, NXpm, Pidx, Midx
implicit none

complex*16, intent(in) :: FpropX(NR,NR), BpropX(NR,NR)
integer                :: i, j
complex*16,intent(out) :: PropPM(NXpm,NXpm)

PropPM = 0.d0

do i = 1, NXpm
   do j = 1, NXpm
      PropPM(i,j) = FpropX(Pidx(i),Pidx(j))*BpropX(Midx(j),Midx(i))
   end do
end do

end subroutine STPropPMbasis
!##############################################################
!       Filtered prop of time length eps
!###############################################################
subroutine FilteredProp(PropPM,FilProp)
use parameters, only : NXpm, cFil, dFil, GridY, GridZ
implicit none

complex*16, intent(in)  :: PropPM(NXpm,NXpm)
integer                 :: i, j, k
real*8                  :: fj
complex*16, intent(out) :: FilProp(NXpm,NXpm)

FilProp = 0.d0

!{Manually multiply 2 stps with filinov filters incorporated into each element
do i = 1, NXpm
   do j = 1, NXpm
      do k = 1, NXpm
         call fjCalc(GridY(i),GridY(k),GridY(j),fj)
         FilProp(i,j)=FilProp(i,j)+PropPM(i,k)*PropPM(k,j)*dexp(-5.d-1*(cFil*GridZ(k)**2+dFil*fj**2))
      end do
   end do  
end do
!}
Filprop = Filprop*(1.d0+(cFil*dFil))**(0.5d0) ! Fil prefactor

end subroutine FilteredProp
!##############################################################
! Calculate the correlation function from all the matrices
!###############################################################
subroutine cf0time(Amat,Bmat,cf0t)
use parameters, only : NXpm
implicit none

real*8, intent(in)      :: Amat(NXpm), Bmat(NXpm)
complex*16,intent(out)  :: cf0t

cf0t = dot_product(Amat,Bmat)

end subroutine cf0time
!##############################################################
!       Calculates the fwd and bck prop for full time t
!###############################################################
subroutine Propst(ist,FilPropdt,FilPropst)
use parameters, only : NXpm
implicit none

integer,intent(in)      :: ist
complex*16,intent(in)   :: FilPropdt(NXpm,NXpm)
integer                 :: i
complex*16,intent(out)  :: FilPropst(NXpm,NXpm)

FilPropst    = 0.d0

do i = 1, NXpm
   FilPropst(i,i) = 1.d0
end do
!}

!{Perform istart matmuls, each with an Eps prop
do i = 0, ist-1
   FilPropst = matmul(FilPropst,FilPropdt)
end do
!}
end subroutine Propst
!
!##############################################################
!       Calculates the fwd and bck prop for full time t
!###############################################################
subroutine Fullprop(FilProp,Prop)
use parameters, only : NXpm
implicit none

complex*16,intent(in)   :: FilProp(NXpm,NXpm)
complex*16,intent(inout):: Prop(NXpm,NXpm)

Prop = matmul(Prop,FilProp)

end subroutine Fullprop
!##############################################################
! Calculate fj from yj, yj+1, yj-1 
!###############################################################
subroutine fjCalc(yjp,yj,yjm,fj)
use parameters, only : Mass, Eps
implicit none

real*8, intent(in)      :: yjp, yj, yjm
real*8, intent(out)     :: fj
real*8                  :: Vj, dV, dt

dt = Eps/2.d0

call VofR(yj,Vj,dV)
fj = Mass*(yjp-2.d0*yj+yjm)/dt**2 + dV

end subroutine

!##############################################################
! Calculate the correlation function from all the matrices
!###############################################################
subroutine Corrf(Amat,Bmat,PropPM,cf)
use parameters, only : NXpm, Eps, Nts, cFil, dFil
implicit none

real*8, intent(in)      :: Amat(NXpm), Bmat(NXpm)
complex*16,intent(in)   :: PropPM(NXpm,NXpm)
complex*16              :: temp(NXpm)
real*8                  :: filpre
integer                 :: i
complex*16,intent(out)  :: cf

cf      = 0.d0

temp    = matmul(Amat,PropPM)
cf      = dot_product(temp,Bmat)

end subroutine Corrf 

