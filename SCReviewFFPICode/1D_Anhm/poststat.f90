!###################################################################
!       This module contains subroutines that get averages and std.
!       dev. of data
!###################################################################

module stat
implicit none

contains

!###################################################################
!        Corrf. Avg.
!###################################################################
subroutine corrfdata(FAvgCorrf,FNormInt,MDBreakC,MDTrajC)
use parameters, only : Eps, Nts
implicit none

complex*16,intent(inout) :: FAvgCorrf(3,0:Nts), FNormInt(3,0:Nts)
integer,intent(in)       :: MDBreakC, MDTrajC
complex*16               :: FAvg(3,0:Nts), mean
integer                  :: i, j, t
real*8                   :: sdr, sdi

print*, "MD Broken Trajs.", MDBreakC, MDTrajC, MDBreakC/dble(MDTrajC)*1.d2

do t = 0, Nts
   do i = 1, 3
      FAvg(i,t) = FAvgCorrf(i,t)*dconjg(FNormInt(i,t))/cdabs(FNormInt(i,t))**2
   end do
   mean = sum(Favg(:,t))/3.d0
   call stddev(3,real(Favg(:,t)),real(mean),sdr)   
   call stddev(3,aimag(Favg(:,t)),aimag(mean),sdi)
   write(1,'(E15.6)',advance='no'), t*Eps
   write(1,'(E15.6)',advance='no'), real(mean)
   write(1,'(E15.6)',advance='no'), aimag(mean)
   write(1,'(E15.6)',advance='no'), sdr
   write(1,'(E15.6)'), sdi
end do

end subroutine
!###################################################################
!       Standard Deviation
!###################################################################
subroutine stddev(n,x,mean,sd)
implicit none

integer,intent(in)      :: n
real*8,intent(in)       :: mean, x(n)
real*8,intent(out)      :: sd
integer                 :: i,j

sd = 0.d0

do j = 1, n
        sd = sd + (x(j)-mean)**2
end do
sd = dsqrt(sd/(n-1))

end subroutine

end module stat    
