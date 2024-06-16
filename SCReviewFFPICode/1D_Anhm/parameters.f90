!####################################################################
!       This module defines global variables and the subroutine 
!       input reads the input file to assign values to these vatriables.
!####################################################################
module parameters
implicit none

! Define pi and the imaginary unit
real*8, parameter     :: pi=dacos(-1.d0)
complex*16, parameter :: II=(0.d0,1.d0)

!Degrees of freedom
integer                 :: Ndof, Nts

!Epsilon time slice, and total time
real*8                  :: Eps, Time

!Array for the mass
real*8,allocatable      :: Mass(:), Cfil(:), Dfil(:)

!Initial State parameters: 
real*8,allocatable      :: qIn(:), pIn(:), Width0(:), AlphaY(:), AlphaZ(:)

!# of initial conditions.
integer                 :: NTraj

contains

!####################################################################
!       Input subroutine
!####################################################################
subroutine input()
implicit none

integer         :: i
character*75    :: infostr

open(555,file='input.txt',status='old')

read(555,'(a75)') infostr
read(555,*) Ndof

allocate(Mass(Ndof),qIn(Ndof),pIn(Ndof),Width0(Ndof))
allocate(AlphaY(Ndof),AlphaZ(Ndof),Cfil(Ndof),Dfil(Ndof))
read(555,'(a75)') infostr
read(555,*) Nts

read(555,'(a75)') infostr
read(555,*) Time

Eps = Time/dble(Nts)

read(555,'(a75)') infostr
do i = 1, Ndof
   read(555,*)     Mass(i)
end do

read(555,'(a75)') infostr
do i = 1, Ndof
   read(555,*) qIn(i), pIn(i), Width0(i), AlphaY(i), AlphaZ(i), Cfil(i), Dfil(i)
end do

read(555,'(a75)') infostr
read(555,*) NTraj

end subroutine

end module parameters




