bindir = ../bin/
srcdir = ./src
FC = ifort -O2
FC2= mpif90 -O2
LIB = -llapack -lblas
OPTION = -mcmodel=large
.SUFFIXES: .f90 .o

#--------------------------------------------------
MOD = supply.f90 \
      params.o

SRC  = hamiltonian.f90
SRC2 = timecorr.f90


OBJ  = $(MOD) $(SRC)
OBJ2 = $(MOD) $(SRC2)

# rule to make a .o file from a .f90 file
 .f90.o:
	$(FC) -c $<
 .f.o:
	$(FC) -c $<

dvr.x : $(OBJ) 
	$(FC) $(OPTION) -o $@ $(OBJ) $(LIB)

dyn.x : $(OBJ2)
	$(FC2) $(OPTION) -o $@ $(OBJ2) $(LIB)

all : dyn.x dvr.x 

clear	:
	rm *.o
	rm *.mod
     
