bindir = ../bin/
srcdir = ./src
FC = mpif90 -O2 
LIB = -llapack -lblas
.SUFFIXES: .f90 .o

#--------------------------------------------------

MOD = parameters.f90 \
      potential.f90 \
      gen_MD_traj.f90 \
      estimator.f90 \
      poststat.f90

SRC = main.f90

OBJ = $(MOD) $(SRC)

# rule to make a .o file from a .f90 file
.f90.o:
	$(FC) -c $<
.f.o:
	$(FC) -c $<

dyn.x :  $(OBJ)
	 $(FC) $(OPTION) -o $@ $(OBJ) $(LIB)

clear   :
	rm *.x
	rm *.o
	rm *.mod
