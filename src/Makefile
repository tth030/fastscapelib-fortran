.SUFFIXES:.out .o .s .c .F .f .f90 .F90 .f95 .e .r .y .yr .ye .l .p .sh .csh .h

ifneq ($(findstring geodyn,$(HOSTNAME)),)
  F90 = mpif90
  F90_OPTIONS = -O -Mfreeform -Mdclchk -Mstandard -Mpreprocess
else ifneq ($(findstring pangraph,$(HOSTNAME)),)
  F90 = mpif90
  F90_OPTIONS = -O -cpp
  CC = gcc
else ifneq ($(findstring betzy,$(HOSTNAME)),)
  F90 = ifort
  F90_OPTIONS = -stand f08 -free -implicitnone -convert big_endian -check bounds -warn all -fpp
  CC = gcc
else
  F90 = gfortran
  CC = gcc-10
  F90_OPTIONS = -ffree-form -ffree-line-length-none -fimplicit-none -fbounds-check -Wall -Wno-unused-dummy-argument -cpp -fall-intrinsics -fallow-argument-mismatch -fconvert=big-endian # needed to remove -std=f2008 -ffree-line-length-none is required for long macros
endif

ON_WINDOWS = f

FASTSCAPE_OBJ = \
FastScapeErrorCodes.o \
NumRecipes.o \
Qshep2d_mod.o \
FastScape_ctx.o \
FastScape_api.o \
StreamPowerLaw.o \
Diffusion.o \
LocalMinima.o \
Advect.o \
Advect3d.o \
Advect3d_2.o \
Advect3d_3.o \
Marine.o \
VTK.o \
FlowRouting.o \
Strati.o \
Uplift.o \
TerrainDerivatives.o

.f90.o:
	$(F90) -c $(F90_OPTIONS) $*.f90 -L/usr/lib -llapack -L/usr/lib -lblas

all: $(FASTSCAPE_OBJ) LIB
	make LIB
	make RUNF90
	#make RUNC

LIB: $(FASTSCAPE_OBJ)
	ar rc libFastScape.a $(FASTSCAPE_OBJ)

RUNF90: LIB
	$(F90) -O $(F90_OPTIONS) FastScape_run.f90 -o FastScape_run_f90 -L. -lFastScape -L/usr/lib -llapack -L/usr/lib -lblas

RUNC: LIB
	$(CC) -O FastScape_run.c -o FastScape_run_c -L. -lFastScape -lgfortran -L/usr/lib -llapack -L/usr/lib -lblas

clean:
	rm *.o *.mod *.a
