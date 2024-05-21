
# Master Makefile for making the executable of the generalized biological model.

#.SUFFIXES:
#.SUFFIXES: .F
#.SUFFIXES: .o .f90 .F .F90 
	SHELL	= /bin/sh
#	DEF_FLAGS     = -P -C -traditional 
#	NETCDFINCDIR	= /usr/local/include
# 	ETCDFLIB       = /usr/local/lib -lnetcdf

#         CPP      = /usr/bin/cpp 
#         CPPFLAGS = $(DEF_FLAGS) #-DINTEL   
#         FC    = ifort -checkall
#         FC       = mpif90 -static-libcxa -i-static

BIOMODULE	= \
libbiomodule.a(mod_1D.o)                       \
libbiomodule.a(eut.o)\
libbiomodule.a(mod_parameter.o)                 \
libbiomodule.a(bio_mixing.o)                    \
libbiomodule.a(states.o)\
libbiomodule.a(flags.o)\
libbiomodule.a(Read_InputFile_GD.o)\
libbiomodule.a(gd_light.o)\
libbiomodule.a(sinking.o)\
libbiomodule.a(bio_geo_chem.o)\
libbiomodule.a(zoo.o)\
libbiomodule.a(diatoms.o)\
libbiomodule.a(greens.o)\
libbiomodule.a(carbon.o)\
libbiomodule.a(phosph.o)\
libbiomodule.a(silica.o)\
libbiomodule.a(nitrog.o)\
libbiomodule.a(dissolved_oxygen.o)

LIBS	=	libbiomodule.a

FFLAGS_BIO = -O3 
#FFLAGS_BIO = -g -debug -check all -traceback

libbiomodule.a: $(BIOMODULE)

.SUFFIXES: .o .f90 .F .F90 

.F.o:
	$(CPP) $(CPPFLAGS) $(CPPARGS) $*.F > $*.f90
#	$(FC) -c $(FFLAGS) $*.f90 > $*.o
	$(FC) -c $(FFLAGS_BIO) $*.f90 > $*.o

#	$(AR) $(ARFLGS) $@ $<	
	ar rv libbiomodule.a *.o
	\rm $*.f90

clean:
	/bin/rm -f lib*.a  *.mod *.o *.f90

clobber:
		make clean
		/bin/rm -f *.f90

Makefiles::

includes::
include ../make.inc
# DO NOT DELETE
#
#
#
#
