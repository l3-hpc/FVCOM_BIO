
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
libbiomodule.a(grid.o)                       \
libbiomodule.a(mod_1D.o)                       \
libbiomodule.a(schism_glbl.o) \
libbiomodule.a(mod_parameter.o)                 \
libbiomodule.a(cgem.o) \
libbiomodule.a(mod_nutrient.o)                 \
libbiomodule.a(bio_mixing.o)                    \
libbiomodule.a(date_time.o) \
libbiomodule.a(grid_setup.o) \
libbiomodule.a(cgem_setup.o) \
libbiomodule.a(calcAgrow.o) \
libbiomodule.a(cgem_utils.o) \
libbiomodule.a(cgem_sink.o) \
libbiomodule.a(surface_flux.o) \
libbiomodule.a(calc_solar_zenith.o) \
libbiomodule.a(get_solar.o) \
libbiomodule.a(call_iop_par.o) \
libbiomodule.a(IOP_PARattenuation.o) \
libbiomodule.a(moc_src/singledouble.o) \
libbiomodule.a(moc_src/sw_adtg.o) \
libbiomodule.a(moc_src/sw_ptmp.o) \
libbiomodule.a(moc_src/sw_temp.o) \
libbiomodule.a(moc_src/p80.o) \
libbiomodule.a(moc_src/depth2press.o) \
libbiomodule.a(moc_src/f2pCO2.o) \
libbiomodule.a(moc_src/rho.o) \
libbiomodule.a(moc_src/constants.o) \
libbiomodule.a(moc_src/phsolvers.o) \
libbiomodule.a(moc_src/varsolver.o)\
libbiomodule.a(moc_src/vars.o) \
libbiomodule.a(moc_src/p2fCO2.o) \
libbiomodule.a(moc_src/gasx.o) \
libbiomodule.a(moc_src/rhoinsitu.o) \
libbiomodule.a(cgem_step.o) \
libbiomodule.a(bio_geo_chem.o)

LIBS	=	libbiomodule.a

FFLAGS_BIO = -O3 -fp-model strict
#FFLAGS_BIO = -g -debug -check all -traceback
#CPPFLAGS = -DDEBUG

libbiomodule.a: $(BIOMODULE)

.SUFFIXES: .o .f90 .F .F90 

.F.o:
	$(CPP) $(CPPFLAGS) $(CPPARGS) $*.F > $*.f90
#	$(FC) -c $(FFLAGS) $*.f90 > $*.o
	$(FC) -c $(FFLAGS_BIO) $*.f90 > $*.o
#	$(AR) $(ARFLGS) $@ $<	
	ar rv libbiomodule.a *.o

.F90.o:
	$(CPP) $(CPPFLAGS) $(CPPARGS) $*.F90 > $*.f90
	$(FC) -c $(FFLAGS_BIO) $*.f90 > $*.o
	ar rv libbiomodule.a *.o



clean:
	/bin/rm -f lib*.a  *.mod *.o *.f90

clobber:
		make clean
		/bin/rm -f *.f90

Makefiles::

includes::
include ../make.inc
# DO NOT DELETE
