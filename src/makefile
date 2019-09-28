#****************************************************************************
# C-K SPHERE makefile
# Pablo Daniel Mininni - 28/11/2005
#****************************************************************************

COMP     = gfortran -c
LINK     = gfortran -O4
LMPI     = -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi
IDIR     = -I/opt/openmpi-4.0.1/include -I/opt/openmpi-4.0.1/lib
LDIR     = -L/opt/openmpi-4.0.1/lib -Wl,-flat_namespace -Wl,-commons,use_dylibs

module.o:
	$(COMP) $(IDIR) module.f90

tables.o:
	$(COMP) bessel.f90
	$(COMP) legendre.f90
	$(COMP) integral.f90
	$(COMP) sckfunctions.f90
	$(COMP) tables.f90

sphere_aux.o:
	$(COMP) sphere_aux.f90

parallel.o:
	$(COMP) $(IDIR) parallel.f90

create:
	@$(MAKE) module.o
	@$(MAKE) tables.o
	@$(MAKE) parallel.o
	$(COMP) create_tables.f90
	$(LINK) bessel.o integral.o legendre.o module.o tables.o \
         parallel.o create_tables.o $(LMPI) $(LDIR) -o ../bin/create

test:
	@$(MAKE) module.o
	@$(MAKE) tables.o
	@$(MAKE) parallel.o
	$(COMP) test.f90
	$(LINK) bessel.o integral.o legendre.o module.o tables.o \
         parallel.o test.o $(LMPI) $(LDIR) -o ../bin/test

shd:
	@$(MAKE) module.o
	@$(MAKE) tables.o
	@$(MAKE) sphere_aux.o
	@$(MAKE) parallel.o
	$(COMP) sphere_hd.f90
	$(LINK) module.o legendre.o bessel.o sphere_aux.o sckfunctions.o \
         parallel.o sphere_hd.o $(LMPI) $(LDIR) -o ../bin/shd

smhd:
	@$(MAKE) module.o
	@$(MAKE) tables.o
	@$(MAKE) sphere_aux.o
	@$(MAKE) parallel.o
	$(COMP) sphere_mhd.f90
	$(LINK) module.o legendre.o bessel.o sphere_aux.o sckfunctions.o \
         parallel.o sphere_mhd.o $(LMPI) $(LDIR) -o ../bin/smhd

clean:
	rm -f *.o *.d *.s *.int *.inf *.mod work.pc* ifc*

dist:
	@$(MAKE) clean
	rm -f ../bin/create ../bin/shd ../bin/smhd
