## libxc?

# compiler selection, (un)comment as appropriate
### The Intel(R) fortran compiler (ifort)
LIBXC_DIR=/opt/libxc_ifort
ifeq ($(DEBUG),1)
  FC = ifort
  FCFLAGS = -g -CU -C -traceback -fpe0 -debug -openmp
  LDFLAGS = -openmp
else
  FC = ifort
  FCFLAGS = -O2 -openmp
  LDFLAGS = -O2 -openmp
endif

# ### The GNU fortran compiler (gfortran)
# # LIBXC_DIR=/opt/libxc_gfortran
# ifeq ($(DEBUG),1)
#  FC = gfortran
#  FCFLAGS = -O -g -fbounds-check -Wall -Wunused-parameter -ffpe-trap=invalid -fbacktrace -fdump-core -fopenmp
#  LDFLAGS = -fopenmp
# else
#  FC = gfortran
#  FCFLAGS = -O -fopenmp
#  LDFLAGS = -fopenmp
# endif

ifdef LIBXC_DIR
 INCLUDE_LIBXC=-I$(LIBXC_DIR)/include
 OBJECT_LIBXC=$(LIBXC_DIR)/lib/libxc.a
 FCFLAGS+=-DHAVE_LIBXC
endif

#### name, objects, libraries, includes
NAME=postg2
OBJS= promolmod.o param.o types.o tools.o tools_math.o tools_linpack.o reader.o wfnmod.o meshmod.o io.o sandbox.o postg2.o 
LIBS=$(OBJECT_LIBXC)
INCLUDE=$(INCLUDE_LIBXC)
####

BINS=$(NAME)
BINS_dbg=$(NAME)_dbg

%.o: %.F90
	$(FC) -c $(FCFLAGS) $(INCLUDE) -o $@ $<

%.o: %.f90
	$(FC) -c $(FCFLAGS) $(INCLUDE) -o $@ $<

%.o: %.f
	$(FC) -c $(FCFLAGS) $(INCLUDE) -o $@ $<

%.mod: %.o
	@if [ ! -f $@ ]; then rm $< ; $(MAKE) $< ; fi

# Targets
all: $(BINS)

debug: 
	DEBUG=1 $(MAKE) $(BINS_dbg)

clean:
	rm -f core *.mod *.o 

veryclean:
	rm -f core *.mod *.o

mrproper:
	rm -f core *.mod *.o $(BINS) $(BINS_dbg)

$(NAME): $(OBJS) $(LIBS)
	$(FC) -o $(NAME) $(LDFLAGS) $(OBJS) $(LIBS)

$(NAME)_dbg: $(OBJS) $(LIBS)
	$(FC) -o $(NAME)_dbg $(LDFLAGS) $(OBJS) $(LIBS)

# sandbox dependency
sandbox.o : sandbox.F90 sandbox/*.f90
	$(FC) -c $(FCFLAGS) $(INCLUDE) -o sandbox.o sandbox.F90

# Object dependencies
sandbox.o reader.o : io.mod
sandbox/energy.o sandbox/atomicb.o wfnmod.o postg2.o : meshmod.mod
sandbox.o wfnmod.o types.o tools_math.o tools.o reader.o promolmod.o postg2.o meshmod.o io.o : param.mod
postg2.o : promolmod.mod
postg2.o : reader.mod
sandbox.o wfnmod.o : tools.mod
tools_math.o meshmod.o : tools_linpack.mod
sandbox/atomicb.o promolmod.o postg2.o meshmod.o : tools_math.mod
sandbox/points.o sandbox/plane.o sandbox/line.o sandbox/cube_libxc.o sandbox/cube.o sandbox/atomicb.o sandbox.o wfnmod.o reader.o promolmod.o postg2.o meshmod.o io.o : types.mod
sandbox/atomicb.o sandbox.o postg2.o : wfnmod.mod
postg2.o : sandbox.mod
