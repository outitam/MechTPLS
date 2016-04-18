TARGET	= ../TPLS

F90	= mpif90 

LDFLAGS = -O3 
FFLAGS  = -O3
LIBS    = -llapack 

SRC90	= tpls_constants.f90 \
	  tpls_mpi.f90 \
	  tpls_maths.f90 \
	  tpls_io.f90 \
	  tpls_userchk.f90 \
	  tpls_levelset.f90 \
	  tpls_configuration.f90 \
	  tpls_selective_frequency.f90 \
	  tpls_fluid.f90 \
	  tpls_pressure_solver.f90 \
	  jacobi_iteration_allflux.f90 \
	  tpls_helmholtz.f90 \
	  tpls_solver.f90 \
	  main.f90

OBJS90 	= $(SRC90:.f90=.o)

all:  $(TARGET)

$(TARGET): $(OBJS90)
	$(F90) $(LDFLAGS) -o $@ $(OBJS90) $(LIBS)
	@echo " ----------- Exectuable $@ cree ----------- "

.SUFFIXES:
.SUFFIXES: .f90 .o
.f90.o:
	$(F90) $(FFLAGS) $(INC) -c $<

clean:
	rm -f -R $(TARGET) $(OBJS90) *.out *.L *.mod $(REP)
