include Makefile.inc
CODE = thermal_conductivity
PROG = ../../build/$(CODE)/$(CODE)
OBJECT_PATH=../../build/$(CODE)/

OBJS = \
$(OBJECT_PATH)main.o\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)scattering.o\
$(OBJECT_PATH)kappa.o\
$(OBJECT_PATH)cumulative_kappa.o

LPATH = -L../../lib $(blaslapackLPATH) $(incLPATHmpi) $(incLPATHhdf)
IPATH = -I../../inc/libolle -I../../inc/libflap $(blaslapackIPATH) $(incIPATHmpi) $(incIPATHhdf)
LIBS = ../../lib/libolle.a -lflap $(blaslapackLIBS) $(incLIBSmpi) $(incLIBShdf)

#OPT = -O0 -fbacktrace -fcheck=all -finit-real=nan -finit-derived
F90 = $(FC) $(LPATH) $(IPATH) $(MODULE_FLAG) $(OBJECT_PATH) #$(warnings_gcc)
F90FLAGS = $(OPT) $(MODS) $(LIBS)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) $(OBJECT_PATH)*.mod

$(OBJECT_PATH)main.o: \
main.f90\
$(OBJECT_PATH)options.o\
$(OBJECT_PATH)scattering.o\
$(OBJECT_PATH)kappa.o\
$(OBJECT_PATH)cumulative_kappa.o
	$(F90) $(OPT) $(F90FLAGS) -c $< $(LIBS) -o $@
$(OBJECT_PATH)scattering.o: scattering.f90
	$(F90) $(OPT) $(F90FLAGS) -c $< $(LIBS) -o $@
$(OBJECT_PATH)kappa.o:\
	kappa.f90\
	$(OBJECT_PATH)scattering.o
	$(F90) $(OPT) $(F90FLAGS) -c $< $(LIBS) -o $@
$(OBJECT_PATH)cumulative_kappa.o: cumulative_kappa.f90\
	$(OBJECT_PATH)kappa.o
	$(F90) $(OPT) $(F90FLAGS) -c $< $(LIBS) -o $@
$(OBJECT_PATH)options.o: options.f90
	$(F90) $(OPT) $(F90FLAGS) -c $< $(LIBS) -o $@
