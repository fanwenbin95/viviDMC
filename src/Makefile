EXE = viviDMC
POT = h2opes

FC = ifort
# -CB -g -traceback -debug parallel
FFLAGS = -CB -traceback -debug parallel -g \
	-I${MKLROOT}/include/intel64/lp64 \
	-I${MKLROOT}/include \
	-lmkl_lapack95_lp64 \
	-lstdc++ -O2
LDFLAGS = -lmkl_rt \
	${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
	-Wl,--start-group \
	${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
	${MKLROOT}/lib/intel64/libmkl_intel_thread.a \
	${MKLROOT}/lib/intel64/libmkl_core.a \
	-Wl,--end-group \
	-liomp5 -lpthread -lm -ldl \
	-L${MKLROOT}/lib/intel64 \
#	-lmbxlib -fopenmp -lfftw3 \
#	-L$(MBX_HOME)/install/lib/static \
#	-I$(MBX_HOME)/install/include

MODOBJ = module.o input.o init.o prpgt.o ana.o
OBJ = freq.o main.o \
interface.o math.o util.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<
%.o: %.f
	$(FC) $(FFLAGS) -c $<

deafult: $(MODOBJ) $(OBJ) $(POT).o
	$(FC) $(FFLAGS) $^ -o $(EXE) $(LDFLAGS)
#	./$(EXE)
#	mv fort.999 1.xyz

rm:
	rm -f *.o *.mod
	rm -f default $(EXE)
	rm -f fort*
	rm -f *~


