CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALL_RESULTS := tmul
ALL_RESULTS += tmul_man
ALL_RESULTS += tmul_libgcc
ALL_RESULTS += tmul_man_libgcc

all: $(ALL_RESULTS)

mulq.o:../../mulq_rm.c
	$(CC) $(COPT) -c $< -o $@

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $<

mulq_with_all_rounding_modes.o:mulq_with_all_rounding_modes.c
	$(CC) $(COPT) -c $<

tmul_man.o:tmul_man.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -I../util -c $<

tmul.o:tmul.cpp ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -c $< -o $@

tmul:tmul.o mulq_with_all_rounding_modes.o mulq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man:tmul_man.o mulq_with_all_rounding_modes.o mulq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@



tmul_libgcc:tmul.o mulq_with_all_rounding_modes.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man_libgcc:tmul_man.o mulq_with_all_rounding_modes.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@
