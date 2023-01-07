CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALL_RESULTS := tadd
ALL_RESULTS += tsub
ALL_RESULTS += tadd_man
ALL_RESULTS += tadd_libgcc
ALL_RESULTS += tsub_libgcc
ALL_RESULTS += tadd_man_libgcc

all: $(ALL_RESULTS)

addsubq.o:../../addsubq.c
	$(CC) $(COPT) -c $<

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $<

addsubq_with_all_rounding_modes.o:addsubq_with_all_rounding_modes.c
	$(CC) $(COPT) -c $<

tadd_man.o:tadd_man.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -I../util -c $<

tadd.o:tadd.cpp ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -c $< -o $@

tsub.o:tadd.cpp ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -c $< -o $@ -DSUBTRACT

tadd:tadd.o addsubq_with_all_rounding_modes.o addsubq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsub:tsub.o addsubq_with_all_rounding_modes.o addsubq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tadd_man:tadd_man.o addsubq_with_all_rounding_modes.o addsubq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@



tadd_libgcc:tadd.o addsubq_with_all_rounding_modes.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsub_libgcc:tsub.o addsubq_with_all_rounding_modes.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tadd_man_libgcc:tadd_man.o addsubq_with_all_rounding_modes.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@
