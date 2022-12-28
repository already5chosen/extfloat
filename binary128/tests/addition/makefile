CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native -I../c

all: tadd tsub tadd_man

addsubq.o:../../addsubq.c
	$(CC) $(COPT) -c $<

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $<

tadd_man.o:tadd_man.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -I../util -c $<

tadd.o:tadd.cpp ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -c $< -o $@

tsub.o:tadd.cpp ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -c $< -o $@ -DSUBTRACT

tadd:tadd.o addsubq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsub:tsub.o addsubq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tadd_man:tadd_man.o addsubq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@
