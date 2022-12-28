CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native -I../c

all: tmul tmul_man

mulq.o:../../mulq.c
	$(CC) $(COPT) -c $<

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $<

tmul_man.o:tmul_man.c  ../util/qmx2mpfr.h
	$(CC) $(COPT) -I../util -c $<

tmul.o:tmul.cpp ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -c $<

tmul:tmul.o mulq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man:tmul_man.o mulq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@
