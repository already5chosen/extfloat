CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

all: tmuls

mulsq.o:../../vector/mulsq.c ../../vector/mulsq.h
	$(CC) $(COPT) -c $<

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $<

tmuls.o:tmuls.cpp ../../vector/mulsq.h ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -I../../vector -c $<

tmuls:tmuls.o mulsq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@
