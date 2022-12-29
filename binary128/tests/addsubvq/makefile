CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

all: taddv tsubv

addsubvq.o:../../vector/addsubvq.c ../../vector/addsubvq.h
	$(CC) $(COPT) -c $<

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $<

taddv.o:taddv.cpp ../../vector/addsubvq.h ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -I../../vector -c $< -o $@

tsubv.o:taddv.cpp ../../vector/addsubvq.h ../util/qmx2mpfr.h
	$(CPP) $(COPT) -I../util -I../../vector -c $< -o $@ -DSUBTRACT

taddv:taddv.o addsubvq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsubv:tsubv.o addsubvq.o qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

