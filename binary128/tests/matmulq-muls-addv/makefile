CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native


ALLEXE :=
ALLEXE += mm_speed_ma
ALLEXE += mm_prec_ma

all: ${ALLEXE}

mulsq.o:../../vector/mulsq.c ../../vector/mulsq.h
	$(CC) $(COPT) -c $<

addsubvq.o:../../vector/addsubvq.c ../../vector/addsubvq.h
	$(CC) $(COPT) -c $<

matmulq.o:matmulq.c ../matmulq/matmulq.h
	$(CC) $(COPT) -I../matmulq -I../../vector -c $<

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $<

mm_speed.o:../matmulq/mm_speed.cpp ../matmulq/matmulq.h
	$(CPP) $(COPT) -c $<

mm_prec.o:../matmulq/mm_prec.cpp ../matmulq/matmulq.h ../util/qmx2mpfr.h
	$(CPP) $(COPT)  -I../util -c $<

mm_speed_ma:mm_speed.o matmulq.o addsubvq.o mulsq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=$@.map

mm_prec_ma:mm_prec.o qmx2mpfr.o matmulq.o addsubvq.o mulsq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=$@.map
