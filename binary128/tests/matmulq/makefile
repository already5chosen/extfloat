CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native


ALLEXE :=
ALLEXE += mm_speed
ALLEXE += mm_speed_ma
ALLEXE += mm_prec
ALLEXE += mm_prec_ma

all: ${ALLEXE}


mulq.o:../../mulq.c
	$(CC) $(COPT) -c $< -o $@

addsubq.o:../../addsubq.c
	$(CC) $(COPT) -c $<

matmulq.o:matmulq.c matmulq.h
	$(CC) $(COPT) -c $<

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $<

mm_speed.o:mm_speed.cpp matmulq.h
	$(CPP) $(COPT) -c $<

mm_prec.o:mm_prec.cpp matmulq.h ../util/qmx2mpfr.h
	$(CPP) $(COPT)  -I../util -c $<

mm_speed:mm_speed.o matmulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=$@.map

mm_speed_ma:mm_speed.o matmulq.o addsubq.o mulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=$@.map

mm_prec:mm_prec.o qmx2mpfr.o matmulq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=$@.map

mm_prec_ma:mm_prec.o qmx2mpfr.o matmulq.o addsubq.o mulq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=$@.map
