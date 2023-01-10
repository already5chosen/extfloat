CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native


ALLEXE :=
ALLEXE += mm_speed
ALLEXE += mm_speed_afmf
ALLEXE += mm_speed_armf
ALLEXE += mm_speed_afmr
ALLEXE += mm_speed_armr
# ALLEXE += mm_prec
# ALLEXE += mm_prec_ma

all: ${ALLEXE}


mulq.o:../../mulq.c
	$(CC) $(COPT) -c $< -o $@

addsubq.o:../../addsubq.c
	$(CC) $(COPT) -c $< -o $@


mulq_rm.o:../../mulq_rm.c
	$(CC) $(COPT) -c $< -o $@

addsubq_rm.o:../../addsubq_rm.c
	$(CC) $(COPT) -c $< -o $@


matmulq.o:../matmulq/matmulq.c ../matmulq/matmulq.h
	$(CC) $(COPT) -c $< -o $@

qmx2mpfr.o:../util/qmx2mpfr.c ../util/qmx2mpfr.h
	$(CC) $(COPT) -c $< -o $@

mm_speed.o:../matmulq/mm_speed.cpp ../matmulq/matmulq.h
	$(CPP) $(COPT) -c $< -o $@

mm_prec.o:../matmulq/mm_prec.cpp ../matmulq/matmulq.h ../util/qmx2mpfr.h
	$(CPP) $(COPT)  -I../util -c $<


mm_speed:mm_speed.o matmulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=$@.map

mm_speed_afmf:mm_speed.o matmulq.o addsubq.o mulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=$@.map

mm_speed_armf:mm_speed.o matmulq.o addsubq_rm.o mulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=$@.map

mm_speed_afmr:mm_speed.o matmulq.o addsubq.o mulq_rm.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=$@.map

mm_speed_armr:mm_speed.o matmulq.o addsubq_rm.o mulq_rm.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=$@.map


mm_prec:mm_prec.o qmx2mpfr.o matmulq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=$@.map

mm_prec_ma:mm_prec.o qmx2mpfr.o matmulq.o addsubq.o mulq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=$@.map
