CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native


ALLEXE :=
ALLEXE += mm_speed
ALLEXE += mm_speed_ma
ALLEXE += mm_prec
ALLEXE += mm_prec_ma

all: ${ALLEXE}

# compilation

obj/mulq.o:../../mulq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/addsubq.o:../../addsubq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/matmulq.o:matmulq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mm_speed.o:mm_speed.cpp
	$(CPP) $(COPT) -c -MD $< -o $@

obj/mm_prec.o:mm_prec.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

-include obj/mulq.d
-include obj/addsubq.d
-include obj/matmulq.d
-include obj/qmx2mpfr.d
-include obj/mm_speed.d
-include obj/mm_prec.d

# link

mm_speed:obj/mm_speed.o obj/matmulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_ma:obj/mm_speed.o obj/matmulq.o obj/addsubq.o obj/mulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_prec:obj/mm_prec.o obj/qmx2mpfr.o obj/matmulq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=./map/$@.map

mm_prec_ma:obj/mm_prec.o obj/qmx2mpfr.o obj/matmulq.o obj/addsubq.o obj/mulq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=./map/$@.map


$(info $(shell mkdir -p obj map))

clean:
	@rm -f obj/*  map/*.map
