CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALL_RESULTS := tadd
ALL_RESULTS += tsub
ALL_RESULTS += tadd_man

all: $(ALL_RESULTS)

# compilation

obj/addsubq.o:../../addsubq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/tadd_man.o:tadd_man.c
	$(CC) $(COPT) -c -MD $< -o $@ -I../util

obj/tadd.o:tadd.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

obj/tsub.o:tadd.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util -DSUBTRACT

-include obj/addsubq.d
-include obj/qmx2mpfr.d
-include obj/tadd_man.d
-include obj/tadd.d
-include obj/tsub.d

# link

tadd:obj/tadd.o obj/addsubq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsub:obj/tsub.o obj/addsubq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tadd_man:obj/tadd_man.o obj/addsubq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@


$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
