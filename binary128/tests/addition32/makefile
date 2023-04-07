CC=gcc -m32
CPP=g++ -m32
COPT=-Wall -O2 -march=native

ALL_RESULTS := tadd
ALL_RESULTS += tsub
ALL_RESULTS += tadd_man
ALL_RESULTS += tadd_man_gcclib

all: $(ALL_RESULTS)

# compilation

obj/addsubq32.o:../../addsubq32.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mpfr_strtofr_clipped_to_float128.o:../util/mpfr_strtofr_clipped_to_float128.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/tadd_man.o:../addition/tadd_man.c
	$(CC) $(COPT) -c -MD $< -o $@ -I../util

obj/tadd.o:../addition/tadd.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

obj/tsub.o:../addition/tadd.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util -DSUBTRACT


-include obj/addsubq32.d
-include obj/qmx2mpfr.d
-include obj/mpfr_strtofr_clipped_to_float128.d
-include obj/tadd_man.d
-include obj/tadd.d
-include obj/tsub.d

# link

tadd:obj/tadd.o obj/addsubq32.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsub:obj/tsub.o obj/addsubq32.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tadd_man:obj/tadd_man.o obj/addsubq32.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@

tadd_man_gcclib:obj/tadd_man.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@


$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
