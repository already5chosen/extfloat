CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALL_RESULTS := tmul
ALL_RESULTS += tmul_man
ALL_RESULTS += tmul_gcclib
ALL_RESULTS += tmul_man_gcclib

all: $(ALL_RESULTS)

# compilation

obj/mulq.o:../../mulq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mpfr_strtofr_clipped_to_float128.o:../util/mpfr_strtofr_clipped_to_float128.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/tmul_man.o:tmul_man.c
	$(CC) $(COPT) -c -MD $< -o $@ -I../util

obj/tmul.o:tmul.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

-include obj/mulq.d
-include obj/qmx2mpfr.d
-include obj/mpfr_strtofr_clipped_to_float128.d
-include obj/tmul_man.d
-include obj/tmul.d

# link

tmul:obj/tmul.o obj/mulq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man:obj/tmul_man.o obj/mulq.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@

tmul_gcclib:obj/tmul.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man_gcclib:obj/tmul_man.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@


$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
