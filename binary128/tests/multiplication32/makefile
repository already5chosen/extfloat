CC=gcc -m32
CPP=g++ -m32
COPT=-Wall -O2 -march=native

ALL_RESULTS := tmul
ALL_RESULTS += tmul_man

all: $(ALL_RESULTS)

# compilation

obj/mulq32.o:../../mulq32.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mpfr_strtofr_clipped_to_float128.o:../util/mpfr_strtofr_clipped_to_float128.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/tmul_man.o:../multiplication/tmul_man.c
	$(CC) $(COPT) -c -MD $< -o $@ -I../util

obj/tmul.o:../multiplication/tmul.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

-include obj/mulq32.d
-include obj/qmx2mpfr.d
-include obj/mpfr_strtofr_clipped_to_float128.d
-include obj/tmul_man.d
-include obj/tmul.d

# link

tmul:obj/tmul.o obj/mulq32.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man:obj/tmul_man.o obj/mulq32.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@


$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
