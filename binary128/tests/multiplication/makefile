CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALL_RESULTS := tmul
ALL_RESULTS += tmul_man

all: $(ALL_RESULTS)

# compilation

obj/mulq.o:../../mulq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/tmul_man.o:tmul_man.c
	$(CC) $(COPT) -c -MD $< -o $@ -I../util

obj/tmul.o:tmul.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

-include obj/mulq.d
-include obj/qmx2mpfr.d
-include obj/tmul_man.d
-include obj/tmul.d

# link

tmul:obj/tmul.o obj/mulq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man:obj/tmul_man.o obj/mulq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@


$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
