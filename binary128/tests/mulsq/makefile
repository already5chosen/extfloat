CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALLEXE :=
ALLEXE += tmuls

all: ${ALLEXE}

# compilation

obj/mulsq.o:../../vector/mulsq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/tmuls.o:tmuls.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util -I../../vector

-include obj/mulsq.d
-include obj/qmx2mpfr.d
-include obj/tmuls.d

# link

tmuls:obj/tmuls.o obj/mulsq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@


$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
