CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALL_RESULTS := taddv
ALL_RESULTS += tsubv

all: $(ALL_RESULTS)

# compilation

obj/addsubvq.o:../../vector/addsubvq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/taddv.o:taddv.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util -I../../vector

obj/tsubv.o:taddv.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util -I../../vector -DSUBTRACT

-include obj/addsubvq.d
-include obj/qmx2mpfr.d
-include obj/taddv.d
-include obj/tsubv.d

# link

taddv:obj/taddv.o obj/addsubvq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsubv:obj/tsubv.o obj/addsubvq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@


$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
