CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native


ALLEXE :=
ALLEXE += mm_speed
ALLEXE += mm_speed_afmf
ALLEXE += mm_speed_armf
ALLEXE += mm_speed_afmr
ALLEXE += mm_speed_armr
ALLEXE += mm_speed_afmr_mocks
ALLEXE += mm_speed_armr_mocks
ALLEXE += mm_speed_armf_mocks
ALLEXE += mm_speed_armr_mocks2
# ALLEXE += mm_prec
# ALLEXE += mm_prec_ma

all: ${ALLEXE}


# compilation

obj/mulq.o:../../mulq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/addsubq.o:../../addsubq.c
	$(CC) $(COPT) -c -MD $< -o $@


obj/mulq_rm.o:../../mulq_rm.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/addsubq_rm.o:../../addsubq_rm.c
	$(CC) $(COPT) -c -MD $< -o $@


obj/mulq_rm_mocks.o:../../mulq_rm_mocks.c
	$(CC) $(COPT) -c -MD $< -o $@


obj/addsubq_rm_mocks.o:../../addsubq_rm_mocks.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/matmulq.o:../matmulq/matmulq.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mm_speed.o:../matmulq/mm_speed.cpp
	$(CPP) $(COPT) -c -MD $< -o $@

obj/mm_prec.o:../matmulq/mm_prec.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

-include obj/mulq.d
-include obj/addsubq.d
-include obj/mulq_rm.d
-include obj/addsubq_rm.d
-include obj/mulq_rm_mocks.d
-include obj/addsubq_rm_mocks.d
-include obj/matmulq.d
-include obj/qmx2mpfr.d
-include obj/mm_speed.d
-include obj/mm_prec.d

# link

mm_speed:obj/mm_speed.o obj/matmulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_afmf:obj/mm_speed.o obj/matmulq.o obj/addsubq.o obj/mulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_armf:obj/mm_speed.o obj/matmulq.o obj/addsubq_rm.o obj/mulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_afmr:obj/mm_speed.o obj/matmulq.o obj/addsubq.o obj/mulq_rm.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_armr:obj/mm_speed.o obj/matmulq.o obj/addsubq_rm.o obj/mulq_rm.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_afmr_mocks:obj/mm_speed.o obj/matmulq.o obj/addsubq.o obj/mulq_rm_mocks.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_armr_mocks:obj/mm_speed.o obj/matmulq.o obj/addsubq_rm.o obj/mulq_rm_mocks.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_armf_mocks:obj/mm_speed.o obj/matmulq.o obj/addsubq_rm_mocks.o obj/mulq.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map

mm_speed_armr_mocks2:obj/mm_speed.o obj/matmulq.o obj/addsubq_rm_mocks.o obj/mulq_rm_mocks.o
	$(CPP) -s $^ -lm -lquadmath -o $@ -Xlinker -Map=./map/$@.map


mm_prec:obj/mm_prec.o obj/qmx2mpfr.o obj/matmulq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=./map/$@.map

mm_prec_ma:obj/mm_prec.o obj/qmx2mpfr.o obj/matmulq.o obj/addsubq.o obj/mulq.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@ -Xlinker -Map=./map/$@.map


$(info $(shell mkdir -p obj map))

clean:
	@rm -f obj/*  map/*.map
