CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALL_RESULTS := tadd
ALL_RESULTS += tsub
ALL_RESULTS += tadd_libgcc
ALL_RESULTS += tsub_libgcc
ALL_RESULTS := tadd_mocks
ALL_RESULTS += tsub_mocks
ALL_RESULTS += tadd_man
ALL_RESULTS += tadd_man_libgcc
ALL_RESULTS += tadd_man_mocks

all: $(ALL_RESULTS)

# compilation

obj/addsubq.o:../../addsubq_rm.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/addsubq_rm_mocks.o:../../addsubq_rm_mocks.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mpfr_strtofr_clipped_to_float128.o:../util/mpfr_strtofr_clipped_to_float128.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/addsubq_with_all_rounding_modes.o:addsubq_with_all_rounding_modes.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/tadd_man.o:tadd_man.c
	$(CC) $(COPT) -c -MD $< -o $@ -I../util

obj/tadd.o:tadd.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

obj/tsub.o:tadd.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util -DSUBTRACT

-include obj/addsubq.d
-include obj/addsubq_rm_mocks.d
-include obj/qmx2mpfr.d
-include obj/mpfr_strtofr_clipped_to_float128.d
-include obj/tadd_man.d
-include obj/tadd.d
-include obj/tsub.d

# link

tadd:obj/tadd.o obj/addsubq_with_all_rounding_modes.o obj/addsubq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsub:obj/tsub.o obj/addsubq_with_all_rounding_modes.o obj/addsubq.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tadd_man:obj/tadd_man.o obj/addsubq_with_all_rounding_modes.o obj/addsubq.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@



tadd_libgcc:obj/tadd.o obj/addsubq_with_all_rounding_modes.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsub_libgcc:obj/tsub.o obj/addsubq_with_all_rounding_modes.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tadd_man_libgcc:obj/tadd_man.o obj/addsubq_with_all_rounding_modes.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@



tadd_mocks:obj/tadd.o obj/addsubq_with_all_rounding_modes.o obj/addsubq_rm_mocks.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tsub_mocks:obj/tsub.o obj/addsubq_with_all_rounding_modes.o obj/addsubq_rm_mocks.o obj/qmx2mpfr.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tadd_man_mocks:obj/tadd_man.o obj/addsubq_with_all_rounding_modes.o obj/addsubq_rm_mocks.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@



$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
