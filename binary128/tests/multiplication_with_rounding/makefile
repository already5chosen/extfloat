CC=gcc
CPP=g++
COPT=-Wall -O2 -march=native

ALL_RESULTS := tmul
ALL_RESULTS += tmul_man
ALL_RESULTS += tmul_libgcc
ALL_RESULTS += tmul_man_libgcc
ALL_RESULTS += tmul_mocks
ALL_RESULTS += tmul_man_mocks

all: $(ALL_RESULTS)

# compilation

obj/mulq.o:../../mulq_rm.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mulq_rm_mocks.o:../../mulq_rm_mocks.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/qmx2mpfr.o:../util/qmx2mpfr.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mpfr_strtofr_clipped_to_float128.o:../util/mpfr_strtofr_clipped_to_float128.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mulq_test_generator.o:../util/mulq_test_generator.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/mulq_with_all_rounding_modes.o:mulq_with_all_rounding_modes.c
	$(CC) $(COPT) -c -MD $< -o $@

obj/tmul_man.o:tmul_man.c
	$(CC) $(COPT) -c -MD $< -o $@ -I../util

obj/tmul.o:tmul.cpp
	$(CPP) $(COPT) -c -MD $< -o $@ -I../util

-include obj/mulq.d
-include obj/mulq_rm_mocks.d
-include obj/qmx2mpfr.d
-include obj/mpfr_strtofr_clipped_to_float128.d
-include obj/mulq_test_generator.d
-include obj/mulq_with_all_rounding_modes.d
-include obj/tmul_man.d
-include obj/tmul.d

# link

tmul:obj/tmul.o obj/mulq_with_all_rounding_modes.o obj/mulq.o obj/qmx2mpfr.o obj/mulq_test_generator.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man:obj/tmul_man.o obj/mulq_with_all_rounding_modes.o obj/mulq.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@


tmul_mocks:obj/tmul.o obj/mulq_with_all_rounding_modes.o obj/mulq_rm_mocks.o obj/qmx2mpfr.o obj/mulq_test_generator.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man_mocks:obj/tmul_man.o obj/mulq_with_all_rounding_modes.o obj/mulq_rm_mocks.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@

tmul_libgcc:obj/tmul.o obj/mulq_with_all_rounding_modes.o obj/qmx2mpfr.o obj/mulq_test_generator.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -lpthread -o $@

tmul_man_libgcc:obj/tmul_man.o obj/mulq_with_all_rounding_modes.o obj/qmx2mpfr.o obj/mpfr_strtofr_clipped_to_float128.o
	$(CPP) -s $^ -lm -lquadmath -lmpfr -lgmp -o $@


$(info $(shell mkdir -p obj))

clean:
	@rm -f obj/*
