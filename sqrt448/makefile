CC = gcc
AS = clang
CPP = g++
COPT = -Wall -O2

all: prec_test.exe speed_test.exe prec_test_a.exe speed_test_a.exe

prec_test.o: prec_test.cpp sqrt448.h
	${CPP} ${COPT} -c $<

speed_test.o: speed_test.cpp sqrt448.h
	${CPP} ${COPT} -c $<

sqrt448.o: sqrt448.cpp sqrt448.h
	${CPP} ${COPT} -c $<

sqrt448_a.o: sqrt448_a.cpp sqrt448.h
	${CPP} ${COPT} -c $<

asm_sqrt448.o: asm_sqrt448.s
	${AS} ${COPT} -c $<

prec_test.exe : prec_test.o sqrt448.o
	${CPP} $+ -o $@

speed_test.exe : speed_test.o sqrt448.o
	${CPP} $+ -o $@

prec_test_a.exe : prec_test.o sqrt448_a.o asm_sqrt448.o
	${CPP} $+ -o $@

speed_test_a.exe : speed_test.o sqrt448_a.o asm_sqrt448.o
	${CPP} $+ -o $@
