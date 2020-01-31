CC = gcc
CPP = g++
COPT = -Wall -O2

all: prec_test.exe speed_test.exe

prec_test.o: prec_test.cpp sqrt448.h
	${CPP} ${COPT} -c $<

speed_test.o: speed_test.cpp sqrt448.h
	${CPP} ${COPT} -c $<

sqrt448.o: sqrt448.cpp sqrt448.h
	${CPP} ${COPT} -c $<

prec_test.exe : prec_test.o sqrt448.o
	${CPP} $+ -o $@

speed_test.exe : speed_test.o sqrt448.o
	${CPP} $+ -o $@
