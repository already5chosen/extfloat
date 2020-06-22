CC = gcc
AS = clang
CPP = g++
COPT = -Wall -O2

ALL :=
ALL += prec_test.exe speed_test.exe
ALL += prec_test.x2.exe speed_test.x2.exe
ALL += prec_test.w.exe speed_test.w.exe
ALL += prec_test.q.exe speed_test.q.exe
ALL += prec_test.q2.exe speed_test.q2.exe
ALL += prec_test.q2s.exe speed_test.q2s.exe
ALL += prec_test.q2m.exe speed_test.q2m.exe
ALL += prec_test_a.exe speed_test_a.exe
ALL += prec_test_a.x2.exe speed_test_a.x2.exe
ALL += prec_test_a.w.exe speed_test_a.w.exe
ALL += prec_test_a.q.exe speed_test_a.q.exe
ALL += prec_test_a.q2.exe speed_test_a.q2.exe
ALL += prec_test_a.q2l.exe speed_test_a.q2l.exe
ALL += prec_test_a.q2m.exe speed_test_a.q2m.exe
all: ${ALL}

prec_test.o: prec_test.cpp sqrt448.h
	${CPP} ${COPT} -c $<

speed_test.o: speed_test.cpp sqrt448.h
	${CPP} ${COPT} -c $<

sqrt448.o: sqrt448.cpp sqrt448.h
	${CPP} ${COPT} -c $<

sqrt448_%.o: sqrt448_%.cpp sqrt448.h
	${CPP} ${COPT} -c $<

asm_%.o: asm_%.s
	${AS} ${COPT} -c $<

prec_test.exe : prec_test.o sqrt448.o
	${CPP} $+ -o $@

speed_test.exe : speed_test.o sqrt448.o
	${CPP} $+ -o $@

prec_test.x2.exe : prec_test.o sqrt448_x2.o
	${CPP} $+ -o $@

speed_test.x2.exe : speed_test.o sqrt448_x2.o
	${CPP} $+ -o $@

prec_test.w.exe : prec_test.o sqrt448_w.o
	${CPP} $+ -o $@

speed_test.w.exe : speed_test.o sqrt448_w.o
	${CPP} $+ -o $@

prec_test.q.exe : prec_test.o sqrt448_q.o
	${CPP} $+ -o $@

speed_test.q.exe : speed_test.o sqrt448_q.o
	${CPP} $+ -o $@

prec_test.q2.exe : prec_test.o sqrt448_q2.o
	${CPP} $+ -o $@

prec_test.q2s.exe : prec_test.o sqrt448_q2s.o
	${CPP} $+ -o $@

speed_test.q2s.exe : speed_test.o sqrt448_q2s.o
	${CPP} $+ -o $@

prec_test.q2m.exe : prec_test.o sqrt448_q2m.o
	${CPP} $+ -o $@

speed_test.q2m.exe : speed_test.o sqrt448_q2m.o
	${CPP} $+ -o $@

speed_test.q2.exe : speed_test.o sqrt448_q2.o
	${CPP} $+ -o $@

prec_test_a.exe : prec_test.o sqrt448_a.o asm_sqrt448.o
	${CPP} $+ -o $@

speed_test_a.exe : speed_test.o sqrt448_a.o asm_sqrt448.o
	${CPP} $+ -o $@

prec_test_a.x2.exe : prec_test.o sqrt448_a.o asm_sqrt448_x2.o
	${CPP} $+ -o $@

speed_test_a.x2.exe : speed_test.o sqrt448_a.o asm_sqrt448_x2.o
	${CPP} $+ -o $@

prec_test_a.w.exe : prec_test.o sqrt448_a.o asm_sqrt448_w.o
	${CPP} $+ -o $@

speed_test_a.w.exe : speed_test.o sqrt448_a.o asm_sqrt448_w.o
	${CPP} $+ -o $@

prec_test_a.q.exe : prec_test.o sqrt448_a.o asm_sqrt448_q.o
	${CPP} $+ -o $@

speed_test_a.q.exe : speed_test.o sqrt448_a.o asm_sqrt448_q.o
	${CPP} $+ -o $@

prec_test_a.q2.exe : prec_test.o sqrt448_a.o asm_sqrt448_q2.o
	${CPP} $+ -o $@

speed_test_a.q2.exe : speed_test.o sqrt448_a.o asm_sqrt448_q2.o
	${CPP} $+ -o $@

prec_test_a.q2l.exe : prec_test.o sqrt448_a.o asm_sqrt448_q2l.o
	${CPP} $+ -o $@

speed_test_a.q2l.exe : speed_test.o sqrt448_a.o asm_sqrt448_q2l.o
	${CPP} $+ -o $@

prec_test_a.q2m.exe : prec_test.o sqrt448_a.o asm_sqrt448_q2m.o
	${CPP} $+ -o $@

speed_test_a.q2m.exe : speed_test.o sqrt448_a.o asm_sqrt448_q2m.o
	${CPP} $+ -o $@
