CC   = gcc
CP   = g++
#COPT = -Ofast -Wall -mavx
#CPPOPT = -std=c++11 -I"/c/Users/user/modular-boost"
COPT = -O3 -Wall -mavx
# COPT = -O3 -Wall
CPPOPT = -std=c++11 -I"/c/data/download/boost_1_62_x"

.c.o:
	${CC} -c ${COPT} $<
	${CC} -MM ${COPT} -MF $(patsubst %.o,%.d,$@) $<

.cpp.o:
	${CP} -c  ${COPT} ${CPPOPT} $<
	${CP} -MM ${COPT} ${CPPOPT} -MF $(patsubst %.o,%.d,$@) $<

# all: div_test.exe mul_test.exe add_test.exe uuu.exe cnv_to_test.exe cnvf_to_test.exe cnv_from_test.exe sqrt_test.exe
ALL_EXE :=
ALL_EXE += tst_cnv_to_d.exe
ALL_EXE += tst_cnv_to_f.exe
ALL_EXE += tst_cnv_from.exe
ALL_EXE += tst_mul.exe
ALL_EXE += tst_fma.exe
ALL_EXE += tst_div.exe
ALL_EXE += tst_add.exe
ALL_EXE += tst_sqrt.exe
ALL_EXE += tst_rsqrt.exe
ALL_EXE += tst_mod_pow2.exe
ALL_EXE += tst_trunc.exe
ALL_EXE += tst_sin.exe
ALL_EXE += tst_cos.exe
ALL_EXE += tst_sinPi.exe
ALL_EXE += tst_cosPi.exe
ALL_EXE += tst_atan.exe
ALL_EXE += tst_atan2.exe
ALL_EXE += tst_asin.exe
ALL_EXE += tst_acos.exe
ALL_EXE += mk_sin_tabs.exe
ALL_EXE += tst_cmp.exe
ALL_EXE += tst_accadd.exe
ALL_EXE += tst_accmadd.exe
all: ${ALL_EXE}

OBJLIST :=
OBJLIST += extfloat128.o
OBJLIST += extfloat128_acc.o
OBJLIST += extfloat128_sincos.o
OBJLIST += extfloat128_atan.o
OBJLIST += tst_div.o
OBJLIST += tst_mul.o
OBJLIST += tst_add.o
OBJLIST += tst_fma.o
OBJLIST += tst_sqrt.o
OBJLIST += tst_rsqrt.o
OBJLIST += tst_cnv_to_d.o
OBJLIST += tst_cnv_to_f.o
OBJLIST += tst_cnv_from.o
OBJLIST += tst_mod_pow2.o
OBJLIST += tst_trunc.o
OBJLIST += tst_sin.o
OBJLIST += tst_cos.o
OBJLIST += tst_sinPi.o
OBJLIST += tst_cosPi.o
OBJLIST += tst_atan.o
OBJLIST += tst_atan2.o
OBJLIST += tst_asin.o
OBJLIST += tst_acos.o
OBJLIST += mk_sin_tabs.o
OBJLIST += tst_cmp.o
OBJLIST += tst_accadd.o
OBJLIST += tst_accmadd.o

-include $(OBJLIST:.o=.d)

tst_cnv_to_d.exe : tst_cnv_to_d.o extfloat128.o
	${CP} $+ -o$@

tst_cnv_to_f.exe : tst_cnv_to_f.o extfloat128.o
	${CP} $+ -o$@

tst_cnv_from.exe : tst_cnv_from.o extfloat128.o
	${CP} $+ -o$@

tst_mul.exe : tst_mul.o extfloat128.o
	${CP} $+ -o$@

tst_div.exe : tst_div.o extfloat128.o
	${CP} $+ -o$@

tst_add.exe : tst_add.o extfloat128.o
	${CP} $+ -o$@

tst_sqrt.exe : tst_sqrt.o extfloat128.o
	${CP} $+ -o$@

tst_rsqrt.exe : tst_rsqrt.o extfloat128.o
	${CP} $+ -o$@

tst_fma.exe : tst_fma.o extfloat128.o
	${CP} $+ -o$@

tst_mod_pow2.exe : tst_mod_pow2.o extfloat128.o
	${CP} $+ -o$@

tst_trunc.exe : tst_trunc.o extfloat128.o
	${CP} $+ -o$@

tst_sin.exe : tst_sin.o extfloat128.o extfloat128_acc.o extfloat128_sincos.o
	${CP} $+ -o$@

tst_cos.exe : tst_cos.o extfloat128.o extfloat128_acc.o extfloat128_sincos.o
	${CP} $+ -o$@

tst_sinPi.exe : tst_sinPi.o extfloat128.o extfloat128_acc.o extfloat128_sincos.o
	${CP} $+ -o$@

tst_cosPi.exe : tst_cosPi.o extfloat128.o extfloat128_acc.o extfloat128_sincos.o
	${CP} $+ -o$@

tst_atan.exe : tst_atan.o extfloat128.o extfloat128_atan.o
	${CP} $+ -o$@

tst_atan2.exe : tst_atan2.o extfloat128.o extfloat128_atan.o
	${CP} $+ -o$@

tst_asin.exe : tst_asin.o extfloat128.o extfloat128_atan.o
	${CP} $+ -o$@

tst_acos.exe : tst_acos.o extfloat128.o extfloat128_acc.o extfloat128_atan.o extfloat128_sincos.o
	${CP} $+ -o$@

mk_sin_tabs.exe : mk_sin_tabs.o extfloat128.o
	${CP} $+ -o$@

tst_cmp.exe : tst_cmp.o extfloat128.o
	${CP} $+ -o$@

tst_accadd.exe : tst_accadd.o extfloat128.o extfloat128_acc.o
	${CP} $+ -o$@

tst_accmadd.exe : tst_accmadd.o extfloat128.o extfloat128_acc.o
	${CP} $+ -o$@

clean:
	rm *.o *.d $(ALL_EXE)