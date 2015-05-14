# Compile Data Collection Rating programs
# Compile age adjustment programs            (makefile)
#
all: ageadjs.o ageadj.o aiplage.o bestpred.o bestpred_fmt4.o bestpred_parm.o  bestpred_log.o bestpred

ageadjs.o: ageadjs.c
	gcc -ggdb -O2 -L. -I. -c ageadjs.c

aiplage.o: aiplage.c aiplage.h
	gcc -ggdb -O2 -L. -I. -c aiplage.c

ageadj.o: ageadj.c ageadj.h
	gcc -ggdb -O2 -L. -I. -c ageadj.c

bestpred_parm.o: bestpred_parm.f90
	gfortran -ggdb -O2 -L. -I. -c bestpred_parm.f90 \
		-fno-automatic -fno-underscoring

bestpred_log.o: bestpred_log.f90
	gfortran -ggdb -O2 -L. -I. -c bestpred_log.f90 \
		-fno-automatic -fno-underscoring

bestpred.o: bestpred.f90 bestpred_log.o
	gfortran -ggdb -O2 -L. -I. -c bestpred.f90 \
		-fno-automatic -fno-underscoring \
		bestpred_log.o

bestpred_fmt4.o: bestpred_fmt4.f90 ageadjs.o aiplage.o bestpred.o
	gfortran -ggdb -O2 -L. -I. -c bestpred_fmt4.f90 \
		-fno-automatic -fno-underscoring \
		ageadjs.o \
		aiplage.o \
		bestpred.o

bestpred: bestpred_main.f90 bestpred_parm.o bestpred_log.o bestpred_fmt4.o bestpred.o ageadjs.o aiplage.o
	  gfortran -v -ggdb -O2 -L. -I. bestpred_main.f90 -o bestpred \
		-fno-underscoring -fno-automatic \
		bestpred_fmt4.o	\
		ageadjs.o \
                aiplage.o \
                bestpred.o \
		bestpred_log.o \
		bestpred_parm.o
