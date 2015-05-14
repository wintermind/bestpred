###############################################################################
# NAME:         makefile for GNU Fortran compiler
# VERSION:      2.0 rc 7
# RELEASED:     24 NOVEMBER 2014
# AUTHORS:      Paul M. VanRaden (paul@aipl.arsusda.gov)
#               John B. Cole (john.cole@ars.usda.gov)
# DESCRIPTION:  makefile for the bestpred package from AIPL.
###############################################################################
all: ageadjs.o aiplage.o bestpred.o bestpred_fmt4.o bestpred_parm.o bestpred_log.o bestpred

ageadjs.o: ageadjs.c
	gcc -L. -I. -c ageadjs.c

aiplage.o: aiplage.c aiplage.h
	gcc -L. -I. -c aiplage.c

#ageadj.o: ageadj.c ageadj.h
#	gcc -c ageadj.c

bestpred_parm.o: bestpred_parm.f90
	gfortran -v -O3 -ggdb -L. -I. -c bestpred_parm.f90 \
		-fno-underscoring \
		-fno-automatic \
		-fimplicit-none

bestpred_log.o: bestpred_log.f90
	gfortran -v -O3 -ggdb -L. -I. -c bestpred_log.f90 \
		-fno-underscoring \
		-fno-automatic \
		-fimplicit-none

bestpred.o: bestpred.f90 bestpred_parm.o bestpred_log.o
	gfortran -v -O3 -ggdb -L. -I. -c bestpred.f90 \
		-fno-underscoring \
		-fno-automatic \
                -fimplicit-none \
		bestpred_parm.o \
		bestpred_log.o

bestpred_fmt4.o: bestpred_fmt4.f90 bestpred.o ageadjs.o aiplage.o
	gfortran -v -O3 -ggdb -L. -I. -c bestpred_fmt4.f90 \
		-fno-underscoring \
		-fno-automatic \
		-fimplicit-none \
		ageadjs.o \
		aiplage.o

bestpred: bestpred_main.f90 bestpred_parm.o bestpred_log.o bestpred_fmt4.o bestpred.o ageadjs.o aiplage.o
	gfortran -v -O3 -ggdb -L. -I. bestpred_main.f90 -o bestpred \
		-fno-underscoring \
		-fno-automatic \
		-fimplicit-none \
		bestpred_fmt4.o \
		bestpred.o \
		bestpred_parm.o \
		bestpred_log.o \
		aiplage.o \
		ageadjs.o
