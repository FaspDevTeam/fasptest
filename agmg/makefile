#=========================================================================#
#  SPECS (Serial and Parallel, Efficient, Cunning Solvers) (c) 2010-2011  #
#             School of Mathematics and Computational Science             #
#        Xiaoqiang Yue, Zhiyang Zhou, Zheng Lee and Chunsheng Feng        #
#                           Xiangtan University                           #
#=========================================================================#

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Suffixes
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.SUFFIXES:.c .f90 .f .for

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compiler
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
F90C = gfortran
MyCC = gcc

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Lib
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
MATHLIB = -lm
BLASLIB = -lblas
AGMG_LIB = -L /usr/lib
LAPACKLIB = -llapack
EXTENSELIB = ${MATHLIB} ${BLASLIB} ${LAPACKLIB}
SPECSLIBFLAGS = ${AGMG_LIB} ${EXTENSELIB}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compile Option
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C_COMPILE_FLAGS = -O2 -Wall
F_COMPILE_FLAGS = -O2 -Wall
CFLAGS = ${C_COMPILE_FLAGS}
FFLAGS = ${F_COMPILE_FLAGS}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  Rules for Compilation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
.c.o:
	${MyCC} -o $@ -c ${CFLAGS} $<
.f.o:
	${F90C} -o $@ -c ${FFLAGS} $<
.for.o:
	${F90C} -o $@ -c ${FFLAGS} $<
.f90.o:
	${F90C} -o $@ -c ${FFLAGS} $<

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Object
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
agmgdir = /home/stone/Tests/AGMG_3.2.0-aca/SRC

OBJS = driver.o io.o $(agmgdir)/dagmg.o $(agmgdir)/dagmg_mumps.o

test: ${OBJS}
	${F90C} -o stone ${OBJS} ${SPECSLIBFLAGS}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Clean
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
clean:
	-rm -f *.o
	-rm -f *~
	-rm -f stone
