#
##
## Makefile for the mpi-parallel SeqQuest code
##
## Last revised:  5Dec12-PAS - chama libraries
## Configured for parallel: openmpi, mkl with integrated scalapack/blacs
##
#

######
###### Choose compiler (pick one and only one):
######
## general:
# FC = f77
## g77:
# FC = g77
## pgi compiler:
# FC = /opt/pgi/linux86-64/7.2-1/bin/pgf77
## f77+mpi:

# Here is the parallel compiler
#
FC = mpif77

######
###### Set compiler flags (pick one and only one):
######
## general/preferred:
# FFLAGS = -O -fast
## general/when -fast is not available, or with pgi compiler:
# FFLAGS = -O
## for some g77 type-change (type-casting?) causes compile failure:
# FFLAGS = -O -fno-globals
## on SGI machines:
# FFLAGS = -O3 -mips4 -OPT:IEEE_arithmetic=3:alias=restrict
FFLAGS = -O0 -g -fbacktrace -std=legacy

## everyone else but chama
LFLAGS =
## chama
#LFLAGS = -mcmodel medium -shared-intel
LIBS = 

######
###### Math libraries (to pick up mpi, basics for math and links, lapack/blas ... )
######
#
# Select appropriate math library to link in.
# Internally, code default is to make use of blas/3 routines,
# an option that is flag-controlled (doblas3).  With an
# optimized blas library, this gives the best performance,
# but without an optimized library, using blas/3 routines
# gives much slower execution.  If your platform does not
# have an optimized library and you must use the explicit
# source provided below, edit subroutine "dat0opt" in the
# file lcaosubs.f to set default "doblas3" to .false. for
# most efficient execution.
#
## dec math libraries:
# MATH = -ldxml
## compaq renamed math libraries:
# MATH = -lcxml
## SGI math libraries:
# MATH = -lcomplib.sgimath
## linux:
# MATH = -llapack -lblas
## linux/ibm:
# MATH = -static -llapack -lblas
## explicit source as provided:
# MATH = lapack.o blas.o
## pgi compiler: likely overcomplete, but it works
# MATH = -L/opt/pgi/linux86-64/7.2-1/lib/ -llapack -lblas -lpgftnrtl -lnspgc -lpgc -lm -lstd -lacml

###### Math libraries for parallel

## MPI link (simple):
 MATH = -llapack -lblas -lmpi

## ICC/liberty math libraries:
## MATH = -L/apps/intel/mkl70/lib/32 -lmkl_ia32 -lguide -lmkl_lapack
# MATH = -L/apps/intel/mkl70/lib/32 -lmkl -lmkl_ia32 -lguide -lmkl_lapack \
#       -L/apps/intel/fce/lib -lsvml
## NWCC/spirit-old math libraries:
# MATH = -L/apps/intel/mkl701/lib/em64t -lmkl -lguide -lmkl_lapack64
## NWCC/spirit-openmpi math libraries:
# MATH = -L/apps/x86_64/mpi/openmpi/intel-9.1-f040-c045/openmpi-1.2.2_mx/lib -lmpi \
#        -L/projects/global/x86_64/libraries/l_mkl_p_9.1.018/lib/em64t \
#        -lmkl -lguide -lmkl_lapack

##
## >>>>> openmpi-1.2.4/mkl8 (had problems with mkl9), with separate scalapack
##
#MATH = -L/apps/x86_64/mpi/openmpi/intel-10.0-f025-c025/openmpi-1.2.4_ofed/lib -lmpi \
#       -L/projects/global/x86_64/libraries/l_mkl_p_8.0.1.006/lib/em64t \
#       -lmkl -lguide -lmkl_lapack64
##
## >>>>> OPENMPI with MKL10 with integrated scalapack 
##
#MKL_LIB = /projects/global/x86_64/libraries/l_mkl_10.0.3.020/lib/em64t
#MATH = -L/apps/x86_64/mpi/openmpi/intel-10.0-f025-c025/openmpi-1.2.4_ofed/lib -lmpi \
#       -L/projects/global/x86_64/libraries/l_mkl_10.0.3.020/lib/em64t \
#       -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide
##
## >>>>> MKL with integrated scalapack, redsky
##
#MKL_LIB = /projects/global/x86_64/compilers/intel/intel-11.1-cprof-064/mkl/lib/em64t
#MATH = -L/apps/x86_64/mpi/openmpi/intel-11.1-f064-c064/openmpi-1.4.2_static_nomemmgr_oobpr/lib/ -lmpi \
#       -L/projects/global/x86_64/compilers/intel/intel-11.1-cprof-064/mkl/lib/em64t \
#       -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lguide
##       -lmkl_lapack -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
##


######
###### scalapack parallel math libraries
######
##
####  First, define patch to scalapack and blacs libraries
##
## downloaded and built scalapack parallel math libraries in the home directory

 PMATHHOME = /opt/packages/scalapack/2.1.0

# PBLACSHOME = /apps/libraries/blacs
## NWCC/spirit-openmpi
# PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.2.2_mx_intel-9.1-f040-c045/scalapack-1.7.4/lib
# PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.2.2_mx_intel-9.1-f040-c045/BLACS-1.1/lib
## tbird: mva-0.9.8
#PMATHHOME = /apps/x86_64/libraries/scalapack/mvapich-0.9.8-ofed-intel-9.1/scalapack-1.7.4/lib
#PBLACSHOME = /apps/x86_64/libraries/BLACS/mvapich-0.9.8-ofed-intel-9.1/BLACS-1.1/lib
## tbird: openmpi-1.1.2
#PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.1.2-ofed-intel-9.1/scalapack-1.7.4/lib
#PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.1.2-ofed-intel-9.1/BLACS-1.1/lib
## tbird: openmpi-1.2.4/mkl8 (had problems with mkl9)
#PMATHHOME = /apps/x86_64/libraries/scalapack/openmpi-1.2.4_ofed_intel-10.0-f025-c025/scalapack-1.7.4/lib
#PBLACSHOME = /apps/x86_64/libraries/BLACS/openmpi-1.2.4_ofed_intel-10.0-f025-c025/BLACS-1.1/lib
##
#### >>>>> build pmath to libraries
##
PMATH =  \
-L ${PMATHHOME} -lscalapack
#${PBLACSHOME}/blacsCinit_MPI-LINUX-0.a \
#${PBLACSHOME}/blacsF77init_MPI-LINUX-0.a \
#${PBLACSHOME}/blacs_MPI-LINUX-0.a
##
## >>>>> MKL with integrated scalapack libraries
##
#chama#PMATH = -L${MKL_LIB} -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64

#########################################################################
#########################################################################
#########################################################################

######
###### Source files
######

qvrs = 265a

## main program
#MVRS=lcaomain.o
MVRSTP=lcaomain_tp.o

## program subroutines
SUBS=lcaosubs.o
TPSUBS=lcaosubs_tp.o
#SUBS=subs/subs.a
#TPSUBS=subs_tp/subs_tp.a

######
###### Platform-dependent configuration file (serial or mpi-parallel)
######

## general serial code
# SUBSU = utl_gen.o
## general MPI code
SUBSU = utl_mpi.o


######
###### Where it all actually happens
######

quest: quest_${qvrs}.x

help:
	@echo "Make SeqQuest program - openmpi parallel"
	@echo "Usage: make <platform>, only platform available: quest"
	@echo "For parallel builds, you will have to specify library locations in makefile"
	@echo ""

clean:
	rm -f *.o quest_${qvrs}.x

######
###### The actual build
######

quest_${qvrs}.x : ${MVRSTP} ${SUBS} ${TPSUBS} ${SUBSU}
	$(FC) $(LFLAGS) -o $@  ${MVRSTP} ${SUBSU} ${TPSUBS} ${SUBS} ${PMATH} ${MATH}


