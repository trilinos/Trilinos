## Set the C, Fortran compilers,flags etc. for PGI compilers (Redhat 6.X)
CC_pgilinux = pgcc 
FC_pgilinux = pgf77
LD_pgilinux = pgcc
# flag for the gcc compiler. 
CFLAGS_pgilinux = -O3 -DHAVE_FORTRAN_UNDERSCORE
#flags for the Fortran compiler.
FFLAGS_pgilinux = -O3
SYSLIBS_pgilinux = -lm -lpgftnrtl -lg2c
AR_pgilinux = ar
RANLIB_pgilinux = ranlib

# Threaded BLAS
#BLAS_LIB_pgilinux = -lblas-2 -lpthread

# Optimized single PE BLAS
#BLAS_LIB_pgilinux = -lblas

# Fortran BLAS
BLAS_LIB_pgilinux = -lblasf

LAPACK_LIB_pgilinux = -llapack
TIMER_OBJ_pgilinux = second_linux.o

## Set the C, Fortran compilers,flags etc. for Gnu PC compilers (Redhat 6.X)
CC_pclinux = cc 
FC_pclinux = g77
LD_pclinux = g77
# flag for the gcc compiler. 
CFLAGS_pclinux = -O0 -g -DHAVE_FORTRAN_UNDERSCORE
#flags for the Fortran compiler.
FFLAGS_pclinux = -O0 -g
SYSLIBS_pclinux = -lm
AR_pclinux = ar
RANLIB_pclinux = ranlib

# Threaded BLAS
#BLAS_LIB_pclinux = -lblas-2 -lpthread

# Optimized single PE BLAS
#BLAS_LIB_pclinux = -lblas

# Fortran BLAS
BLAS_LIB_pclinux = -lblasf

LAPACK_LIB_pclinux = -llapack
TIMER_OBJ_pclinux = second_linux.o

## Set the C, Fortran compilers,flags etc. for a SUN Ultra (Solaris 2.6).
CC_ultra = cc 
FC_ultra = f77
LD_ultra = cc
# flag for the gcc compiler. 
#CFLAGS_ultra = -ansi -Wall -g -msoft-float -DHAVE_FORTRAN_UNDERSCORE
# flags for the C compiler.
#CFLAGS_ultra = -O -Xc -XO3 -xtarget=ultra2 -xcache=16/32/1:4096/64/1
CFLAGS_ultra = -O3 -Xc -v -DHAVE_FORTRAN_UNDERSCORE
#flags for the Fortran compiler.
#FFLAGS_ultra = -O -xtarget=ultra2 -xcache=16/32/1:4096/64/1
FFLAGS_ultra = -O3
# SET the DEC Fortran libs that need to be called with a C linker.
# This works on my digital unix system.
SYSLIBS_ultra = -lF77 -lM77 -lsunmath 
AR_ultra = ar
RANLIB_ultra = ranlib
# Optimized single PE BLAS
BLAS_LIB_ultra =

LAPACK_LIB_ultra = -xlic_lib=sunperf 
TIMER_OBJ_ultra = second_linux.o


## Set the C, Fortran compilers,flags etc. for a DEC alpha.
CC_alpha = cc 
FC_alpha = f77
LD_alpha = cc
# flag for the gcc compiler. 
#CFLAGS_alpha = -ansi -Wall -g -msoft-float -DHAVE_FORTRAN_UNDERSCORE
# flags for the DEC C compiler.
CFLAGS_alpha = -O3 -warnprotos -ifo -std1 -DHAVE_FORTRAN_UNDERSCORE
#flags for the DEC Fortran compiler.
FFLAGS_alpha = -O3 -fpe2 -warn unused -warn declarations -warn argument_checking
# SET the DEC Fortran libs that need to be called with a C linker.
# This works on my digital unix system.
SYSLIBS_alpha = -lfor -lutil -lFutil -lots
AR_alpha = ar
RANLIB_alpha = ranlib
BLAS_LIB_alpha= -lblas
LAPACK_LIB_alpha = -llapack
TIMER_OBJ_alpha = second_alpha.o

## Set the C, Fortran compilers,flags etc. for SGI (IRIX 64 and above)
CC_irix64 = cc 
FC_irix64 = f77
LD_irix64 = cc
# flags for the SGI C Compiler.
CFLAGS_irix64 = -O3 -g -64 -r10000 -DHAVE_FORTRAN_UNDERSCORE
#flags for the SGI Fortran compiler.
FFLAGS_irix64 = -O3 -g -64 -r10000
# SET the SGI Fortran libs that need to be called with a C linker.
# This works on irix64 and above SGI systems.
SYSLIBS_irix64 = -O3 -g -64 -r10000 -lftn
AR_irix64 = ar
RANLIB_irix64 = true
BLAS_LIB_irix64= -lblas
LAPACK_LIB_irix64 =
TIMER_OBJ_irix64 = second_rs6000.o

# Set up the generic flags
FC = $(FC_$(ARCH))
CC = $(CC_$(ARCH))
SYSLIBS = $(SYSLIBS_$(ARCH))
LOADER = $(LD_$(ARCH))
FFLAGS = $(FFLAGS_$(ARCH))
CFLAGS = $(CFLAGS_$(ARCH)) -D$(ARCH)
AR = $(AR_$(ARCH))
RANLIB = $(RANLIB_$(ARCH))
LIB = $(LAPACK_LIB_$(ARCH)) $(BLAS_LIB_$(ARCH))
TIMER_OBJ = $(TIMER_OBJ_$(ARCH))
