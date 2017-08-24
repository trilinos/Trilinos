#### common

INCLUDES= -I../../src

#### hasen : cuda k80

ifneq (,$(findstring hansen,$(HOST)))
#CXX=$(NVCC_WRAPPER)
MKL_CFLAG= 
KOKKOS=/home/kyukim/Work/lib/kokkos/install/cuda/gnu/release
CXXFLAGS=-O3 -g 
LIBS=
endif

#### white : cuda p100

ifneq (,$(findstring white,$(HOST)))
CXX=$(NVCC_WRAPPER)
MKL_CFLAG= 
KOKKOS=/ascldap/users/kyukim/Work/lib/kokkos/install/cuda/gnu/release
CXXFLAGS=-O3 -g -D__KOKKOSKERNELS_NVIDIA_CUBLAS__ -DFAST_MATRIX_FILL
LIBS=-lcublas
endif

#### bowman : knl 1x68x4

ifneq (,$(findstring bowman,$(HOST)))
#CXX=icpc
#MKL_CFLAG = -mkl #:sequential
LIBXSMM_ROOT = /home/kyukim/Work/lib/libxsmm/1.7.1/intel
MKL_ROOT = /home/kyukim/Work/lib/mkl/mkl_compact_lnx32e/mkl
MKL_CFLAG = -I${MKL_ROOT}/include
KOKKOS=/home/kyukim/Work/lib/kokkos/install/knl/intel/release
KOKKOS_USE_TPLS=experimental_memkind
CXXFLAGS=-O3 -g -restrict -qopt-report=5 -fp-model fast=2 -D__KOKKOSKERNELS_INTEL_MKL__ -D__KOKKOSKERNELS_INTEL_MKL_BATCHED__ -D__KOKKOSKERNELS_INTEL_MKL_COMPACT_BATCHED__ -DMIMIC_SPARC -DFAST_MATRIX_FILL -D__KOKKOSKERNELS_LIBXSMM__ -I${LIBXSMM_ROOT}/include 
#LIBS = -lm -lnuma
LIBS = ${LIBXSMM_ROOT}/lib/libxsmm.a -lm -lnuma -Wl,--start-group ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.a ${MKL_ROOT}/lib/intel64/libmkl_intel_thread.a ${MKL_ROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 
endif

#### shepard : haswell 2x16x2

ifneq (,$(findstring shepard,$(HOST)))
#CXX=icpc
MKL_CFLAG = -mkl # -DMKL_DIRECT_CALL # -DMKL_DIRECT_CALL_SEQ #:sequential
KOKKOS = /home/kyukim/Work/lib/kokkos/install/hsw/intel/release
KOKKOS_USE_TPLS=experimental_memkind
CXXFLAGS=-O3 -g -restrict -qopt-report=5 -D__KOKKOSKERNELS_INTEL_MKL__ -D__KOKKOSKERNELS_INTEL_MKL_BATCHED__
LIBS = -lm -lnuma
endif

#### bread : kyungjoo local machine

ifneq (,$(findstring bread,$(HOST)))
CXX=icpc
MKL_CFLAG = -mkl #:sequential 
KOKKOS=/home/kyukim/Work/lib/kokkos/install/intel/debug
CXXFLAGS=-O3 -g -restrict -qopt-report=5 -D__KOKKOSKERNELS_INTEL_MKL__ 
LIBS = -lm 
endif
