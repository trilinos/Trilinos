#### common

INCLUDES= -I../

#### hasen : cuda k80

ifneq (,$(findstring hansen,$(HOST)))
#CXX=$(NVCC_WRAPPER)
MKL_CFLAG= 
KOKKOS=/home/kyukim/Work/lib/kokkos/install/cuda/gnu/release
CXXFLAGS=-O3 -g 
LIBS=
endif

#### bowman : knl 1x68x4

ifneq (,$(findstring bowman,$(HOST)))
#CXX=icpc
MKL_CFLAG = -mkl #:sequential
KOKKOS=/home/kyukim/Work/lib/kokkos/install/knl/intel/release
KOKKOS_USE_TPLS=experimental_memkind
CXXFLAGS=-O3 -g -restrict -qopt-report=5 -D__KOKKOSKERNELS_INTEL_MKL__ -D__KOKKOSKERNELS_INTEL_MKL_BATCHED__
LIBS = -lm -lnuma
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
KOKKOS=/home/kyukim/Work/lib/kokkos/install/intel/release
CXXFLAGS=-O3 -g -restrict -qopt-report=5 -D__KOKKOSKERNELS_INTEL_MKL__ 
LIBS = -lm 
endif
