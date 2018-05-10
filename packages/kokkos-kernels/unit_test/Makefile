
#=======================================================================
#=================== Settings ==========================================
#=======================================================================

KOKKOS_PATH ?= ${HOME}/Kokkos/kokkos
KOKKOS_SRC_PATH ?= ${KOKKOS_PATH}
KOKKOSKERNELS_PATH ?= ${HOME}/Kokkos/kokkos-kernels
KOKKOSKERNELS_SRC_PATH ?= ${KOKKOSKERNELS_PATH}
GTEST_PATH ?= ${KOKKOS_SRC_PATH}/tpls/gtest

KOKKOS_DEVICES=OpenMP
KOKKOS_CUDA_OPTIONS=enable_lambda

KOKKOSKERNELS_SCALARS="'double,complex<double>'"
KOKKOSKERNELS_OPTIONS=eti-only

default: build

ifneq (,$(findstring Cuda,$(KOKKOS_DEVICES)))
CXX = ${KOKKOS_PATH}/config/nvcc_wrapper
else
CXX = g++
endif

CXXFLAGS = -O3 -g
LINK = ${CXX}
LINKFLAGS =  

DEPFLAGS = -M

OBJ =
LIB = -lpthread


#=======================================================================
#=================== Include Kokkos and KokkosKernels ==================
#=======================================================================

#include $(KOKKOS_PATH)/Makefile.kokkos
include $(KOKKOSKERNELS_PATH)/Makefile.kokkos-kernels

#=======================================================================
#=================== Add common source files ===========================
#=======================================================================

INC  = -I${GTEST_PATH} -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/blas
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/sparse
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/graph
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/common
INC += -I$(KOKKOSKERNELS_SRC_PATH)/unit_test/batched




vpath %.cpp ${KOKKOSKERNELS_SRC_PATH}/unit_test

TEST_HEADERS =  $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/blas/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/sparse/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/graph/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/common/*.hpp)
TEST_HEADERS += $(wildcard $(KOKKOSKERNELS_SRC_PATH)/unit_test/batched/*.hpp)



TEST_TARGETS =
TARGETS =

#=======================================================================
#=================== OpenMP BACKEND ====================================
#=======================================================================

vpath %.cpp ${KOKKOSKERNELS_SRC_PATH}/unit_test/openmp
ifeq ($(KOKKOS_INTERNAL_USE_OPENMP), 1)
  ifeq ($(KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_OPENMP), 1)
    KOKKOSKERNELS_INTERNAL_TEST_OPENMP = 1
  else
    ifeq ($(KOKKOSKERNELS_INTERNAL_ETI_ONLY), 1)
      KOKKOSKERNELS_INTERNAL_TEST_OPENMP = 0
    else
      KOKKOSKERNELS_INTERNAL_TEST_OPENMP = 1
    endif
  endif
endif
ifeq ($(KOKKOSKERNELS_INTERNAL_TEST_OPENMP), 1)
  INC += -I${KOKKOSKERNELS_SRC_PATH}/unit_test/openmp
  OBJ_OPENMP = Test_Main.o gtest-all.o
  OBJ_OPENMP += Test_OpenMP_Blas1_abs.o
  OBJ_OPENMP += Test_OpenMP_Blas1_asum.o
  OBJ_OPENMP += Test_OpenMP_Blas1_axpby.o
  OBJ_OPENMP += Test_OpenMP_Blas1_axpy.o
  OBJ_OPENMP += Test_OpenMP_Blas1_dot.o
  OBJ_OPENMP += Test_OpenMP_Blas1_mult.o
  OBJ_OPENMP += Test_OpenMP_Blas1_nrm1.o
  OBJ_OPENMP += Test_OpenMP_Blas1_nrm2.o
  OBJ_OPENMP += Test_OpenMP_Blas1_nrm2_squared.o
  OBJ_OPENMP += Test_OpenMP_Blas1_nrminf.o
  OBJ_OPENMP += Test_OpenMP_Blas1_reciprocal.o
  OBJ_OPENMP += Test_OpenMP_Blas1_scal.o
  OBJ_OPENMP += Test_OpenMP_Blas1_sum.o
  OBJ_OPENMP += Test_OpenMP_Blas1_update.o
  OBJ_OPENMP += Test_OpenMP_Blas2_gemv.o
  OBJ_OPENMP += Test_OpenMP_Blas3_gemm.o
  OBJ_OPENMP += Test_OpenMP_Sparse_spmv.o
  OBJ_OPENMP += Test_OpenMP_Sparse_trsv.o
  OBJ_OPENMP += Test_OpenMP_Sparse_spgemm.o
  OBJ_OPENMP += Test_OpenMP_Sparse_spadd.o
  OBJ_OPENMP += Test_OpenMP_Sparse_gauss_seidel.o
  OBJ_OPENMP += Test_OpenMP_Sparse_CrsMatrix.o
  OBJ_OPENMP += Test_OpenMP_Sparse_findRelOffset.o
  OBJ_OPENMP += Test_OpenMP_Sparse_replaceSumIntoLonger.o
  OBJ_OPENMP += Test_OpenMP_Sparse_replaceSumInto.o
  OBJ_OPENMP += Test_OpenMP_Graph_graph_color.o
  OBJ_OPENMP += Test_OpenMP_Graph_graph_color_d2.o
  OBJ_OPENMP += Test_OpenMP_Common_ArithTraits.o
  OBJ_OPENMP += Test_OpenMP_Common_set_bit_count.o
#  OBJ_OPENMP += Test_OpenMP_Common_float128.o
 # Real 
  OBJ_OPENMP += Test_OpenMP_Batched_SerialMatUtil_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialGemm_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialTrsm_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialLU_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialGemv_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialTrsv_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamMatUtil_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamGemm_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamTrsm_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamLU_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamGemv_Real.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamTrsv_Real.o
 # Complex
  OBJ_OPENMP += Test_OpenMP_Batched_SerialMatUtil_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialGemm_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialTrsm_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialLU_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialGemv_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_SerialTrsv_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamMatUtil_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamGemm_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamTrsm_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamLU_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamGemv_Complex.o
  OBJ_OPENMP += Test_OpenMP_Batched_TeamTrsv_Complex.o
  # Vector
  OBJ_OPENMP += Test_OpenMP_Batched_VectorArithmatic.o
  OBJ_OPENMP += Test_OpenMP_Batched_VectorMath.o
  OBJ_OPENMP += Test_OpenMP_Batched_VectorRelation.o
  OBJ_OPENMP += Test_OpenMP_Batched_VectorLogical.o
  OBJ_OPENMP += Test_OpenMP_Batched_VectorMisc.o
  OBJ_OPENMP += Test_OpenMP_Batched_VectorView.o
  TARGETS += KokkosKernels_UnitTest_OpenMP
  TEST_TARGETS += test-openmp
endif

KokkosKernels_UnitTest_OpenMP: $(OBJ_OPENMP) $(KOKKOS_LINK_DEPENDS) $(KOKKOSKERNELS_LINK_DEPENDS) $(TEST_HEADERS)
	$(LINK) $(EXTRA_PATH) $(OBJ_OPENMP) $(KOKKOSKERNELS_LDFLAGS) $(KOKKOSKERNELS_LIBS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosKernels_UnitTest_OpenMP

test-openmp: KokkosKernels_UnitTest_OpenMP
	./KokkosKernels_UnitTest_OpenMP

#=======================================================================
#=================== Cuda BACKEND ====================================
#=======================================================================

vpath %.cpp ${KOKKOSKERNELS_SRC_PATH}/unit_test/cuda
ifeq ($(KOKKOS_INTERNAL_USE_CUDA), 1)
  ifeq ($(KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_CUDA), 1)
    KOKKOSKERNELS_INTERNAL_TEST_CUDA = 1
  else
    ifeq ($(KOKKOSKERNELS_INTERNAL_ETI_ONLY), 1)
      KOKKOSKERNELS_INTERNAL_TEST_CUDA = 0
    else
      KOKKOSKERNELS_INTERNAL_TEST_CUDA = 1
    endif
  endif
endif
ifeq ($(KOKKOSKERNELS_INTERNAL_TEST_CUDA), 1)
  INC += -I${KOKKOSKERNELS_SRC_PATH}/unit_test/cuda
  OBJ_CUDA = Test_Main.o gtest-all.o
  OBJ_CUDA += Test_Cuda_Blas1_abs.o
  OBJ_CUDA += Test_Cuda_Blas1_asum.o
  OBJ_CUDA += Test_Cuda_Blas1_axpby.o
  OBJ_CUDA += Test_Cuda_Blas1_axpy.o
  OBJ_CUDA += Test_Cuda_Blas1_dot.o
  OBJ_CUDA += Test_Cuda_Blas1_mult.o
  OBJ_CUDA += Test_Cuda_Blas1_nrm1.o
  OBJ_CUDA += Test_Cuda_Blas1_nrm2.o
  OBJ_CUDA += Test_Cuda_Blas1_nrm2_squared.o
  OBJ_CUDA += Test_Cuda_Blas1_nrminf.o
  OBJ_CUDA += Test_Cuda_Blas1_reciprocal.o
  OBJ_CUDA += Test_Cuda_Blas1_scal.o
  OBJ_CUDA += Test_Cuda_Blas1_sum.o
  OBJ_CUDA += Test_Cuda_Blas1_update.o
  OBJ_CUDA += Test_Cuda_Blas2_gemv.o
  OBJ_CUDA += Test_Cuda_Blas3_gemm.o #Not yet ready need to figure out how to handle CUBLAS
  #OBJ_CUDA += Test_Cuda_Sparse_spmv.o
  #OBJ_CUDA += Test_Cuda_Sparse_trsv.o #removing trsv from cuda unit test as it runs only sequential.
  OBJ_CUDA += Test_Cuda_Sparse_spgemm.o
  OBJ_CUDA += Test_Cuda_Sparse_spadd.o
  OBJ_CUDA += Test_Cuda_Sparse_gauss_seidel.o
  OBJ_CUDA += Test_Cuda_Sparse_CrsMatrix.o
 #OBJ_CUDA += Test_Cuda_Sparse_findRelOffset.o #removing findRelOffset from cuda test as the implementation is sequential.
  OBJ_CUDA += Test_Cuda_Sparse_replaceSumIntoLonger.o
  OBJ_CUDA += Test_Cuda_Sparse_replaceSumInto.o
  OBJ_CUDA += Test_Cuda_Graph_graph_color.o
  OBJ_CUDA += Test_Cuda_Graph_graph_color_d2.o
  OBJ_CUDA += Test_Cuda_Common_ArithTraits.o
  OBJ_CUDA += Test_Cuda_Common_set_bit_count.o
  # Real
  OBJ_CUDA += Test_Cuda_Batched_SerialMatUtil_Real.o
  OBJ_CUDA += Test_Cuda_Batched_SerialGemm_Real.o
  OBJ_CUDA += Test_Cuda_Batched_SerialTrsm_Real.o
  OBJ_CUDA += Test_Cuda_Batched_SerialLU_Real.o
  OBJ_CUDA += Test_Cuda_Batched_SerialGemv_Real.o
  OBJ_CUDA += Test_Cuda_Batched_SerialTrsv_Real.o
  OBJ_CUDA += Test_Cuda_Batched_TeamMatUtil_Real.o
  OBJ_CUDA += Test_Cuda_Batched_TeamGemm_Real.o
  OBJ_CUDA += Test_Cuda_Batched_TeamTrsm_Real.o
  OBJ_CUDA += Test_Cuda_Batched_TeamLU_Real.o
  OBJ_CUDA += Test_Cuda_Batched_TeamGemv_Real.o
  OBJ_CUDA += Test_Cuda_Batched_TeamTrsv_Real.o
  # Complex
  OBJ_CUDA += Test_Cuda_Batched_SerialMatUtil_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_SerialGemm_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_SerialTrsm_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_SerialLU_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_SerialGemv_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_SerialTrsv_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_TeamMatUtil_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_TeamGemm_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_TeamTrsm_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_TeamLU_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_TeamGemv_Complex.o
  OBJ_CUDA += Test_Cuda_Batched_TeamTrsv_Complex.o
  TARGETS += KokkosKernels_UnitTest_Cuda
  TEST_TARGETS += test-cuda
endif

KokkosKernels_UnitTest_Cuda: $(OBJ_CUDA) $(KOKKOS_LINK_DEPENDS) $(KOKKOSKERNELS_LINK_DEPENDS) $(TEST_HEADERS)
	$(LINK) $(EXTRA_PATH) $(OBJ_CUDA) $(KOKKOSKERNELS_LDFLAGS) $(KOKKOSKERNELS_LIBS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosKernels_UnitTest_Cuda

test-cuda: KokkosKernels_UnitTest_Cuda
	./KokkosKernels_UnitTest_Cuda

#=======================================================================
#=================== Serial BACKEND ====================================
#=======================================================================

vpath %.cpp ${KOKKOSKERNELS_SRC_PATH}/unit_test/serial
ifeq ($(KOKKOS_INTERNAL_USE_SERIAL), 1)
  ifeq ($(KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_SERIAL), 1)
    KOKKOSKERNELS_INTERNAL_TEST_SERIAL = 1
  else
    ifeq ($(KOKKOSKERNELS_INTERNAL_ETI_ONLY), 1)
      KOKKOSKERNELS_INTERNAL_TEST_SERIAL = 0
    else
      KOKKOSKERNELS_INTERNAL_TEST_SERIAL = 1
    endif
  endif
endif

ifeq ($(KOKKOSKERNELS_INTERNAL_TEST_SERIAL), 1)
  INC += -I${KOKKOSKERNELS_SRC_PATH}/unit_test/serial
  OBJ_SERIAL = Test_Main.o gtest-all.o
  OBJ_SERIAL += Test_Serial_Blas1_abs.o
  OBJ_SERIAL += Test_Serial_Blas1_asum.o
  OBJ_SERIAL += Test_Serial_Blas1_axpby.o
  OBJ_SERIAL += Test_Serial_Blas1_axpy.o
  OBJ_SERIAL += Test_Serial_Blas1_dot.o
  OBJ_SERIAL += Test_Serial_Blas1_mult.o
  OBJ_SERIAL += Test_Serial_Blas1_nrm1.o
  OBJ_SERIAL += Test_Serial_Blas1_nrm2.o
  OBJ_SERIAL += Test_Serial_Blas1_nrm2_squared.o
  OBJ_SERIAL += Test_Serial_Blas1_nrminf.o
  OBJ_SERIAL += Test_Serial_Blas1_reciprocal.o
  OBJ_SERIAL += Test_Serial_Blas1_scal.o
  OBJ_SERIAL += Test_Serial_Blas1_sum.o
  OBJ_SERIAL += Test_Serial_Blas1_update.o
  OBJ_SERIAL += Test_Serial_Blas2_gemv.o
  OBJ_SERIAL += Test_Serial_Blas3_gemm.o
  OBJ_SERIAL += Test_Serial_Sparse_spmv.o
  OBJ_SERIAL += Test_Serial_Sparse_trsv.o
  OBJ_SERIAL += Test_Serial_Sparse_spgemm.o
  OBJ_SERIAL += Test_Serial_Sparse_spadd.o
  OBJ_SERIAL += Test_Serial_Sparse_gauss_seidel.o
  OBJ_SERIAL += Test_Serial_Sparse_CrsMatrix.o
  OBJ_SERIAL += Test_Serial_Sparse_findRelOffset.o
  OBJ_SERIAL += Test_Serial_Sparse_replaceSumIntoLonger.o
  OBJ_SERIAL += Test_Serial_Sparse_replaceSumInto.o
  OBJ_SERIAL += Test_Serial_Graph_graph_color.o
  OBJ_SERIAL += Test_Serial_Graph_graph_color_d2.o
  OBJ_SERIAL += Test_Serial_Common_ArithTraits.o
  OBJ_SERIAL += Test_Serial_Common_set_bit_count.o
#  OBJ_SERIAL += Test_Serial_Common_float128.o
  # Real
  OBJ_SERIAL += Test_Serial_Batched_SerialMatUtil_Real.o
  OBJ_SERIAL += Test_Serial_Batched_SerialGemm_Real.o
  OBJ_SERIAL += Test_Serial_Batched_SerialTrsm_Real.o
  OBJ_SERIAL += Test_Serial_Batched_SerialLU_Real.o
  OBJ_SERIAL += Test_Serial_Batched_SerialGemv_Real.o
  OBJ_SERIAL += Test_Serial_Batched_SerialTrsv_Real.o
  OBJ_SERIAL += Test_Serial_Batched_TeamMatUtil_Real.o
  OBJ_SERIAL += Test_Serial_Batched_TeamGemm_Real.o
  OBJ_SERIAL += Test_Serial_Batched_TeamTrsm_Real.o
  OBJ_SERIAL += Test_Serial_Batched_TeamLU_Real.o
  OBJ_SERIAL += Test_Serial_Batched_TeamGemv_Real.o
  OBJ_SERIAL += Test_Serial_Batched_TeamTrsv_Real.o
  # Complex
  OBJ_SERIAL += Test_Serial_Batched_SerialMatUtil_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_SerialGemm_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_SerialTrsm_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_SerialLU_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_SerialGemv_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_SerialTrsv_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_TeamMatUtil_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_TeamGemm_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_TeamTrsm_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_TeamLU_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_TeamGemv_Complex.o
  OBJ_SERIAL += Test_Serial_Batched_TeamTrsv_Complex.o
  # Vector
  OBJ_SERIAL += Test_Serial_Batched_VectorArithmatic.o
  OBJ_SERIAL += Test_Serial_Batched_VectorMath.o
  OBJ_SERIAL += Test_Serial_Batched_VectorRelation.o
  OBJ_SERIAL += Test_Serial_Batched_VectorLogical.o
  OBJ_SERIAL += Test_Serial_Batched_VectorMisc.o
  OBJ_SERIAL += Test_Serial_Batched_VectorView.o
  TARGETS += KokkosKernels_UnitTest_Serial
  TEST_TARGETS += test-serial
endif

KokkosKernels_UnitTest_Serial: $(OBJ_SERIAL) $(KOKKOS_LINK_DEPENDS) $(KOKKOSKERNELS_LINK_DEPENDS) $(TEST_HEADERS)
	$(LINK) $(EXTRA_PATH) $(OBJ_SERIAL) $(KOKKOSKERNELS_LDFLAGS) $(KOKKOSKERNELS_LIBS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosKernels_UnitTest_Serial

test-serial: KokkosKernels_UnitTest_Serial
	./KokkosKernels_UnitTest_Serial

#=======================================================================
#=================== Threads BACKEND ===================================
#=======================================================================

vpath %.cpp ${KOKKOSKERNELS_SRC_PATH}/unit_test/threads
ifeq ($(KOKKOS_INTERNAL_USE_PTHREADS), 1)
  ifeq ($(KOKKOSKERNELS_INTERNAL_INST_EXECSPACE_THREADS), 1)
    KOKKOSKERNELS_INTERNAL_TEST_THREADS = 1
  else
    ifeq ($(KOKKOSKERNELS_INTERNAL_ETI_ONLY), 1)
      KOKKOSKERNELS_INTERNAL_TEST_THREADS = 0
    else
      KOKKOSKERNELS_INTERNAL_TEST_THREADS = 1
    endif
  endif
endif
ifeq ($(KOKKOSKERNELS_INTERNAL_TEST_THREADS), 1)
  INC += -I${KOKKOSKERNELS_SRC_PATH}/unit_test/threads
  OBJ_THREADS = Test_Main.o gtest-all.o
  OBJ_THREADS += Test_Threads_Blas1_abs.o
  OBJ_THREADS += Test_Threads_Blas1_asum.o
  OBJ_THREADS += Test_Threads_Blas1_axpby.o
  OBJ_THREADS += Test_Threads_Blas1_axpy.o
  OBJ_THREADS += Test_Threads_Blas1_dot.o
  OBJ_THREADS += Test_Threads_Blas1_mult.o
  OBJ_THREADS += Test_Threads_Blas1_nrm1.o
  OBJ_THREADS += Test_Threads_Blas1_nrm2.o
  OBJ_THREADS += Test_Threads_Blas1_nrm2_squared.o
  OBJ_THREADS += Test_Threads_Blas1_nrminf.o
  OBJ_THREADS += Test_Threads_Blas1_reciprocal.o
  OBJ_THREADS += Test_Threads_Blas1_scal.o
  OBJ_THREADS += Test_Threads_Blas1_sum.o
  OBJ_THREADS += Test_Threads_Blas1_update.o
  OBJ_THREADS += Test_Threads_Blas2_gemv.o
  OBJ_THREADS += Test_Threads_Blas3_gemm.o
  OBJ_THREADS += Test_Threads_Sparse_spmv.o
  OBJ_THREADS += Test_Threads_Sparse_trsv.o
  OBJ_THREADS += Test_Threads_Sparse_spgemm.o
  OBJ_THREADS += Test_Threads_Sparse_spadd.o
  OBJ_THREADS += Test_Threads_Sparse_gauss_seidel.o
  OBJ_THREADS += Test_Threads_Sparse_findRelOffset.o
  OBJ_THREADS += Test_Threads_Sparse_replaceSumIntoLonger.o
  OBJ_THREADS += Test_Threads_Sparse_replaceSumInto.o
  OBJ_THREADS += Test_Threads_Sparse_CrsMatrix.o
  OBJ_THREADS += Test_Threads_Graph_graph_color.o
  OBJ_THREADS += Test_Threads_Graph_graph_color_d2.o
  OBJ_THREADS += Test_Threads_Common_ArithTraits.o
  OBJ_THREADS += Test_Threads_Common_set_bit_count.o
#  OBJ_THREADS += Test_Threads_Common_float128.o
  TARGETS += KokkosKernels_UnitTest_Threads
  TEST_TARGETS += test-threads
endif

KokkosKernels_UnitTest_Threads: $(OBJ_THREADS) $(KOKKOS_LINK_DEPENDS) $(KOKKOSKERNELS_LINK_DEPENDS) $(TEST_HEADERS)
	$(LINK) $(EXTRA_PATH) $(OBJ_THREADS) $(KOKKOSKERNELS_LDFLAGS) $(KOKKOSKERNELS_LIBS) $(KOKKOS_LIBS) $(LIB) $(KOKKOS_LDFLAGS) $(LDFLAGS) -o KokkosKernels_UnitTest_Threads

test-threads: KokkosKernels_UnitTest_Threads
	./KokkosKernels_UnitTest_Threads

#=======================================================================
#=================== Make Rules ========================================
#=======================================================================

test: $(TEST_TARGETS)

build: $(TARGETS)

clean: kokkos-clean kokkoskernels-clean
	rm -f *.o

# Compilation rules

%.o:%.cpp $(KOKKOS_CPP_DEPENDS) $(KOKKOSKERNELS_CPP_DEPENDS) $(TEST_HEADERS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOSKERNELS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(INC) $(CXXFLAGS) $(EXTRA_INC) -c $<

gtest-all.o:$(GTEST_PATH)/gtest/gtest-all.cc
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(INC) $(CXXFLAGS) $(EXTRA_INC) -c $(GTEST_PATH)/gtest/gtest-all.cc
