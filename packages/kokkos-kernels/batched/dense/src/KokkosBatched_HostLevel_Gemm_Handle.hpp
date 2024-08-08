//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef __KOKKOSBATCHED_HOSTLEVEL_GEMM_HANDLE_DECL_HPP__
#define __KOKKOSBATCHED_HOSTLEVEL_GEMM_HANDLE_DECL_HPP__

#include "KokkosBatched_Kernel_Handle.hpp"

namespace KokkosBatched {

/// \brief Tpl algorithm types. See BatchedGemmHandle for details.
namespace GemmTplAlgos {
enum GEMM_TPL_ALGOS : int { CUBLAS = N_BASE_ALGOS, MAGMA, N };
}

/// \brief KokkosBatched algorithm types. See BatchedGemmHandle for details.
namespace GemmKokkosBatchedAlgos {
enum GEMM_KOKKOS_BATCHED_ALGOS : int {
  KK_TEAM = GemmTplAlgos::N,
  KK_TEAMVECTOR,
  KK_SERIALSIMD,
  KK_TEAMSIMD,
  KK_SERIAL_RANK0,
  KK_SERIAL_SHMEM,
  KK_DBLBUF,
  N
};
}

#define GEMM_ALGO_STRS                                                                  \
  "GemmTplAlgos::CUBLAS", "GemmTplAlgos::MAGMA", "GemmKokkosBatchedAlgos::KK_TEAM",     \
      "GemmKokkosBatchedAlgos::KK_TEAMVECTOR", "GemmKokkosBatchedAlgos::KK_SERIALSIMD", \
      "GemmKokkosBatchedAlgos::KK_TEAMSIMD", "GemmKokkosBatchedAlgos::KK_SERIAL_RANK0", \
      "GemmKokkosBatchedAlgos::KK_SERIAL_SHMEM", "GemmKokkosBatchedAlgos::KK_DBLBUF"
// clang-format off
/// \brief Handle for selecting runtime behavior of the BatchedGemm interface.
///
/// \param kernelAlgoType  Specifies which algorithm to use for invocation (default, SQUARE).
///
///                        Specifies whether to select optimal invocations based on inputs and
///                        heuristics:
///                          SQUARE select invocations based on square matrix heuristics where M=N
///                          TALL   select invocations based on tall   matrix heuristics where M>N
///                          WIDE   select invocations based on wide   matrix heuristics where M<N
///    
///                        Specifies which cmake-enabled TPL algorithm to invoke:
///                          ARMPL    Invoke the ArmPL TPL interface  (Currently UNSUPPORTED)
///                          MKL      Invoke the MKL TPL interface    (Currently UNSUPPORTED)
///                          CUBLAS   Invoke the CuBLAS TPL interface (Currently UNSUPPORTED)
///                          MAGMA    Invoke the Magma TPL interface  (Currently UNSUPPORTED)
///                        Note: Requires that input views for A, B, and C reside on either host
///                              or device depending on the TPL selected.
///                        Note: If the user selects a TPL, an error will be thrown if:
///                                1. The TPL is not enabled via cmake
///                                2. The input views do not reside on the host/device as needed
///    
///                        Specifies which kokkos-kernels (KK) algorithm to invoke:
///                          KK_SERIAL       Invoke SerialGemm     via RangePolicy(BatchSz)
///                          KK_TEAM         Invoke TeamGemm       via TeamPolicy(BatchSz)
///                          KK_TEAMVECTOR   Invoke TeamVectorGemm via TeamPolicy(BatchSz)
///                          KK_SERIALSIMD   Invoke SerialGemm     via TeamPolicy(BatchSz)
///                          KK_TEAMSIMD     Invoke TeamGemm       via TeamPolicy(BatchSz)
///                          KK_SERIAL_RANK0 Invoke SerialGemm     via RangePolicy(BatchSz*N*M)
///                                          Each thread computes one element of C.
///                          KK_SERIAL_SHMEM Invoke SerialGemm     via TeamPolicy(BatchSz)
///                                          Copies A and B to shared memory before GEMM.
///                                          Each vector lane solves one element of C via SerialGemm.
///                          KK_DBLBUF       Solve GEMM            via TeamPolicy(BatchSz*TILES)
///                                          Uses a tuned functor with tiling and double buffering
///                                          via shared memory and register buffers.
///                                          KK_DBLBUF generally performs better on GPUs when M, N >= 24.
/// \param teamSz      Specifies the team size that will affect any KK algorithm which uses
///                    TeamPolicy (default, Kokkos::AUTO).
///                    Note: Only applied if useAlgo_type == KK_*
/// \param vecLen      Specifies the vector length that will affect any KK algorithm which
///                    uses TeamPolicy and Kokkos::ThreadVectorRange or Kokkos::TeamVectorRange
///                    (default, Kokkos::AUTO).
///                    Note: Only applied if useAlgo_type == KK_*
// clang-format on
class BatchedGemmHandle : public BatchedKernelHandle {
 public:
  BatchedGemmHandle(int kernelAlgoType = BaseHeuristicAlgos::SQUARE, int teamSize = 0, int vecLength = 0)
      : BatchedKernelHandle(kernelAlgoType, teamSize, vecLength) {
#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    if (!_tplParamsSet && kernelAlgoType == GemmTplAlgos::CUBLAS) {
      static cublasHandle_t cublas_handle;
      _tplParamsSingleton.cublas_handle = &cublas_handle;
      _tplParamsSet                     = true;
    }
#endif  // CUBLAS

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
    if (!_tplParamsSet && kernelAlgoType == GemmTplAlgos::MAGMA) {
      static magma_queue_t magma_queue;
      _tplParamsSingleton.magma_queue = &magma_queue;
      _tplParamsSet                   = true;
    }
#endif  // MAGMA
  };

  BatchedGemmHandle(bool tplParamsSet, int kernelAlgoType = BaseHeuristicAlgos::SQUARE, int teamSize = 0,
                    int vecLength = 0)
      : BatchedKernelHandle(kernelAlgoType, teamSize, vecLength) {
    _tplParamsSet = tplParamsSet;
  };

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
  BatchedGemmHandle(cublasHandle_t &cublas_handle, int kernelAlgoType = BaseHeuristicAlgos::SQUARE, int teamSize = 0,
                    int vecLength = 0)
      : BatchedGemmHandle(true, kernelAlgoType, teamSize, vecLength) {
    _tplParamsSingleton.cublas_handle = &cublas_handle;
  };
#endif  // CUBLAS

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
  BatchedGemmHandle(magma_queue_t &magma_queue, int kernelAlgoType = BaseHeuristicAlgos::SQUARE, int teamSize = 0,
                    int vecLength = 0)
      : BatchedGemmHandle(true, kernelAlgoType, teamSize, vecLength) {
    _tplParamsSingleton.magma_queue = &magma_queue;
  };
#endif  // MAGMA

  decltype(auto) get_tpl_params() {
#if _kernelAlgoType == CUBLAS && defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    return _tplParamsSingleton.cublas_handle;
#elif _kernelAlgoType == MAGMA && defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
    return _tplParamsSingleton.magma_queue;
#else
    return this->BatchedKernelHandle::get_tpl_params();
#endif
  }

  std::string get_kernel_algo_type_str() const { return gemm_algo_type_strs[_kernelAlgoType]; }

 private:
  const char *gemm_algo_type_strs[GemmKokkosBatchedAlgos::N] = {BASE_ALGO_STRS, GEMM_ALGO_STRS};
};

}  // namespace KokkosBatched

#endif  // __KOKKOSBATCHED_HOSTLEVEL_GEMM_HANDLE_DECL_HPP__
