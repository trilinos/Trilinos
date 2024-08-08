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

#ifndef KOKKOSKERNELS_KOKKOSBATCHED_KERNEL_HEADER_HPP
#define KOKKOSKERNELS_KOKKOSBATCHED_KERNEL_HEADER_HPP

#include <sstream>
#include "KokkosKernels_Error.hpp"

#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
#include <mkl.h>
#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
#include "armpl.h"
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
#include "cuda_runtime.h"
#include "cublas_v2.h"
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
#include <magma_v2.h>
#include <magma_batched.h>
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

namespace KokkosBatched {

/// \brief Heuristic algorithm types. See BatchedKernelHandle for details.
namespace BaseHeuristicAlgos {
enum BASE_HEURISTIC_ALGOS : int { SQUARE = 0, TALL, WIDE, N };
}

/// \brief Tpl algorithm types. See BatchedKernelHandle for details.
namespace BaseTplAlgos {
enum BASE_TPL_ALGOS : int { ARMPL = BaseHeuristicAlgos::N, MKL, N };
}

/// \brief KokkosBatched algorithm types. See BatchedKernelHandle for details.
namespace BaseKokkosBatchedAlgos {
enum BASE_KOKKOS_BATCHED_ALGOS : int { KK_SERIAL = BaseTplAlgos::N, N };
}

#define N_BASE_ALGOS BaseKokkosBatchedAlgos::N
#define BASE_ALGO_STRS                                                                                         \
  "BaseHeuristicAlgos::SQUARE", "BaseHeuristicAlgos::TALL", "BaseHeuristicAlgos::WIDE", "BaseTplAlgos::ARMPL", \
      "BaseTplAlgosMKL", "BaseKokkosBatchedAlgos::KK_SERIAL"

/// \brief TplParams abstracts underlying handle or execution queue type.
struct TplParams {
  union {
#if defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
    // queue mkl_queue;
    // TODO: Add queue header? Cannot find any declarations in intel-18, let
    // alone oneAPI 2021
#endif  // KOKKOSKERNELS_ENABLE_TPL_MKL

#if defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
    armpl_int_t ninter = 1;
#endif  // KOKKOSKERNELS_ENABLE_TPL_ARMPL

#if defined(KOKKOSKERNELS_ENABLE_TPL_CUBLAS)
    cublasHandle_t *cublas_handle;
#endif  // KOKKOSKERNELS_ENABLE_TPL_CUBLAS

#if defined(KOKKOSKERNELS_ENABLE_TPL_MAGMA)
    magma_queue_t *magma_queue;
#endif  // KOKKOSKERNELS_ENABLE_TPL_MAGMA

    uint8_t empty_union;  // suppress empty union build warning
  };
};

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
///                        Note: If the heuristics indicate SIMD views are required for optimal
///                              performance, notify the user that SIMD views are required for
///                              optimal performance.
///    
///                        Specifies which cmake-enabled TPL algorithm to invoke:
///                          TPL algorithms currently SUPPORTED by all batched routines:
///    
///                          TPL algorithms currently UNSUPPORTED by all batched routines:
///                            ARMPL    Invoke the ArmPL TPL interface
///                            MKL      Invoke the MKL TPL interface
///                            CUBLAS   Invoke the CuBLAS TPL interface
///                            MAGMA    Invoke the Magma TPL interface
///                        Note: See each routine's handle for details about which TPLs are supported.
///                        Note: Requires that input views for A, B, and C reside on either host
///                              or device depending on the TPL selected.
///                        Note: If the user selects a TPL, an error will be thrown if:
///                                1. The TPL is not enabled via cmake
///                                2. The input views do not reside on the host/device as needed
///    
///                        Specifies which kokkos-kernels (KK) algorithm to invoke:
///                          KK algorithms currently SUPPORTED   by all KK batched routines:
///                            KK_SERIAL       Invoke SerialFUNC     via RangePolicy(BatchSz)
///                          KK algorithms currently UNSUPPORTED by all KK batched routines:
///                            KK_TEAM         Invoke TeamFUNC       via TeamPolicy(BatchSz)
///                            KK_TEAMVECTOR   Invoke TeamVectorFUNC via TeamPolicy(BatchSz)
///                            KK_SERIALSIMD   Invoke SerialFUNC     via TeamPolicy(BatchSz)
///                            KK_TEAMSIMD     Invoke TeamFUNC       via TeamPolicy(BatchSz)
///                            KK_SERIAL_RANK0 Invoke SerialFUNC     via RangePolicy(BatchSz*LHS_N*LHS_M)
///                            KK_SERIAL_SHMEM Invoke SerialFUNC     via TeamPolicy(BatchSz)
///                            KK_DBLBUF       Solve FUNC            via TeamPolicy(BatchSz*TILES)
///                                            Uses a tuned functor with tiling and double buffering
///                                            via shared memory and register buffers.
///                        Note: See each routine's handle for details about which KK algorithms are
///                              supported.
/// \param teamSz          Specifies the team size that will affect any KK algorithm which uses
///                        TeamPolicy (default, Kokkos::AUTO).
///                        Note: Only applied if useAlgo_type == KK_*
/// \param vecLen          Specifies the vector length that will affect any KK algorithm which
///                        uses TeamPolicy and Kokkos::ThreadVectorRange or Kokkos::TeamVectorRange
///                        (default, Kokkos::AUTO).
///                        Note: Only applied if useAlgo_type == KK_*
/// \param enabledDebug    toggle debug messages.
// clang-format on
class BatchedKernelHandle {
 public:
  int teamSz       = 0;
  int vecLen       = 0;
  bool enableDebug = false;

  BatchedKernelHandle(int kernelAlgoType = BaseHeuristicAlgos::SQUARE, int teamSize = 0, int vecLength = 0)
      : teamSz(teamSize), vecLen(vecLength), _kernelAlgoType(kernelAlgoType) {
#if !defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL) || ARMPL_BUILD < 1058
    if (_kernelAlgoType == BaseTplAlgos::ARMPL) {
      std::ostringstream os;
      os << "KokkosBatched::BatchedKernelHandle requires "
            "KOKKOSKERNELS_ENABLE_TPL_ARMPL and armpl version 21.0.0+"
         << std::endl;
      KokkosKernels::Impl::throw_runtime_exception(os.str());
    }
#endif  // !defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
  };

  int get_kernel_algo_type() const { return _kernelAlgoType; }

  std::string get_kernel_algo_type_str() const { return algo_type_strs[_kernelAlgoType]; }

  decltype(auto) get_tpl_params() const {
#if _kernelAlgoType == ARMPL && defined(KOKKOSKERNELS_ENABLE_TPL_ARMPL)
    return &_tplParamsSingleton.ninter;
#elif _kernelAlgoType == MKL && defined(KOKKOSKERNELS_ENABLE_TPL_MKL)
    return "BaseTplAlgos::MKL does not support any tpl parameters";
#else
    return "Unsupported kernelAlgoType:" + get_kernel_algo_type_str() + ".";
#endif
  }

  // clang-format off
  /// kernelAlgoType: Specifies which algorithm to use for invocation (default, SQUARE).
  /// _tplParams:     a handle or queue specific to the TPL API.
  ///                 managed internally unless provided by user via
  ///                 constructor overload
  // clang-format on
 protected:
  // Define TPL params singleton as static class method variable
  // Only one instance of tplParamsGlobalStorage may exist per process.
  static TplParams &_get_tpl_params_singleton() {
    static TplParams tplParamsGlobalStorage;
    return tplParamsGlobalStorage;
  }

  int _kernelAlgoType            = BaseHeuristicAlgos::SQUARE;
  TplParams &_tplParamsSingleton = _get_tpl_params_singleton();
  bool _tplParamsSet             = false;

 private:
  const char *algo_type_strs[N_BASE_ALGOS] = {BASE_ALGO_STRS};
};

}  // namespace KokkosBatched

#endif  // KOKKOSKERNELS_KOKKOSBATCHED_KERNEL_HEADER_HPP
