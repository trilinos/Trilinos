// clang-format off
// @HEADER
// *****************************************************************************
//                            Tacho package
//
// Copyright 2022 NTESS and the Tacho contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
// clang-format on
#ifndef __TACHO_GEMM_ON_DEVICE_HPP__
#define __TACHO_GEMM_ON_DEVICE_HPP__

/// \file  Tacho_Gemm_OnDevice.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {

template <typename ArgTransA, typename ArgTransB> struct Gemm<ArgTransA, ArgTransB, Algo::OnDevice> {
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int blas_invoke(const ScalarType alpha, const ViewTypeA &A, const ViewTypeB &B, const ScalarType beta,
                                const ViewTypeC &C) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = C.extent(0), n = C.extent(1),
                       k = (std::is_same<ArgTransB, Trans::NoTranspose>::value ? B.extent(0) : B.extent(1));

    if (m > 0 && n > 0 && k > 0) {
      Blas<value_type>::gemm(ArgTransA::param, ArgTransB::param, m, n, k, value_type(alpha), A.data(), A.stride_1(),
                             B.data(), B.stride_1(), value_type(beta), C.data(), C.stride_1());
    }
    return 0;
  }

#if defined(KOKKOS_ENABLE_CUDA)
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int cublas_invoke(cublasHandle_t &handle, const ScalarType alpha, const ViewTypeA &A,
                                  const ViewTypeB &B, const ScalarType beta, const ViewTypeC &C) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    const ordinal_type m = C.extent(0), n = C.extent(1),
                       k = (std::is_same<ArgTransB, Trans::NoTranspose>::value ? B.extent(0) : B.extent(1));

    int r_val(0);
    if (m > 0 && n > 0 && k > 0) {
      r_val = Blas<value_type>::gemm(handle, ArgTransA::cublas_param, ArgTransB::cublas_param, m, n, k, alpha, A.data(),
                                     A.stride_1(), B.data(), B.stride_1(), beta, C.data(), C.stride_1());
    }
    return r_val;
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int rocblas_invoke(rocblas_handle &handle, const ScalarType alpha, const ViewTypeA &A,
                                   const ViewTypeB &B, const ScalarType beta, const ViewTypeC &C) {
    typedef typename ViewTypeA::non_const_value_type value_type;

    const ordinal_type m = C.extent(0), n = C.extent(1),
                       k = (std::is_same<ArgTransB, Trans::NoTranspose>::value ? B.extent(0) : B.extent(1));

    int r_val(0);
    if (m > 0 && n > 0 && k > 0) {
      r_val = Blas<value_type>::gemm(handle, ArgTransA::rocblas_param, ArgTransB::rocblas_param, m, n, k, alpha,
                                     A.data(), A.stride_1(), B.data(), B.stride_1(), beta, C.data(), C.stride_1());
    }
    return r_val;
  }
#endif

  template <typename MemberType, typename ScalarType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int invoke(MemberType &member, const ScalarType alpha, const ViewTypeA &A, const ViewTypeB &B,
                           const ScalarType beta, const ViewTypeC &C) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeB::non_const_value_type value_type_b;
    typedef typename ViewTypeC::non_const_value_type value_type_c;

    typedef typename ViewTypeA::memory_space memory_space;
    typedef typename ViewTypeB::memory_space memory_space_b;
    typedef typename ViewTypeC::memory_space memory_space_c;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeB::rank == 2, "B is not rank 2 view.");
    static_assert(ViewTypeC::rank == 2, "C is not rank 2 view.");

    static_assert(std::is_same<value_type, value_type_b>::value && std::is_same<value_type_b, value_type_c>::value,
                  "A, B and C do not have the same value type.");

    static_assert(std::is_same<memory_space, memory_space_b>::value &&
                      std::is_same<memory_space_b, memory_space_c>::value,
                  "A, B and C do not have the same memory space.");
    int r_val(0);
    if (std::is_same<memory_space, Kokkos::HostSpace>::value)
      r_val = blas_invoke(alpha, A, B, beta, C);
#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<memory_space, Kokkos::CudaSpace>::value || std::is_same<memory_space, Kokkos::CudaUVMSpace>::value)
      r_val = cublas_invoke(member, alpha, A, B, beta, C);
#endif
#if defined(KOKKOS_ENABLE_HIP)
    if (std::is_same<memory_space, Kokkos::HIPSpace>::value)
      r_val = rocblas_invoke(member, alpha, A, B, beta, C);
#endif
    return r_val;
  }
};

} // namespace Tacho
#endif
