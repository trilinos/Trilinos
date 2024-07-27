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
#ifndef __TACHO_HERK_ON_DEVICE_HPP__
#define __TACHO_HERK_ON_DEVICE_HPP__

/// \file  Tacho_Herk_OnDevice.hpp
/// \brief BLAS hermitian rank-k update
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {

template <typename ArgUplo, typename ArgTrans> struct Herk<ArgUplo, ArgTrans, Algo::OnDevice> {
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeC>
  inline static int blas_invoke(const ScalarType alpha, const ViewTypeA &A, const ScalarType beta, const ViewTypeC &C) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type n = C.extent(0),
                       k = (std::is_same<ArgTrans, Trans::NoTranspose>::value ? A.extent(1) : A.extent(0));

    if (n > 0 && k > 0) {
      Blas<value_type>::herk(ArgUplo::param, ArgTrans::param, n, k, value_type(alpha), A.data(), A.stride_1(),
                             value_type(beta), C.data(), C.stride_1());
    }
    return 0;
  }
#if defined(KOKKOS_ENABLE_CUDA)
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeC>
  inline static int cublas_invoke(cublasHandle_t &handle, const ScalarType alpha, const ViewTypeA &A,
                                  const ScalarType beta, const ViewTypeC &C) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type n = C.extent(0),
                       k = (std::is_same<ArgTrans, Trans::NoTranspose>::value ? A.extent(1) : A.extent(0));

    int r_val(0);
    if (n > 0 && k > 0) {
      if (std::is_same<value_type, float>::value || std::is_same<value_type, double>::value)
        r_val =
            Blas<value_type>::herk(handle, ArgUplo::cublas_param,
                                   std::is_same<ArgTrans, Trans::ConjTranspose>::value ? Trans::Transpose::cublas_param
                                                                                       : ArgTrans::cublas_param,
                                   n, k, alpha, A.data(), A.stride_1(), beta, C.data(), C.stride_1());
      else if (std::is_same<value_type, Kokkos::complex<float>>::value ||
               std::is_same<value_type, Kokkos::complex<double>>::value)
        r_val = Blas<value_type>::herk(handle, ArgUplo::cublas_param, ArgTrans::cublas_param, n, k, alpha, A.data(),
                                       A.stride_1(), beta, C.data(), C.stride_1());
    }
    return r_val;
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  template <typename ScalarType, typename ViewTypeA, typename ViewTypeC>
  inline static int rocblas_invoke(rocblas_handle &handle, const ScalarType alpha, const ViewTypeA &A,
                                   const ScalarType beta, const ViewTypeC &C) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type n = C.extent(0),
                       k = (std::is_same<ArgTrans, Trans::NoTranspose>::value ? A.extent(1) : A.extent(0));

    int r_val(0);
    if (n > 0 && k > 0) {
      if (std::is_same<value_type, float>::value || std::is_same<value_type, double>::value)
        r_val =
            Blas<value_type>::herk(handle, ArgUplo::rocblas_param,
                                   std::is_same<ArgTrans, Trans::ConjTranspose>::value ? Trans::Transpose::rocblas_param
                                                                                       : ArgTrans::rocblas_param,
                                   n, k, alpha, A.data(), A.stride_1(), beta, C.data(), C.stride_1());
      else if (std::is_same<value_type, Kokkos::complex<float>>::value ||
               std::is_same<value_type, Kokkos::complex<double>>::value)
        r_val = Blas<value_type>::herk(handle, ArgUplo::rocblas_param, ArgTrans::rocblas_param, n, k, alpha, A.data(),
                                       A.stride_1(), beta, C.data(), C.stride_1());
    }
    return r_val;
  }
#endif

  template <typename MemberType, typename ScalarType, typename ViewTypeA, typename ViewTypeC>
  inline static int invoke(MemberType &member, const ScalarType alpha, const ViewTypeA &A, const ScalarType beta,
                           const ViewTypeC &C) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeC::non_const_value_type value_type_c;

    typedef typename ViewTypeA::memory_space memory_space;
    typedef typename ViewTypeC::memory_space memory_space_c;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeC::rank == 2, "C is not rank 2 view.");

    static_assert(std::is_same<value_type, value_type_c>::value, "A and C do not have the same value type.");

    static_assert(std::is_same<memory_space, memory_space_c>::value, "A and C do not have the same memory space.");
    int r_val(0);
    if (std::is_same<memory_space, Kokkos::HostSpace>::value)
      r_val = blas_invoke(alpha, A, beta, C);
#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<memory_space, Kokkos::CudaSpace>::value || std::is_same<memory_space, Kokkos::CudaUVMSpace>::value)
      r_val = cublas_invoke(member, alpha, A, beta, C);
#endif
#if defined(KOKKOS_ENABLE_HIP)
    if (std::is_same<memory_space, Kokkos::HIPSpace>::value)
      r_val = rocblas_invoke(member, alpha, A, beta, C);
#endif
    return r_val;
  }
};

} // namespace Tacho
#endif
