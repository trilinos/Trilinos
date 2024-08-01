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
#ifndef __TACHO_CHOL_ON_DEVICE_HPP__
#define __TACHO_CHOL_ON_DEVICE_HPP__

/// \file  Tacho_Chol_OnDevice.hpp
/// \brief BLAS general matrix matrix multiplication
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <typename ArgUplo> struct Chol<ArgUplo, Algo::OnDevice> {
  template <typename ViewTypeA> inline static int lapack_invoke(const ViewTypeA &A) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    if (m > 0) {
      Lapack<value_type>::potrf(ArgUplo::param, m, A.data(), A.stride_1(), &r_val);
    }
    return r_val;
  }

#if defined(KOKKOS_ENABLE_CUDA)
  template <typename ViewTypeA, typename ViewTypeW>
  inline static int cusolver_invoke(cusolverDnHandle_t &handle, const ViewTypeA &A, const ViewTypeW &W) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeW::non_const_value_type work_value_type;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    if (m > 0) {
      int *devInfo = (int *)W.data();
      value_type *workspace = W.data() + 1;
      int lwork = W.span() - 1;
      r_val = Lapack<value_type>::potrf(handle, ArgUplo::cublas_param, m, A.data(), A.stride_1(), workspace, lwork,
                                        devInfo);
    }
    return r_val;
  }

  template <typename ViewTypeA>
  inline static int cusolver_buffer_size(cusolverDnHandle_t &handle, const ViewTypeA &A, int *lwork) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    if (m > 0)
      r_val = Lapack<value_type>::potrf_buffersize(handle, ArgUplo::cublas_param, m, A.data(), A.stride_1(), lwork);
    return r_val;
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  template <typename ViewTypeA, typename ViewTypeW>
  inline static int rocsolver_invoke(rocblas_handle &handle, const ViewTypeA &A, const ViewTypeW &W) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    if (m > 0) {
      int *devInfo = (int *)W.data();
      r_val = Lapack<value_type>::potrf(handle, ArgUplo::rocblas_param, m, A.data(), A.stride_1(), devInfo);
    }
    return r_val;
  }
#endif

  template <typename MemberType, typename ViewTypeA, typename ViewTypeW>
  inline static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeW &W) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeW::non_const_value_type value_type_w;

    typedef typename ViewTypeA::memory_space memory_space;
    typedef typename ViewTypeW::memory_space memory_space_w;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeW::rank == 1, "W is not rank 1 view.");

    static_assert(std::is_same<value_type, value_type_w>::value, "A and W do not have the same value type.");

    static_assert(std::is_same<memory_space, memory_space_w>::value, "A and W do not have the same memory space.");
    int r_val(0);
    if (std::is_same<memory_space, Kokkos::HostSpace>::value) {
      r_val = lapack_invoke(A);
    }

#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<memory_space, Kokkos::CudaSpace>::value ||
        std::is_same<memory_space, Kokkos::CudaUVMSpace>::value) {
      if (W.span() == 0) {
        int lwork;
        r_val = cusolver_buffer_size(member, A, &lwork);
        r_val = lwork + 1;
      } else
        r_val = cusolver_invoke(member, A, W);
    }
#endif

#if defined(KOKKOS_ENABLE_HIP)
    if (std::is_same<memory_space, Kokkos::HIPSpace>::value) {
      if (W.span() == 0) {
        r_val = 2;
      } else
        r_val = rocsolver_invoke(member, A, W);
    }
#endif

    return r_val;
  }
};

} // namespace Tacho
#endif
