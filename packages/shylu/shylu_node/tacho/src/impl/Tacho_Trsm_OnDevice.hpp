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
#ifndef __TACHO_TRSM_ON_DEVICE_HPP__
#define __TACHO_TRSM_ON_DEVICE_HPP__

/// \file  Tacho_Trsm_OnDevice.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

template <typename ArgSide, typename ArgUplo, typename ArgTransA>
struct Trsm<ArgSide, ArgUplo, ArgTransA, Algo::OnDevice> {

  template <typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB>
  inline static int blas_invoke(const DiagType diagA, const ScalarType alpha, const ViewTypeA &A, const ViewTypeB &B) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = B.extent(0);
    const ordinal_type n = B.extent(1);

    if (m > 0 && n > 0)
      Blas<value_type>::trsm(ArgSide::param, ArgUplo::param, ArgTransA::param, diagA.param, m, n, value_type(alpha),
                             A.data(), A.stride(1), B.data(), B.stride(1));
    return 0;
  }

#if defined(KOKKOS_ENABLE_CUDA)
  template <typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB>
  inline static int cublas_invoke(cublasHandle_t &handle, const DiagType diagA, const ScalarType alpha,
                                  const ViewTypeA &A, const ViewTypeB &B) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = B.extent(0);
    const ordinal_type n = B.extent(1);

    int r_val(0);
    if (m > 0 && n > 0)
      Blas<value_type>::trsm(handle, ArgSide::cublas_param, ArgUplo::cublas_param, ArgTransA::cublas_param,
                             diagA.cublas_param, m, n, alpha, A.data(), A.stride(1), B.data(), B.stride(1));
    return r_val;
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  template <typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB>
  inline static int rocblas_invoke(rocblas_handle &handle, const DiagType diagA, const ScalarType alpha,
                                   const ViewTypeA &A, const ViewTypeB &B) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = B.extent(0);
    const ordinal_type n = B.extent(1);

    int r_val(0);
    if (m > 0 && n > 0)
      Blas<value_type>::trsm(handle, ArgSide::rocblas_param, ArgUplo::rocblas_param, ArgTransA::rocblas_param,
                             diagA.rocblas_param, m, n, alpha, A.data(), A.stride(1), B.data(), B.stride(1));
    return r_val;
  }
#endif

  template <typename MemberType, typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB>
  inline static int invoke(MemberType &member, const DiagType diagA, const ScalarType alpha, const ViewTypeA &A,
                           const ViewTypeB &B) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeB::non_const_value_type value_type_b;

    typedef typename ViewTypeA::memory_space memory_space;
    typedef typename ViewTypeB::memory_space memory_space_b;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeB::rank == 2, "B is not rank 2 view.");

    static_assert(std::is_same<value_type, value_type_b>::value, "A and B do not have the same value type.");

    static_assert(std::is_same<memory_space, memory_space_b>::value, "A and B do not have the same memory space.");

    int r_val(0);
    if (std::is_same<memory_space, Kokkos::HostSpace>::value)
      r_val = blas_invoke(diagA, alpha, A, B);
#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<memory_space, Kokkos::CudaSpace>::value || std::is_same<memory_space, Kokkos::CudaUVMSpace>::value)
      r_val = cublas_invoke(member, diagA, alpha, A, B);
#endif
#if defined(KOKKOS_ENABLE_HIP)
    if (std::is_same<memory_space, Kokkos::HIPSpace>::value)
      r_val = rocblas_invoke(member, diagA, alpha, A, B);
#endif
    return r_val;
  }
};

// Trsm with deficient-diagonals (used by LDL no-pivot with zero-diagonal check)
template <typename ArgSide, typename ArgUplo, typename ArgTransA>
struct Trsm_defs<ArgSide, ArgUplo, ArgTransA, Algo::OnDevice> {
  template <typename MemberType, typename DiagType, typename ScalarType, typename ViewTypeA, typename ViewTypeB>
  inline static int invoke(MemberType &member, const DiagType diagA, const ScalarType alpha, const ViewTypeA &A,
                           const ViewTypeB &B) {
    using exec_space = MemberType;
    using policy_type = Kokkos::RangePolicy<exec_space>;
    using value_type = typename ViewTypeA::non_const_value_type;
    const auto &exec_instance = member;
    const value_type zero (0);
    const value_type one (1);

    const ordinal_type m = B.extent(0);
    const ordinal_type n = B.extent(1);
    // Side::Left, Uplo::Upper, Trans::Transpose
    if (ArgSide::param != 'L' || ArgUplo::param != 'U' || ArgTransA::param != 'T')
      printf( " Trsm_defs(%c,%c,%c) not implemented\n",ArgSide::param,ArgUplo::param,ArgTransA::param );
    const auto policy_scale = policy_type(exec_instance, 0, m);
    for (ordinal_type i = 0; i < m; i++) {
      Kokkos::parallel_for(policy_scale, KOKKOS_LAMBDA(const ordinal_type &j) {
        if (A(i, i) == zero ) {
          // if tiny pivot, zero out off-diagonal
          B(i, j) = zero;
        } else {
          if (diagA.param != 'U') {
            // scale
            B(i, j) /= A(i, i);
          }
          // update
          for (ordinal_type k=j+1; k<m; k++) B(k, j) -= A(k, i) * B(i, j);
        }
      });
      // reset zero-pivot with one
      // TODO: move it out, reset after TRSM with off-diagonal blocks
      //if (A(i, i) == zero ) A(i, i) = one;
    }
    return 0;
  }
};

} // namespace Tacho
#endif
