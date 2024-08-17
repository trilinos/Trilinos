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
#ifndef __TACHO_LDL_ON_DEVICE_HPP__
#define __TACHO_LDL_ON_DEVICE_HPP__

/// \file  Tacho_LDL_OnDevice.hpp
/// \brief LDL device solver
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_LDL_External.hpp"

namespace Tacho {

template <> struct LDL<Uplo::Lower, Algo::OnDevice> {
  template <typename ViewTypeA, typename ViewTypeP, typename ViewTypeW>
  inline static int lapack_invoke(const ViewTypeA &A, const ViewTypeP &P, const ViewTypeW &W) {
    return LDL<Uplo::Lower, Algo::External>::invoke(A, P, W);
  }

#if defined(KOKKOS_ENABLE_CUDA)
  template <typename ViewTypeA, typename ViewTypeP, typename ViewTypeW>
  inline static int cusolver_invoke(cusolverDnHandle_t &handle, const ViewTypeA &A, const ViewTypeP &P,
                                    const ViewTypeW &W) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeW::non_const_value_type work_value_type;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    if (m > 0) {
      int *devInfo = (int *)W.data();
      work_value_type *workspace = W.data() + 1;
      int lwork = (W.span() - 1);
      r_val = Lapack<value_type>::sytrf(handle, CUBLAS_FILL_MODE_LOWER, m, A.data(), A.stride_1(), P.data(), workspace,
                                        lwork, devInfo);
    }
    return r_val;
  }

  template <typename ViewTypeA>
  inline static int cusolver_buffer_size(cusolverDnHandle_t &handle, const ViewTypeA &A, int *lwork) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    if (m > 0)
      r_val = Lapack<value_type>::sytrf_buffersize(handle, m, A.data(), A.stride_1(), lwork);
    return r_val;
  }
#endif

#if defined(KOKKOS_ENABLE_HIP)
  template <typename ViewTypeA, typename ViewTypeP, typename ViewTypeW>
  inline static int rocsolver_invoke(rocblas_handle &handle, const ViewTypeA &A, const ViewTypeP &P,
                                     const ViewTypeW &W) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    if (m > 0) {
      int *devInfo = (int *)W.data();
      r_val = Lapack<value_type>::sytrf(handle, rocblas_fill_lower, m, A.data(), A.stride_1(), P.data(), devInfo);
    }
    return r_val;
  }
#endif

  template <typename MemberType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeW>
  inline static int invoke(MemberType &member, const ViewTypeA &A, const ViewTypeP &P, const ViewTypeW &W) {
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
      if (W.span() == 0) {
        int lwork = A.extent(0) * 32;
        return lwork;
      } else {
        r_val = lapack_invoke(A, P, W);
      }
    }

#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<memory_space, Kokkos::CudaSpace>::value ||
        std::is_same<memory_space, Kokkos::CudaUVMSpace>::value) {
      if (W.span() == 0) {
        int lwork;
        r_val = cusolver_buffer_size(member, A, &lwork);
        r_val = lwork + 1;
      } else
        r_val = cusolver_invoke(member, A, P, W);
    }
#endif

#if defined(KOKKOS_ENABLE_HIP)
    if (std::is_same<memory_space, Kokkos::HIPSpace>::value) {
      if (W.span() == 0) {
        r_val = 2;
      } else
        r_val = rocsolver_invoke(member, A, P, W);
    }
#endif
    return r_val;
  }

  template <typename ViewTypeA, typename ViewTypeP, typename ViewTypeD>
  inline static int lapack_modify(const ViewTypeA &A, const ViewTypeP &P, const ViewTypeD &D) {
    return LDL<Uplo::Lower, Algo::External>::modify(A, P, D);
  }

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  template <typename ExecSpaceType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeD>
  inline static int device_modify(ExecSpaceType &exec_instance, const ViewTypeA &A, const ViewTypeP &P,
                                  const ViewTypeD &D) {
    using exec_space = ExecSpaceType;
    typedef typename ViewTypeA::non_const_value_type value_type;
    const ordinal_type m = A.extent(0);

    int r_val(0);
    if (m > 0) {
      value_type *KOKKOS_RESTRICT Aptr = A.data();
      ordinal_type *KOKKOS_RESTRICT ipiv = P.data();
      ordinal_type *KOKKOS_RESTRICT fpiv = ipiv + m;
      ordinal_type *KOKKOS_RESTRICT perm = fpiv + m;
      ordinal_type *KOKKOS_RESTRICT peri = perm + m;

      const value_type one(1), zero(0);
      Kokkos::RangePolicy<exec_space> range_policy(exec_instance, 0, m);
      Kokkos::parallel_for(
          "PermutationSet", range_policy, KOKKOS_LAMBDA(const ordinal_type i) { perm[i] = i; });
      exec_instance.fence();
      Kokkos::parallel_for(
          "ExtractDiagonalsAndPostProcessing", range_policy, KOKKOS_LAMBDA(const ordinal_type j) {
            const bool single = (j == (m - 1));
            for (ordinal_type i = 0; i < m; ++i) {
              if (ipiv[i] < 0) {
                {
                  // first pivot
                  if (single) {
                    ipiv[i] = 0; /// invalidate this pivot
                    fpiv[i] = 0;

                    D(i, 0) = A(i, i);
                    D(i, 1) = A(i + 1, i); /// symmetric
                    A(i, i) = one;
                  }
                }
                {
                  // second pivot
                  i ++;
                  const ordinal_type fla_pivot = -ipiv[i] - i - 1;
                  if (single) {
                    fpiv[i] = fla_pivot;
                  }
                  if (fla_pivot) {
                    value_type *KOKKOS_RESTRICT src = Aptr + i;
                    value_type *KOKKOS_RESTRICT tgt = src + fla_pivot;
                    if (j < (i - 1)) {
                      const ordinal_type idx = j * m;
                      swap(src[idx], tgt[idx]);
                    }
                  }

                  if (single) {
                    D(i, 0) = A(i, i - 1);
                    D(i, 1) = A(i, i);
                    A(i, i - 1) = zero;
                    A(i, i) = one;
                  }
                }
              } else {
                const ordinal_type fla_pivot = ipiv[i] - i - 1;
                if (single) {
                  fpiv[i] = fla_pivot;
                }
                if (fla_pivot) {
                  value_type *src = Aptr + i;
                  value_type *tgt = src + fla_pivot;
                  if (j < i) {
                    const ordinal_type idx = j * m;
                    swap(src[idx], tgt[idx]);
                  }
                }

                if (single) {
                  D(i, 0) = A(i, i);
                  A(i, i) = one;
                }
              }

              /// apply pivots to perm vector
              if (single) {
                if (fpiv[i]) {
                  const ordinal_type pidx = i + fpiv[i];
                  swap(perm[i], perm[pidx]);
                }
              }
            }
          });
      exec_instance.fence();
      Kokkos::parallel_for(
          "PermutationInverse", range_policy, KOKKOS_LAMBDA(const ordinal_type i) { peri[perm[i]] = i; });
    }
    return r_val;
  }
#endif
  template <typename MemberType, typename ViewTypeA, typename ViewTypeP, typename ViewTypeD>
  inline static int modify(MemberType &member, const ViewTypeA &A, const ViewTypeP &P, const ViewTypeD &D) {
    typedef typename ViewTypeA::non_const_value_type value_type;
    typedef typename ViewTypeD::non_const_value_type value_type_d;

    typedef typename ViewTypeA::memory_space memory_space;
    typedef typename ViewTypeP::memory_space memory_space_p;
    typedef typename ViewTypeD::memory_space memory_space_d;

    static_assert(ViewTypeA::rank == 2, "A is not rank 2 view.");
    static_assert(ViewTypeP::rank == 1, "P is not rank 1 view.");
    static_assert(ViewTypeD::rank == 2, "D is not rank 2 view.");

    static_assert(std::is_same<value_type, value_type_d>::value, "A and D do not have the same value type.");

    static_assert(std::is_same<memory_space, memory_space_p>::value, "A and P do not have the same memory space.");

    static_assert(std::is_same<memory_space, memory_space_d>::value, "A and D do not have the same memory space.");

    int r_val(0);
    if (std::is_same<memory_space, Kokkos::HostSpace>::value) {
      r_val = lapack_modify(A, P, D);
    }

#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<memory_space, Kokkos::CudaSpace>::value ||
        std::is_same<memory_space, Kokkos::CudaUVMSpace>::value) {
      r_val = device_modify(member, A, P, D);
    }
#endif
#if defined(KOKKOS_ENABLE_HIP)
    if (std::is_same<memory_space, Kokkos::HIPSpace>::value) {
      r_val = device_modify(member, A, P, D);
    }
#endif
    return r_val;
  }
};
} // namespace Tacho
#endif
