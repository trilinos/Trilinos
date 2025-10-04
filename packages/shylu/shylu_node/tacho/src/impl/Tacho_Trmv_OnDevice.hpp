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
#ifndef __TACHO_TRMV_ON_DEVICE_HPP__
#define __TACHO_TRMV_ON_DEVICE_HPP__

/// \file  Tacho_Trmv_OnDevice.hpp
/// \brief BLAS triangular solve matrix
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Blas_External.hpp"

namespace Tacho {

template <typename ArgUplo, typename ArgTransA> struct Trmv<ArgUplo, ArgTransA, Algo::OnDevice> {

  template <typename DiagType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int blas_invoke(const DiagType diagA, const ViewTypeA &A, const ViewTypeB &B, const ViewTypeC &C) {
    using value_type =  typename ViewTypeA::non_const_value_type;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    static_assert(ArgUplo::param == 'U' || ArgUplo::param == 'u', "A is not upper-triangular.");

    const ordinal_type mC = C.extent(0), nC = C.extent(1);
    const ordinal_type mA = A.extent(0), nA = A.extent(1);

    if (mC > 0 && nC > 0) {
      const value_type one(1);
      const value_type zero(0);
      const bool transA = (ArgTransA::param != 'N' && ArgTransA::param != 'n');
      const ordinal_type mn = (mA < nA ? mA : nA);

      const auto dB = Kokkos::subview(B, range_type(0, mn), Kokkos::ALL());
      const auto dC = Kokkos::subview(C, range_type(0, mn), Kokkos::ALL());

      if (nC == 1) {
        Blas<value_type>::trmv(ArgUplo::param, ArgTransA::param, diagA.param,
                               mA, A.data(), A.stride(1),
                                   C.data(), C.stride(0));
        if (!transA) {
          if (mA > nA) {
            Blas<value_type>::gemv(ArgTransA::param, mA-nA, nA,
                                   one,  &A(nA,0), A.stride(1),
                                         B.data(), B.stride(0), 
                                   zero, &C(nA,0), C.stride(0));
          } else if (nA > mA) {
            Blas<value_type>::gemv(ArgTransA::param, mA, nA-mA,
                                   one, &A(0,mA), A.stride(1),
                                        &B(mA,0), B.stride(0), 
                                   one, C.data(), C.stride(0));
          }
        } else {
          if (mA > nA) {
            Blas<value_type>::gemv(ArgTransA::param, mA-nA, nA,
                                   one, &A(nA,0), A.stride(1),
                                        &B(nA,0), B.stride(0), 
                                   one, C.data(), C.stride(0));
          } else if (nA > mA) {
            Blas<value_type>::gemv(ArgTransA::param, mA, nA-mA,
                                   one,  &A(0,mA), A.stride(1),
                                         B.data(), B.stride(0), 
                                   zero, &C(mA,0), C.stride(0));
          }
        }
      } else {
        for (ordinal_type j = 0; j < nC; j++) {
          Blas<value_type>::trmv(ArgUplo::param, ArgTransA::param, diagA.param,
                                 mA, A.data(), A.stride(1),
                                     &C(0,j),  C.stride(0));
        }
        //Blas<value_type>::trmm(Side::Left::param, ArgUplo::param, ArgTransA::param, diagA.param,
        //                       mA, nC,
        //                       one, A.data(), A.stride(1),
        //                            C.data(), C.stride(1));
        if (ArgTransA::param == 'N' || ArgTransA::param == 'n') {
          if (mA > nA) {
            Blas<value_type>::gemm(ArgTransA::param, 'N',
                                   mA-nA, nC, nA,
                                   one,  &A(nA,0), A.stride(1),
                                         B.data(), B.stride(1),
                                   zero, &C(nA,0), C.stride(1));
          } else if (nA > mA) {
            Blas<value_type>::gemm(ArgTransA::param, 'N',
                                   mA, nC, nA-mA,
                                   one, &A(0,mA), A.stride(1),
                                        &B(mA,0), B.stride(1),
                                   one, C.data(), C.stride(1));
          }
        } else {
          if (mA > nA) {
            Blas<value_type>::gemm(ArgTransA::param, 'N',
                                   mA, nC, mA-nA,
                                   one, &A(nA,0), A.stride(1),
                                        &B(nA,0), B.stride(1), 
                                   one, C.data(), C.stride(1));
          } else if (nA > mA) {
            Blas<value_type>::gemm(ArgTransA::param, 'N',
                                   nA-mA, nC, mA, 
                                   one,  &A(0,mA), A.stride(1),
                                         B.data(), B.stride(1),  
                                   zero, &C(mA,0), C.stride(1));
          }
        }
      }
    }
    return 0;
  }

#if defined(KOKKOS_ENABLE_CUDA)
  template <typename DiagType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int cublas_invoke(cublasHandle_t &handle, const DiagType diagA, const ViewTypeA &A,
                                  const ViewTypeB &B, const ViewTypeC &C) {
    using value_type = typename ViewTypeA::non_const_value_type;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    const value_type one(1);
    const value_type zero(0);

    const ordinal_type mB = B.extent(0), nB = B.extent(1);
    const ordinal_type mA = A.extent(0), nA = A.extent(1);
    const ordinal_type mC = C.extent(0), nC = C.extent(1);
    const bool transA = (ArgTransA::param != 'N' && ArgTransA::param != 'n');

    int r_val(0);
    if (mB > 0 && nB > 0) {
      const ordinal_type mn = (mA < nA ? mA : nA);
      const auto dB = Kokkos::subview(B, range_type(0, mn), Kokkos::ALL());
      const auto dC = Kokkos::subview(C, range_type(0, mn), Kokkos::ALL());
      Kokkos::deep_copy(dC, dB); // TODO: Can we skip this?
      if (nB == 1) {
        r_val = Blas<value_type>::trmv(handle, ArgUplo::cublas_param, ArgTransA::cublas_param, diagA.cublas_param, mA,
                                       A.data(), A.stride(1), C.data(), C.stride(0));
        if (nA > mA) {
          if (!transA) {
            const auto A12 = Kokkos::subview(A, Kokkos::ALL(), range_type(mA, nA));
            const auto B21 = Kokkos::subview(B, range_type(mA, nA), Kokkos::ALL());
            Blas<value_type>::gemv(handle, ArgTransA::cublas_param, mA, nA-mA,
                                   one, A12.data(), A12.stride(1),
                                        B21.data(), B21.stride(0), 
                                   one, C.data(),   C.stride(0));
          } else {
            const auto A12 = Kokkos::subview(A, Kokkos::ALL(), range_type(mA, nA));
            const auto C21 = Kokkos::subview(C, range_type(mA, mC), Kokkos::ALL());
            Blas<value_type>::gemv(handle, ArgTransA::cublas_param, mA, nA-mA,
                                   one,  A12.data(), A12.stride(1),
                                         B.data(),   B.stride(0), 
                                   zero, C21.data(), C21.stride(0));
          }
        } else if (nA != mA) {
          TACHO_TEST_FOR_ABORT(true, "Tall-skinny TRMV not implemented");
        }
      } else {
        for (ordinal_type j = 0; j < nB; j++) {
          const auto Cj = Kokkos::subview(C, Kokkos::ALL(), range_type(j, j + 1));
          r_val = Blas<value_type>::trmv(handle, ArgUplo::cublas_param, ArgTransA::cublas_param, diagA.cublas_param, mA,
                                         A.data(), A.stride(1), Cj.data(), Cj.stride(0));
        }
        if (nA > mA) {
          if (!transA) {
            const auto A12 = Kokkos::subview(A, Kokkos::ALL(), range_type(mA, nA));
            const auto B21 = Kokkos::subview(B, range_type(mA, nA), Kokkos::ALL());
            Blas<value_type>::gemm(handle, ArgTransA::cublas_param, Trans::NoTranspose::cublas_param,
                                   mC, nC, nA-mA,
                                   one, A12.data(), A12.stride(1),
                                        B21.data(), B21.stride(1),
                                   one, C.data(),   C.stride(1));
          } else {
            const auto A12 = Kokkos::subview(A, Kokkos::ALL(), range_type(mA, nA));
            const auto C21 = Kokkos::subview(C, range_type(mA, mC), Kokkos::ALL());
            Blas<value_type>::gemm(handle, ArgTransA::cublas_param, Trans::NoTranspose::cublas_param,
                                   nA-mA, nC, mA,
                                   one,  A12.data(), A12.stride(1),
                                         B.data(),   B.stride(1),
                                   zero, C21.data(), C21.stride(1));
          }
        } else if (nA != mA) {
          TACHO_TEST_FOR_ABORT(true, "Tall-skinny TRMV not implemented");
        }
      }
    }
    return r_val;
  }
#endif
#if defined(KOKKOS_ENABLE_HIP)
  template <typename DiagType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int rocblas_invoke(rocblas_handle &handle, const DiagType diagA, const ViewTypeA &A,
                                   const ViewTypeB &B, const ViewTypeC &C) {
    using value_type = typename ViewTypeA::non_const_value_type;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

    const value_type one(1);
    const value_type zero(0);

    const ordinal_type mB = B.extent(0), nB = B.extent(1);
    const ordinal_type mA = A.extent(0), nA = A.extent(1);
    const ordinal_type mC = C.extent(0), nC = C.extent(1);
    const bool transA = (ArgTransA::param != 'N' && ArgTransA::param != 'n');

    int r_val(0);
    if (mB > 0 && nB > 0) {
      const ordinal_type mn = (mA < nA ? mA : nA);
      const auto dB = Kokkos::subview(B, range_type(0, mn), Kokkos::ALL());
      const auto dC = Kokkos::subview(C, range_type(0, mn), Kokkos::ALL());
      Kokkos::deep_copy(dC, dB); // TODO: Can we skip this?
      if (nB == 1) {
        r_val = Blas<value_type>::trmv(handle, ArgUplo::rocblas_param,
                                       ArgTransA::rocblas_param, diagA.rocblas_param, mA,
                                       A.data(), A.stride(1), C.data(), C.stride(0));
        if (nA > mA) {
          if (!transA) {
            const auto A12 = Kokkos::subview(A, Kokkos::ALL(), range_type(mA, nA));
            const auto B21 = Kokkos::subview(B, range_type(mA, nA), Kokkos::ALL());
            Blas<value_type>::gemv(handle, ArgTransA::rocblas_param, mA, nA-mA,
                                   one, A12.data(), A12.stride(1),
                                        B21.data(), B21.stride(0), 
                                   one, C.data(),   C.stride(0));
          } else {
            const auto A12 = Kokkos::subview(A, Kokkos::ALL(), range_type(mA, nA));
            const auto C21 = Kokkos::subview(C, range_type(mA, mC), Kokkos::ALL());
            Blas<value_type>::gemv(handle, ArgTransA::rocblas_param, mA, nA-mA,
                                   one,  A12.data(), A12.stride(1),
                                         B.data(),   B.stride(0), 
                                   zero, C21.data(), C21.stride(0));
          }
        } else if (nA != mA) {
          TACHO_TEST_FOR_ABORT(true, "Tall-skinny TRMV not implemented");
        }
      } else {
        for (ordinal_type j = 0; j < nB; j++) {
          const auto Cj = Kokkos::subview(C, Kokkos::ALL(), range_type(j, j + 1));
          r_val = Blas<value_type>::trmv(handle, ArgUplo::rocblas_param,
                                         ArgTransA::rocblas_param, diagA.rocblas_param, mA,
                                         A.data(), A.stride(1), Cj.data(), Cj.stride(0));
        }
        if (nA > mA) {
          if (!transA) {
            const auto A12 = Kokkos::subview(A, Kokkos::ALL(), range_type(mA, nA));
            const auto B21 = Kokkos::subview(B, range_type(mA, nA), Kokkos::ALL());
            Blas<value_type>::gemm(handle, ArgTransA::rocblas_param,
                                   Trans::NoTranspose::rocblas_param,
                                   mC, nC, nA-mA,
                                   one, A12.data(), A12.stride(1),
                                        B21.data(), B21.stride(1), 
                                   one, C.data(),   C.stride(1));
          } else {
            const auto A12 = Kokkos::subview(A, Kokkos::ALL(), range_type(mA, nA));
            const auto C21 = Kokkos::subview(C, range_type(mA, mC), Kokkos::ALL());
            Blas<value_type>::gemm(handle, ArgTransA::rocblas_param,
                                   Trans::NoTranspose::rocblas_param,
                                   nA-mA, nC, mA,
                                   one,  A12.data(), A12.stride(1),
                                         B.data(),   B.stride(1), 
                                   zero, C21.data(), C21.stride(1));
          }
        } else if (nA != mA) {
          TACHO_TEST_FOR_ABORT(true, "Tall-skinny TRMV not implemented");
        }
      }
    }
    return r_val;
  }
#endif

  template <typename MemberType, typename DiagType, typename ViewTypeA, typename ViewTypeB, typename ViewTypeC>
  inline static int invoke(MemberType &member, const DiagType diagA, const ViewTypeA &A, const ViewTypeB &B, const ViewTypeC &C) {
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
      r_val = blas_invoke(diagA, A, B, C);
#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<memory_space, Kokkos::CudaSpace>::value || std::is_same<memory_space, Kokkos::CudaUVMSpace>::value)
      r_val = cublas_invoke(member, diagA, A, B, C);
#endif
#if defined(KOKKOS_ENABLE_HIP)
    if (std::is_same<memory_space, Kokkos::HIPSpace>::value)
      r_val = rocblas_invoke(member, diagA, A, B, C);
#endif

    return r_val;
  }
};

} // namespace Tacho
#endif
