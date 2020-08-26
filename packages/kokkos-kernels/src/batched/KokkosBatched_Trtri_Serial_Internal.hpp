/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef __KOKKOSBATCHED_TRTRI_SERIAL_INTERNAL_HPP__
#define __KOKKOSBATCHED_TRTRI_SERIAL_INTERNAL_HPP__

#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Trmm_Serial_Internal.hpp"

namespace KokkosBatched {

  template<typename AlgoType>
  struct SerialTrtriInternalLower {
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int 
    invoke(const bool use_unit_diag,
           const int am, const int an, 
           ValueType *__restrict__ A, const int as0, const int as1);
  };

  template<typename AlgoType>
  struct SerialTrtriInternalUpper {
    template<typename ValueType>
    KOKKOS_INLINE_FUNCTION
    static int 
    invoke(const bool use_unit_diag,
           const int am, const int an, 
           ValueType *__restrict__ A, const int as0, const int as1);
  };

  template<>
  template<typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  SerialTrtriInternalLower<Algo::Trtri::Unblocked>::
  invoke(const bool use_unit_diag,
         const int am, const int an,
         ValueType *__restrict__ A, const int as0, const int as1) {
    ValueType one(1.0), zero(0.0), A_ii;
    if (!use_unit_diag) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      // Check for singularity
      for (int i = 0; i < am; ++i)
        if (A[i*as0 + i*as1] == zero)
          return i+1;
    }

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = am - 1; i >= 0; --i) {
      A[i*as0 + i*as1] = one / A[i*as0 + i*as1];

      if (i < am - 1) {
        if (use_unit_diag)
          A_ii = -one;
        else
          A_ii = -A[i*as0 + i*as1];

        ValueType *__restrict__ A_subblock = &A[(i+1)*as0 + (i+1)*as1];
        int A_subblock_m = am - i - 1, 
            A_subblock_n = am - i - 1;
        ValueType *__restrict__ A_col_vec  = &A[(i+1)*as0 + i*as1];
        int A_col_vec_m  = am - i - 1,
            A_col_vec_n  = 1;
        // TRMV/TRMM −− x=Ax
        // A((j+1):n,j) = A((j+1):n,(j+1):n) ∗ A((j+1):n,j) ;
        SerialTrmmInternalLeftLower<Algo::Trmm::Unblocked>::invoke(use_unit_diag,
                                                                    false,
                                                                    A_subblock_m, A_subblock_n,
                                                                    A_col_vec_m,  A_col_vec_n,
                                                                    one,
                                                                    A_subblock, as0, as1,
                                                                    A_col_vec,  as0, as1);
        
        // SCAL -- x=ax
        // A((j+1):n,j) = A_ii * A((j+1):n,j)
        SerialScaleInternal::invoke(A_col_vec_m, A_col_vec_n, A_ii, A_col_vec, as0, as1);
      }
    }
    return 0;
  }

  template<>
  template<typename ValueType>
  KOKKOS_INLINE_FUNCTION
  int
  SerialTrtriInternalUpper<Algo::Trtri::Unblocked>::
  invoke(const bool use_unit_diag,
         const int am, const int an,
         ValueType *__restrict__ A, const int as0, const int as1) {
    ValueType one(1.0), zero(0.0), A_ii;


    if (!use_unit_diag) {
#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
      // Check for singularity
      for (int i = 0; i < am; ++i)
        if (A[i*as0 + i*as1] == zero)
          return i+1;
    }

#if defined(KOKKOS_ENABLE_PRAGMA_UNROLL)
#pragma unroll
#endif
    for (int i = 0; i < am; ++i) {
      A[i*as0 + i*as1] = one / A[i*as0 + i*as1];

      if (i > 0) {
        if (use_unit_diag)
          A_ii = -one;
        else
          A_ii = -A[i*as0 + i*as1];

        ValueType *__restrict__ A_subblock = &A[0*as0 + 0*as1];
        int A_subblock_m = i,
            A_subblock_n = i;
        ValueType *__restrict__ A_col_vec  = &A[0*as0 + i*as1];
        int A_col_vec_m  = i,
            A_col_vec_n  = 1;
        // TRMV/TRMM −− x=Ax
        // A(1:(j-1),j) = A(1:(j-1),1:(j-1)) ∗ A(1:(j-1),j) ;
        //SerialTrmm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::NoUnit,Algo::Trmm::Unblocked>
        SerialTrmmInternalLeftUpper<Algo::Trmm::Unblocked>::invoke(use_unit_diag,
                                                                    false,
                                                                    A_subblock_m, A_subblock_n,
                                                                    A_col_vec_m,  A_col_vec_n,
                                                                    one,
                                                                    A_subblock, as0, as1,
                                                                    A_col_vec,  as0, as1);
        
        // SCAL -- x=ax
        // A((j+1):n,j) = A_ii * A((j+1):n,j)
        SerialScaleInternal::invoke(A_col_vec_m, A_col_vec_n, A_ii, A_col_vec, as0, as1);
      }
    }
    return 0;
  }
} // namespace KokkosBatched
#endif // __KOKKOSBATCHED_TRTRI_SERIAL_INTERNAL_HPP__
