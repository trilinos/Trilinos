/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels: Linear Algebra and Graph Kernels
//                 Copyright 2016 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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





#ifndef __KOKKOSKERNELS_GEMM_DECL_HPP__
#define __KOKKOSKERNELS_GEMM_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <immintrin.h>

namespace KokkosKernels {

  ///
  /// Serial Gemm 
  ///

  namespace Serial {
    template<typename ArgTransA,
             typename ArgTransB,
             typename ArgAlgo>
    struct Gemm {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType,
               typename CViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType A,
             const BViewType B,
             const ScalarType beta,
             /**/  CViewType C);
    };
  }

  ///
  /// Team Gemm
  ///

  namespace Team {
    template<typename MemberType,
             typename ArgTransA,
             typename ArgTransB,
             typename ArgAlgo>
    struct Gemm {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType,
               typename CViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const MemberType &member,
             const ScalarType alpha,
             const AViewType A,
             const BViewType B,
             const ScalarType beta,
             /**/  CViewType C);
    };
  }
  
  // specialized for different m and n
  // C(mxn) += alpha * A(mxk) B(kxn)
  template<int mb, int nb>
  struct InnerRankUpdate {
    const int _as0, _as1, _bs0, _bs1, _cs0, _cs1;
    
    KOKKOS_INLINE_FUNCTION
    InnerRankUpdate(const int as0, const int as1, 
                    const int bs0, const int bs1,
                    const int cs0, const int cs1)
      : _as0(as0), _as1(as1), 
        _bs0(bs0), _bs1(bs1), 
        _cs0(cs0), _cs1(cs1) {}
    
    // serial rank update
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int serial_invoke(const ScalarType alpha,
                      const ValueType *__restrict__ A,
                      const ValueType *__restrict__ B,
                      const int k,
                      /**/  ValueType *__restrict__ C);
    
    // serial rank update for remainder
    template<typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int serial_invoke(const ScalarType alpha,
                      const ValueType *__restrict__ A,
                      const ValueType *__restrict__ B,
                      const int m, const int n, const int k,
                      /**/  ValueType *__restrict__ C);

    // team rank update
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int team_invoke_var1(const MemberType &member,
                         const ScalarType alpha,
                         const ValueType *__restrict__ A,
                         const ValueType *__restrict__ B,
                         const int k,
                         /**/  ValueType *__restrict__ C);

    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int team_invoke_var2(const MemberType &member,
                         const ScalarType alpha,
                         const ValueType *__restrict__ A,
                         const ValueType *__restrict__ B,
                         const int k,
                         /**/  ValueType *__restrict__ C);

    // team rank update for remainder
    template<typename MemberType,
             typename ScalarType,
             typename ValueType>
    KOKKOS_INLINE_FUNCTION
    int team_invoke(const MemberType &member,
                    const ScalarType alpha,
                    const ValueType *__restrict__ A,
                    const ValueType *__restrict__ B,
                    const int m, const int n, const int k,
                    /**/  ValueType *__restrict__ C);
    
  };

}

#endif
