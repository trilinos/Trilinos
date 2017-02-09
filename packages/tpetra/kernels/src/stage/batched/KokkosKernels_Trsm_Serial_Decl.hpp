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





#ifndef __KOKKOSKERNELS_TRSM_SERIAL_DECL_HPP__
#define __KOKKOSKERNELS_TRSM_SERIAL_DECL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include <immintrin.h>

namespace KokkosKernels {

  namespace Serial {

    // specialized for different m and n
    // Solve L(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTriSolveLeftLowerUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTriSolveLeftLowerUnitDiag(const int as0, const int as1,
                                     const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int n,
                 /**/  ValueType *__restrict__ B);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };

    // specialized for different m and n
    // Solve L(m x m) X(m x n) = B(m x n)
    template<int bmn>
    struct InnerTriSolveLeftLowerNonUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTriSolveLeftLowerNonUnitDiag(const int as0, const int as1,
                                        const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int n,
                 /**/  ValueType *__restrict__ B);

      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };

    // specialized for different m and n
    // Solve  X(m x n) U(n x n) = B(m x n)
    template<int bmn>
    struct InnerTriSolveRightUpperUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTriSolveRightUpperUnitDiag(const int as0, const int as1,
                                      const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m,
                 /**/  ValueType *__restrict__ B);
      
      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };    

    // specialized for different m and n
    // Solve  X(m x n) U(n x n) = B(m x n)
    template<int bmn>
    struct InnerTriSolveRightUpperNonUnitDiag {
      const int _as0, _as1, _bs0, _bs1;

      KOKKOS_INLINE_FUNCTION
      InnerTriSolveRightUpperNonUnitDiag(const int as0, const int as1,
                                         const int bs0, const int bs1)
        : _as0(as0), _as1(as1),
          _bs0(bs0), _bs1(bs1) {}
      
      // trisolve
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m,
                 /**/  ValueType *__restrict__ B);
      
      // for remainder
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      int invoke(const ValueType *__restrict__ A,
                 const int m, const int n, 
                 /**/  ValueType *__restrict__ B);
    };    

    template<typename ArgSide,
             typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag,
             typename ArgAlgo>
    struct Trsm {

      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType A,
             /**/  BViewType B);
      
    };

    template<typename ArgSide,
             typename ArgUplo,
             typename ArgTrans,
             typename ArgDiag>
    struct Trsm<ArgSide,ArgUplo,ArgTrans,ArgDiag,void> {

      template<typename ScalarType,
               typename AViewType,
               typename BViewType>
      KOKKOS_INLINE_FUNCTION
      static int
      invoke(const ScalarType alpha,
             const AViewType A,
             /**/  BViewType B) {
        return 0;
      }
      
    };

  }

}

#endif
