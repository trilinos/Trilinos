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





#ifndef __KOKKOSKERNELS_GEMM_SERIAL_IMPL_HPP__
#define __KOKKOSKERNELS_GEMM_SERIAL_IMPL_HPP__


/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "KokkosKernels_Util.hpp"

namespace KokkosKernels {

  namespace Serial {

    template<>
    template<typename ScalarType, 
             typename AViewType, 
             typename BViewType, 
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int 
    Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Triple>::
    invoke(const ScalarType alpha,
           const AViewType A,
           const BViewType B,
           const ScalarType beta,
           /**/  CViewType C,
           /**/  void *aux) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      
      typedef typename CViewType::value_type value_type;

      if      (beta == 0) Util::set  (C, value_type(0)   );
      else if (beta != 1) Util::scale(C, value_type(beta));
      
      if (alpha != 0) {
        const int 
          m = C.dimension(0), 
          n = C.dimension(1), 
          k = A.dimension(1);

        for (int i=0;i<m;++i) {
          value_type 
            *__restrict__ pC = &C(i,0);
          for (int j=0;j<n;++j,pC+=C.stride_1()) {
            value_type 
              *__restrict__ pA = &A(i,0),
              *__restrict__ pB = &B(0,j);
            for (int p=0;p<k;++p,pA+=A.stride_1(),pB+=B.stride_0()) 
              *pC += alpha * *pA * *pB;
          }
        }
      }
      return 0;
    }


    template<>
    template<typename ScalarType, 
             typename AViewType, 
             typename BViewType, 
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int 
    Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::Blocked>::
    invoke(const ScalarType alpha,
           const AViewType A,
           const BViewType B,
           const ScalarType beta,
           /**/  CViewType C,
           /**/  void *aux) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      
      constexpr int mc = 32, kc = 32, nc = 32;

      typedef typename CViewType::value_type value_type;

      if      (beta == 0) Util::set  (C, value_type(0)   );
      else if (beta != 1) Util::scale(C, value_type(beta));
      
      if (alpha != 0) {
        const int 
          m = C.dimension(0), 
          n = C.dimension(1), 
          k = A.dimension(1);

        for (int j=0;j<n;j+=nc) {
          const int tj = n-j, jb = (tj < nc ? tj : nc);
          for (int p=0;p<k;p+=kc) {
            const int tp = k-p, pb = (tp < kc ? tp : kc);
            for (int i=0;i<m;i+=mc) {
              const int ti = m-i, ib = (ti < mc ? ti : mc);
              {
                value_type 
                  *__restrict__ pB = &B(p,j);
                for (int jj=0;jj<jb;++jj,pB+=B.stride_1()) {
                  value_type 
                    *__restrict__ pA = &A(i,p),
                    *__restrict__ sC = &C(i,j+jj);
                  for (int ii=0;ii<ib;++ii,sC+=C.stride_0(),pA+=A.stride_0()) {
                    value_type
                      *__restrict__ sA = pA,
                      *__restrict__ sB = pB;
                    for (int pp=0;pp<pb;++pp,sA+=A.stride_1(),sB+=B.stride_0())
                      *sC += alpha * *sA * *sB;
                  }
                }
              }

            }
          }
        }
      }
      return 0;
    }

    template<>
    template<typename ScalarType, 
             typename AViewType, 
             typename BViewType, 
             typename CViewType>
    KOKKOS_INLINE_FUNCTION
    int 
    Gemm<Trans::NoTranspose,Trans::NoTranspose,Algo::Gemm::BlockedWithPacking>::
    invoke(const ScalarType alpha,
           const AViewType A,
           const BViewType B,
           const ScalarType beta,
           /**/  CViewType C,
           /**/  void *aux) {
      // C = beta C + alpha A B
      // C (m x n), A(m x k), B(k x n)
      
      // more precisely, 
      // VL*KC*KC ~ L1/2
      // VL*MC*KC ~ L2/2
      // VL*NC*KC ~ L3/2
      constexpr int mc = 32, kc = 32, nc = 32;

      typedef typename CViewType::value_type value_type;

      if      (beta == 0) Util::set  (C, value_type(0)   );
      else if (beta != 1) Util::scale(C, value_type(beta));
      
      if (alpha != 0) {
        const int 
          m = C.dimension(0), 
          n = C.dimension(1), 
          k = A.dimension(1);

        // wB mapped to L3 cache, wA mapped to L2 cache, a strip of wB mapped to L1
        value_type *__restrict__ wA = (value_type*)aux; aux += (mc*kc*sizeof(value_type));
        value_type *__restrict__ wB = (value_type*)aux;

        for (int j=0;j<n;j+=nc) {
          const int tj = n-j, jb = (tj < nc ? tj : nc);
          for (int p=0;p<k;p+=kc) {
            const int tp = k-p, pb = (tp < kc ? tp : kc);
            // pack B columwise
            {
              value_type 
                *__restrict__ pTo = &wB[0];
              for (int jj=0;jj<jb;++jj) {
                value_type
                  *__restrict__ pFrom = &B(p,j);
                for (int ii=0;ii<pb;++ii,pFrom+=B.stride_0()) 
                  *pTo++ = *pFrom;
              }
            }
            for (int i=0;i<m;i+=mc) {
              const int ti = m-i, ib = (ti < mc ? ti : mc);
              // pack A rowwise
              {
                value_type 
                  *__restrict__ pTo = &wA[0]; 
                for (int ii=0;ii<ib;++ii) {
                  value_type
                    *__restrict__ pFrom = &A(i+ii,p);
                  for (int jj=0;jj<pb;++jj,pFrom+=A.stride_1()) 
                    *pTo++ = *pFrom; 
                }
              }
              {
                value_type 
                  *__restrict__ pB = &wB[0];
                for (int jj=0;jj<jb;++jj,pB+=pb) {
                  value_type 
                    *__restrict__ pA = &wA[0],
                    *__restrict__ sC = &C(i,j+jj);
                  for (int ii=0;ii<ib;++ii,sC+=C.stride_0()) {
                    value_type
                      *__restrict__ sB = pB;
                    for (int pp=0;pp<pb;++pp)
                      *sC += alpha * *pA++ * *sB++;
                  }
                }
              }

            }
          }
        }
      }
      return 0;
    }



  }

}

#endif
