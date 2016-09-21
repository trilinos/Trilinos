// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_Kernels.hpp
    \brief  Header file for small functions used in Intrepid2.
    \author Created by Kyungjoo Kim
*/

#ifndef __INTREPID2_KERNELS_HPP__
#define __INTREPID2_KERNELS_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Kokkos_Core.hpp"

namespace Intrepid2 {

  template<int N, typename ViewType>
  class ViewAdapter {
  private:
    const ordinal_type (&_i)[N];
    const ViewType &_v;

  public:
    KOKKOS_INLINE_FUNCTION
    ViewAdapter(const ordinal_type (&idx_)[N],
                const ViewType &v_)
      : _i(idx_), _v(v_) {}

    KOKKOS_FORCEINLINE_FUNCTION
    typename ViewType::reference_type
    operator()() const {
      switch (N) {
      case 1: return _v(_i[0]                     );
      case 2: return _v(_i[0], _i[1]              );
      case 3: return _v(_i[0], _i[1], _i[2]       );
      case 4: return _v(_i[0], _i[1], _i[2], _i[3]);
      default:
        INTREPID2_TEST_FOR_ABORT( true, ">>> ERROR (Vi): N > 4 is not supported" );
      }
    }

    KOKKOS_FORCEINLINE_FUNCTION
    typename ViewType::reference_type
    operator()(const ordinal_type i0) const {
      switch (N) {
      case 1: return _v(_i[0],                      i0);
      case 2: return _v(_i[0], _i[1],               i0);
      case 3: return _v(_i[0], _i[1], _i[2],        i0);
      case 4: return _v(_i[0], _i[1], _i[2], _i[3], i0);
      default:
        INTREPID2_TEST_FOR_ABORT( true, ">>> ERROR (Vi): N > 4 is not supported" );
      }
    }

    KOKKOS_FORCEINLINE_FUNCTION
    typename ViewType::reference_type
    operator()(const ordinal_type i0, 
               const ordinal_type i1) const {
      switch (N) {
      case 1: return _v(_i[0],                      i0, i1);
      case 2: return _v(_i[0], _i[1],               i0, i1);
      case 3: return _v(_i[0], _i[1], _i[2],        i0, i1);
      case 4: return _v(_i[0], _i[1], _i[2], _i[3], i0, i1);
      default:
        INTREPID2_TEST_FOR_ABORT( true, ">>> ERROR (Vi): N > 4 is not supported" );
      }
    }
    
    KOKKOS_FORCEINLINE_FUNCTION
    ordinal_type
    dimension(ordinal_type i) {
      return _v.dimension(i + N);
    }
  };
  
  namespace Kernels {
    
    // y = Ax
    template<typename yViewType,
             typename AViewType,
             typename xViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    void
    matvec_trans_product_d2( /**/  yViewType &y,
                             const AViewType &A,
                             const xViewType &x ) {
      y(0) = A(0,0)*x(0) + A(1,0)*x(1);
      y(1) = A(0,1)*x(0) + A(1,1)*x(1);
    }

    template<typename yViewType,
             typename AViewType,
             typename xViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    void
    matvec_trans_product_d3( /**/  yViewType &y,
                             const AViewType &A,
                             const xViewType &x ) {
      y(0) = A(0,0)*x(0) + A(1,0)*x(1) + A(2,0)*x(2);
      y(1) = A(0,1)*x(0) + A(1,1)*x(1) + A(2,1)*x(2);
      y(2) = A(0,2)*x(0) + A(1,2)*x(1) + A(2,2)*x(2);
    }

    // y = Ax
    template<typename yViewType,
             typename AViewType,
             typename xViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    void
    matvec_product_d2( /**/  yViewType &y,
                       const AViewType &A,
                       const xViewType &x ) {
      y(0) = A(0,0)*x(0) + A(0,1)*x(1);
      y(1) = A(1,0)*x(0) + A(1,1)*x(1);
    }

    template<typename yViewType,
             typename AViewType,
             typename xViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    void
    matvec_product_d3( /**/  yViewType &y,
                       const AViewType &A,
                       const xViewType &x ) {
      y(0) = A(0,0)*x(0) + A(0,1)*x(1) + A(0,2)*x(2);
      y(1) = A(1,0)*x(0) + A(1,1)*x(1) + A(1,2)*x(2);
      y(2) = A(2,0)*x(0) + A(2,1)*x(1) + A(2,2)*x(2);
    }

    template<typename AViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    typename AViewType::value_type
    det_d1( const AViewType A ) {
      return A(0,0);
    }

    template<typename AViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    typename AViewType::value_type
    det_d2( const AViewType A ) {
      return ( A(0,0) * A(1,1) -
               A(0,1) * A(1,0) );
    }
    
    template<typename AViewType>
    typename AViewType::value_type
    det_d3( const AViewType A ) {
      return ( A(0,0) * A(1,1) * A(2,2) +
               A(1,0) * A(2,1) * A(0,2) +
               A(2,0) * A(0,1) * A(1,2) -
               A(2,0) * A(1,1) * A(0,2) -
               A(0,0) * A(2,1) * A(1,2) -
               A(1,0) * A(0,1) * A(2,2) );
    }

    template<typename xViewType,
             typename yViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    typename xViewType::value_type
    dot( const xViewType x,
         const yViewType y ) {
      typename xViewType::value_type r_val(0);
      ordinal_type i = 0, iend = x.dimension(0);
      for (;i<iend;i+=4) 
        r_val += ( x(i  )*y(i  ) + 
                   x(i+1)*y(i+1) + 
                   x(i+2)*y(i+2) + 
                   x(i+3)*y(i+3) );
      for (;i<iend;++i)
        r_val += x(i)*y(i);
      
      return r_val;
    }

    template<typename xViewType,
             typename yViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    typename xViewType::value_type
    dot_d2( const xViewType x,
            const yViewType y ) {
      return ( x(0)*y(0) + x(1)*y(1) );
    }
    
    template<typename xViewType,
             typename yViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    typename xViewType::value_type
    dot_d3( const xViewType x,
            const yViewType y ) {
      return ( x(0)*y(0) + x(1)*y(1) + x(2)*y(2) );
    }

    template<typename AViewType,
             typename alphaScalarType>
    KOKKOS_FORCEINLINE_FUNCTION
    void
    scale_mat( /**/  AViewType &A,
               const alphaScalarType alpha ) {
      const ordinal_type
        iend = A.dimension(0),
        jend = A.dimension(1);

      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          A(i,j) *= alpha;
    }

    template<typename AViewType,
             typename alphaScalarType>
    KOKKOS_FORCEINLINE_FUNCTION
    void
    inv_scale_mat( /**/  AViewType &A,
                   const alphaScalarType alpha ) {
      const ordinal_type
        iend = A.dimension(0),
        jend = A.dimension(1);
      
      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          A(i,j) /= alpha;
    }

    template<typename AViewType,
             typename alphaScalarType,
             typename BViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    void
    scalar_mult_mat( /**/  AViewType &A,
                     const alphaScalarType alpha,
                     const BViewType &B ) {
      const ordinal_type
        iend = A.dimension(0),
        jend = A.dimension(1);

      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          A(i,j) = alpha*B(i,j);
    }

    template<typename AViewType,
             typename alphaScalarType,
             typename BViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    void
    inv_scalar_mult_mat( /**/  AViewType &A,
                         const alphaScalarType alpha,
                         const BViewType &B ) {
      const ordinal_type
        iend = A.dimension(0),
        jend = A.dimension(1);
      
      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          A(i,j) = B(i,j)/alpha;
    }

  }
}

#endif
