// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_Kernels.hpp
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

  namespace Kernels {
    
    struct Serial {
      template<typename ScalarType,
               typename AViewType,
               typename BViewType,
               typename CViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      gemm_trans_notrans(const ScalarType alpha,
                         const AViewType &A,
                         const BViewType &B,
                         const ScalarType beta,
                         const CViewType &C) {
        //C = beta*C + alpha * A'*B
        const ordinal_type 
          m = C.extent(0),
          n = C.extent(1),
          k = B.extent(0);
        
        for (ordinal_type i=0;i<m;++i)
          for (ordinal_type j=0;j<n;++j) {
            C(i,j) *= beta;
            for (ordinal_type l=0;l<k;++l)
              C(i,j) += alpha*A(l,i)*B(l,j);
          }
      }

      template<typename ScalarType,
               typename AViewType,
               typename BViewType,
               typename CViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      gemm_notrans_trans(const ScalarType alpha,
                         const AViewType &A,
                         const BViewType &B,
                         const ScalarType beta,
                         const CViewType &C) {
        //C = beta*C + alpha * A*B'
        const ordinal_type 
          m = C.extent(0),
          n = C.extent(1),
          k = A.extent(1);
        
        for (ordinal_type i=0;i<m;++i)
          for (ordinal_type j=0;j<n;++j) {
            C(i,j) *= beta;
            for (ordinal_type l=0;l<k;++l)
              C(i,j) += alpha*A(i,l)*B(j,l);
          }
      }
      
      template<typename ScalarType,
               typename AViewType,
               typename xViewType,
               typename yViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      gemv_trans(const ScalarType alpha,
                 const AViewType &A,
                 const xViewType &x,
                 const ScalarType beta,
                 const yViewType &y) {
        //y = beta*y + alpha * A'*x
        const ordinal_type 
          m = y.extent(0),
          n = x.extent(0);
        
        for (ordinal_type i=0;i<m;++i) {
          y(i) *= beta;
          for (ordinal_type j=0;j<n;++j) 
            y(i) += alpha*A(j,i)*x(j);
        }
      }

      template<typename ScalarType,
               typename AViewType,
               typename xViewType,
               typename yViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      gemv_notrans(const ScalarType alpha,
                   const AViewType &A,
                   const xViewType &x,
                   const ScalarType beta,
                   const yViewType &y) {
        //y = beta*y + alpha * A*x
        const ordinal_type 
          m = y.extent(0),
          n = x.extent(0);
        
        for (ordinal_type i=0;i<m;++i) {
          y(i) *= beta;
          for (ordinal_type j=0;j<n;++j) 
            y(i) += alpha*A(i,j)*x(j);
        }
      }

      template<typename matViewType>
      KOKKOS_INLINE_FUNCTION
      static typename matViewType::non_const_value_type
      determinant(const matViewType &mat) {
        INTREPID2_TEST_FOR_ABORT(mat.extent(0) != mat.extent(1), "mat should be a square matrix.");
        INTREPID2_TEST_FOR_ABORT(mat.extent(0) > 3, "Higher dimensions (> 3) are not supported.");

        typename matViewType::non_const_value_type r_val(0);
        const int m = mat.extent(0);
        switch (m) {
        case 1: 
          r_val =  mat(0,0);
          break;
        case 2:
          r_val =  ( mat(0,0) * mat(1,1) -
                     mat(0,1) * mat(1,0) );
          break;
        case 3:
          r_val =  ( mat(0,0) * mat(1,1) * mat(2,2) +
                     mat(1,0) * mat(2,1) * mat(0,2) +
                     mat(2,0) * mat(0,1) * mat(1,2) -
                     mat(2,0) * mat(1,1) * mat(0,2) -
                     mat(0,0) * mat(2,1) * mat(1,2) -
                     mat(1,0) * mat(0,1) * mat(2,2) );
          break;
        }
        return r_val;
      }
      
      template<typename matViewType, 
               typename invViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      inverse(const invViewType &inv,
              const matViewType &mat) {
        INTREPID2_TEST_FOR_ABORT(mat.extent(0) != mat.extent(1), "mat should be a square matrix.");
        INTREPID2_TEST_FOR_ABORT(inv.extent(0) != inv.extent(1), "inv should be a square matrix.");
        INTREPID2_TEST_FOR_ABORT(mat.extent(0) != inv.extent(0), "mat and inv must have the same dimension.");
        INTREPID2_TEST_FOR_ABORT(mat.extent(0) > 3, "Higher dimensions (> 3) are not supported.");
        INTREPID2_TEST_FOR_ABORT(mat.data() == inv.data(), "mat and inv must have different data pointer.");
        
        const auto val = determinant(mat);
        const int m = mat.extent(0);
        switch (m) {
        case 1: {
          inv(0,0) = 1.0/mat(0,0);
          break;
        }
        case 2: {
          inv(0,0) =   mat(1,1)/val;
          inv(1,1) =   mat(0,0)/val;

          inv(1,0) = - mat(1,0)/val;
          inv(0,1) = - mat(0,1)/val;
          break;
        }
        case 3: {
          {
            const auto val0 =   mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2);
            const auto val1 = - mat(1,0)*mat(2,2) + mat(2,0)*mat(1,2);
            const auto val2 =   mat(1,0)*mat(2,1) - mat(2,0)*mat(1,1);
            
            inv(0,0) = val0/val;
            inv(1,0) = val1/val;
            inv(2,0) = val2/val;
          }
          {
            const auto val0 =   mat(2,1)*mat(0,2) - mat(0,1)*mat(2,2);
            const auto val1 =   mat(0,0)*mat(2,2) - mat(2,0)*mat(0,2);
            const auto val2 = - mat(0,0)*mat(2,1) + mat(2,0)*mat(0,1);
            
            inv(0,1) = val0/val;
            inv(1,1) = val1/val;
            inv(2,1) = val2/val;
          }
          {
            const auto val0 =   mat(0,1)*mat(1,2) - mat(1,1)*mat(0,2);
            const auto val1 = - mat(0,0)*mat(1,2) + mat(1,0)*mat(0,2);
            const auto val2 =   mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);

            inv(0,2) = val0/val;
            inv(1,2) = val1/val;
            inv(2,2) = val2/val;
          }
          break;
        }
        }
      }

      template<typename ScalarType,
               typename xViewType,
               typename yViewType,
               typename zViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      z_is_axby(const zViewType &z,
                const ScalarType alpha,
                const xViewType &x,
                const ScalarType beta,
                const yViewType &y) {
        //y = beta*y + alpha*x
        const ordinal_type 
          m = z.extent(0);
        
        for (ordinal_type i=0;i<m;++i) 
          z(i) = alpha*x(i) + beta*y(i);
      }

      template<typename AViewType>
      KOKKOS_INLINE_FUNCTION
      static double
      norm(const AViewType &A, const ENorm normType) {
        typedef typename AViewType::non_const_value_type value_type;
        const ordinal_type m = A.extent(0), n = A.extent(1);
        double r_val = 0;
        switch(normType) {
        case NORM_TWO:{
          for (ordinal_type i=0;i<m;++i)
            for (ordinal_type j=0;j<n;++j)
              r_val += A.access(i,j)*A.access(i,j);
          r_val = sqrt(r_val);
          break;
        }
        case NORM_INF:{
          for (ordinal_type i=0;i<m;++i)
            for (ordinal_type j=0;j<n;++j) {
              const value_type current = Util<value_type>::abs(A.access(i,j));
              r_val = (r_val < current ? current : r_val);
            }
          break;
        }
        case NORM_ONE:{
          for (ordinal_type i=0;i<m;++i)
            for (ordinal_type j=0;j<n;++j)
              r_val += Util<value_type>::abs(A.access(i,j));
          break;
        }
        default: {
          Kokkos::abort("norm type is not supported");
          break;
        }
        }
        return r_val;
      }

      template<typename dstViewType,
               typename srcViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      copy(const dstViewType &dst, const srcViewType &src) { 
        if (dst.data() != src.data()) {
          const ordinal_type m = dst.extent(0), n = dst.extent(1);
          for (ordinal_type i=0;i<m;++i) 
            for (ordinal_type j=0;j<n;++j) 
              dst.access(i,j) = src.access(i,j);
        }
      }

      // y = Ax
      template<typename yViewType,
               typename AViewType,
               typename xViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      static void
      matvec_trans_product_d2( const yViewType &y,
                               const AViewType &A,
                               const xViewType &x ) {
        y(0) = A(0,0)*x(0) + A(1,0)*x(1);
        y(1) = A(0,1)*x(0) + A(1,1)*x(1);
      }

      template<typename yViewType,
               typename AViewType,
               typename xViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      static void
      matvec_trans_product_d3( const yViewType &y,
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
      static void
      matvec_product_d2( const yViewType &y,
                         const AViewType &A,
                         const xViewType &x ) {
        y(0) = A(0,0)*x(0) + A(0,1)*x(1);
        y(1) = A(1,0)*x(0) + A(1,1)*x(1);
      }

      template<typename yViewType,
               typename AViewType,
               typename xViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      static void
      matvec_product_d3( const yViewType &y,
                         const AViewType &A,
                         const xViewType &x ) {
        y(0) = A(0,0)*x(0) + A(0,1)*x(1) + A(0,2)*x(2);
        y(1) = A(1,0)*x(0) + A(1,1)*x(1) + A(1,2)*x(2);
        y(2) = A(2,0)*x(0) + A(2,1)*x(1) + A(2,2)*x(2);
      }

      template<typename yViewType,
               typename AViewType,
               typename xViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      static void
      matvec_product( const yViewType &y,
                      const AViewType &A,
                      const xViewType &x ) {
        switch (y.extent(0)) {
        case 2: matvec_product_d2(y, A, x); break;
        case 3: matvec_product_d3(y, A, x); break;
        default: {
          INTREPID2_TEST_FOR_ABORT(true, "matvec only support dimension 2 and 3 (consider to use gemv interface).");
          break;
        }
        }
      }

      template<typename zViewType,
               typename xViewType,
               typename yViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      static void
      vector_product_d2( const zViewType &z,
                         const xViewType &x,
                         const yViewType &y ) {
        z(0) = x(0)*y(1) - x(1)*y(0);
      }
    
      template<typename zViewType,
               typename xViewType,
               typename yViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      static void
      vector_product_d3( const zViewType &z,
                         const xViewType &x,
                         const yViewType &y ) {
        z(0) = x(1)*y(2) - x(2)*y(1);
        z(1) = x(2)*y(0) - x(0)*y(2);
        z(2) = x(0)*y(1) - x(1)*y(0);      
      }


    };
      


    template<typename xViewType,
             typename yViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    static typename xViewType::value_type
    dot( const xViewType x,
         const yViewType y ) {
      typename xViewType::value_type r_val(0);
      ordinal_type i = 0, iend = x.extent(0);
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
    static typename xViewType::value_type
    dot_d2( const xViewType x,
            const yViewType y ) {
      return ( x(0)*y(0) + x(1)*y(1) );
    }
    
    template<typename xViewType,
             typename yViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    static typename xViewType::value_type
    dot_d3( const xViewType x,
            const yViewType y ) {
      return ( x(0)*y(0) + x(1)*y(1) + x(2)*y(2) );
    }

    template<typename AViewType,
             typename alphaScalarType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void
    scale_mat(       AViewType &A,
               const alphaScalarType alpha ) {
      const ordinal_type
        iend = A.extent(0),
        jend = A.extent(1);

      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          A(i,j) *= alpha;
    }

    template<typename AViewType,
             typename alphaScalarType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void
    inv_scale_mat(       AViewType &A,
                   const alphaScalarType alpha ) {
      const ordinal_type
        iend = A.extent(0),
        jend = A.extent(1);
      
      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          A(i,j) /= alpha;
    }

    template<typename AViewType,
             typename alphaScalarType,
             typename BViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void
    scalar_mult_mat(       AViewType &A,
                     const alphaScalarType alpha,
                     const BViewType &B ) {
      const ordinal_type
        iend = A.extent(0),
        jend = A.extent(1);

      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          A(i,j) = alpha*B(i,j);
    }

    template<typename AViewType,
             typename alphaScalarType,
             typename BViewType>
    KOKKOS_FORCEINLINE_FUNCTION
    static void
    inv_scalar_mult_mat(       AViewType &A,
                         const alphaScalarType alpha,
                         const BViewType &B ) {
      const ordinal_type
        iend = A.extent(0),
        jend = A.extent(1);
      
      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          A(i,j) = B(i,j)/alpha;
    }

  }
}

#endif
