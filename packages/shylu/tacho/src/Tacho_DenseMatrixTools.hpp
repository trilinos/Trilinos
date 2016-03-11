#pragma once
#ifndef __TACHO_DENSE_MATRIX_TOOLS_HPP__
#define __TACHO_DENSE_MATRIX_TOOLS_HPP__

/// \file Tacho_DenseMatrixTools.hpp
/// \brief Generic utility function for dense matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)  

#include "Tacho_Util.hpp"

namespace Tacho { 

  class DenseMatrixTools {
  public:
    /// Elementwise copy
    /// ------------------------------------------------------------------
    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (o)
    /// - Blocking with fence (o)

    /// \brief elementwise copy 
    template<typename DenseMatrixTypeA, 
             typename DenseMatrixTypeB>
    KOKKOS_INLINE_FUNCTION
    static void
    copy(DenseMatrixTypeA &A,
         const DenseMatrixTypeB &B) {
      static_assert( Kokkos::Impl
                     ::is_same<
                     typename DenseMatrixTypeA::space_type,
                     typename DenseMatrixTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename DenseMatrixTypeA::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;
      
      space_type::execution_space::fence();      

      Kokkos::parallel_for( range_policy(0, B.NumCols()), 
                            [&](const ordinal_type j) 
                            {
#pragma unroll
                              for (auto i=0;i<B.NumRows();++i)
                                A.Value(i,j) = B.Value(i,j);
                            } );

      space_type::execution_space::fence();
    }

    /// \brief elementwise copy of lower/upper triangular of matrix 
    template<typename DenseMatrixTypeA, 
             typename DenseMatrixTypeB>
    KOKKOS_INLINE_FUNCTION
    static void
    copy(DenseMatrixTypeA &A,
         const int uplo,
         const DenseMatrixTypeB &B) {
      static_assert( Kokkos::Impl
                     ::is_same<
                     typename DenseMatrixTypeA::space_type,
                     typename DenseMatrixTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename DenseMatrixTypeA::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;
      
      space_type::execution_space::fence();      

      switch (uplo) {
      case Uplo::Lower: {
        Kokkos::parallel_for( range_policy(0, B.NumCols()),
                              [&](const ordinal_type j)
                              {
#pragma unroll
                                for (ordinal_type i=j;i<B.NumRows();++i)
                                  A.Value(i, j) = B.Value(i, j);
                              } );
        break;
      }
      case Uplo::Upper: {
        Kokkos::parallel_for( range_policy(0, B.Numcols()),
                              [&](const ordinal_type j)
                              {
#pragma unroll
                                for (ordinal_type i=0;i<(j+1);++i)
                                  A.Value(i, j) = B.Value(i, j);
                              } );
        break;
      }
      }

      space_type::execution_space::fence();
    }

//     /// \brief elementwise copy of matrix b
//     template<typename DenseMatrixTypeA, 
//              typename DenseMatrixTypeB,
//              typename OrdinalArrayTypeIp,
//              typename OrdinalArrayTypeJp>
//     KOKKOS_INLINE_FUNCTION
//     static void
//     copy(DenseMatrixTypeA &A,
//          const DenseMatrixTypeB &B,
//          const OrdinalArrayTypeIp &ip,         
//          const OrdinalArrayTypeJp &jp) { 
      
//       static_assert( Kokkos::Impl::is_same<DenseMatrixTypeA::space_type,DenseMatrixTypeB::space_type>::value,
//                      "Space type of input matrices does not match" );
      
//       typedef DenseMatrixTypeA::space_type space_type;
//       typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;
      
//       const int idx = ip.is_null() * 10 + jp.is_null();

//       space_type::execution_space::fence();      

//       switch (idx) {
//       case 11: { // ip is null and jp is null: no permutation
//         Kokkos::parallel_for( range_policy(0, B.NumCols()), 
//                               KOKKOS_LAMBDA(const ordinal_type j) 
//                               {
// #pragma unroll
//                                 for (auto i=0;i<B.NumRows();++i)
//                                   A.Value(i,j) = B.Value(i,j);
//                               } );
//         break;
//       }
//       case 0: { // ip is not null and jp is not null: row/column permutation 
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               KOKKOS_LAMBDA(const ordinal_type j) 
//                               {
// #pragma unroll 
//                                 for (auto i=0;i<B._m;++i)
//                                   A.Value(i, j) = B.Value(ip(i), jp(j));
//                               } );
//         break;
//       }
//       case 10: { // ip is not null and jp is null: row permutation
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               [&](const ordinal_type j) 
//                               {
// #pragma unroll 
//                                 for (auto i=0;i<B._m;++i)
//                                   A.Value(i, j) = B.Value(ip(i), j);
//                               } );
//         break;
//       } 
//       case 1: { // ip is null and jp is not null: column permutation
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               [&](const ordinal_type j) 
//                               {
//                                 const ordinal_type jj = jp(j);
// #pragma unroll 
//                                 for (auto i=0;i<B._m;++i)
//                                   A.Value(i, j) = B.Value(i, jj);
//                               } );
//         break;
//       }
//       }

//       space_type::execution_space::fence();
//     }

//     /// \brief elementwise copy of lower/upper triangular of matrix b
//     template<typename VT,
//              typename OT,
//              typename ST>
//     KOKKOS_INLINE_FUNCTION
//     void
//     copy(const int uplo, 
//          const DenseMatrixBase<VT,OT,ST,space_type> &b) { 

//       typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Dynamic> > range_policy;

//       space_type::execution_space::fence();

//       switch (uplo) {
//       case Uplo::Lower: {
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               [&](const ordinal_type j) 
//                               { 
// #pragma unroll 
//                                 for (ordinal_type i=j;i<B._m;++i) 
//                                   A.Value(i, j) = B.Value(i, j);
//                               } );
//         break;
//       }
//       case Uplo::Upper: {
//         Kokkos::parallel_for( range_policy(0, B._n), 
//                               [&](const ordinal_type j) 
//                               { 
// #pragma unroll 
//                                 for (ordinal_type i=0;i<(j+1);++i) 
//                                   A.Value(i, j) = B.Value(i, j);
//                               } );
//         break;
//       }
//       }

//       space_type::execution_space::fence();
//     }
    /// ------------------------------------------------------------------


  };

}

#endif
