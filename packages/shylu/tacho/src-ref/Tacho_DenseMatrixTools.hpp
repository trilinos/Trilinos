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
    /// - Callable in KokkosFunctors (x), no team interface
    /// - Blocking with fence (o)

    /// \brief elementwise copy 
    template<typename DenseMatrixTypeA, 
             typename DenseMatrixTypeB>
    KOKKOS_INLINE_FUNCTION
    static void
    copy(DenseMatrixTypeA &A,
         const DenseMatrixTypeB &B) {
      static_assert( Kokkos::Impl::is_same<
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
                              //#pragma unroll
                              for (auto i=0;i<B.NumRows();++i)
                                A.Value(i,j) = B.Value(i,j);
                            } );

      space_type::execution_space::fence();
    }

    /// \brief elementwise copy of lower/upper triangular of matrix including diagonals
    template<typename DenseMatrixTypeA, 
             typename DenseMatrixTypeB>
    KOKKOS_INLINE_FUNCTION
    static void
    copy(DenseMatrixTypeA &A,
         const int uplo,
         const int offset,
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
                                //#pragma unroll
                                for (ordinal_type i=(j+offset);i<B.NumRows();++i)
                                  A.Value(i, j) = B.Value(i, j);
                              } );
        break;
      }
      case Uplo::Upper: {
        Kokkos::parallel_for( range_policy(0, B.NumCols()),
                              [&](const ordinal_type j)
                              {
                                //#pragma unroll
                                for (ordinal_type i=0;i<=(j-offset);++i)
                                  A.Value(i, j) = B.Value(i, j);
                              } );
        break;
      }
      }

      space_type::execution_space::fence();
    }

    /// Flat to Hier
    /// ------------------------------------------------------------------

    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (o)

    /// \brief compute dimension of hier matrix when the flat is divided by mb x nb
    template<typename DenseMatrixFlatType,
             typename OrdinalType>
    KOKKOS_INLINE_FUNCTION
    static void
    getDimensionOfHierMatrix(OrdinalType &hm,
                             OrdinalType &hn,
                             const DenseMatrixFlatType &flat,
                             const OrdinalType mb,
                             const OrdinalType nb) {
      const auto fm = flat.NumRows();
      const auto fn = flat.NumCols();

      const auto mbb = (mb == 0 ? fm : mb);
      const auto nbb = (nb == 0 ? fn : nb);

      hm = fm/mbb + (fm%mbb > 0);
      hn = fn/nbb + (fn%nbb > 0);
    }
      
    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (o)
    /// - Blocking with fence (o)

    /// \brief fill hier matrix 
    template<typename DenseMatrixHierType,
             typename DenseMatrixFlatType,
             typename OrdinalType>
    KOKKOS_INLINE_FUNCTION
    static void
    getHierMatrix(DenseMatrixHierType &hier,
                  const DenseMatrixFlatType &flat,
                  const OrdinalType mb,
                  const OrdinalType nb) {
      static_assert( Kokkos::Impl::is_same<
                     typename DenseMatrixHierType::space_type,
                     typename DenseMatrixFlatType::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename DenseMatrixHierType::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;

      const OrdinalType hm = hier.NumRows(), hn = hier.NumCols();
      const OrdinalType fm = flat.NumRows(), fn = flat.NumCols();

      space_type::execution_space::fence();

      Kokkos::parallel_for( range_policy(0, hn),
                            [&](const ordinal_type j)
                            {
                              const OrdinalType offn = nb*j;
                              const OrdinalType ntmp = offn + nb; 
                              const OrdinalType n    = ntmp < fn ? nb : (fn - offn); 

                              //#pragma unroll
                              for (ordinal_type i=0;i<hm;++i) {
                                const OrdinalType offm = mb*i;
                                const OrdinalType mtmp = offm + mb; 
                                const OrdinalType m    = mtmp < fm ? mb : (fm - offm); 
                                hier.Value(i, j).setView(flat, offm, m,
                                                         /**/  offn, n);
                              }
                            } );

      space_type::execution_space::fence();
    }

    /// \brief create hier matrix 
    template<typename DenseMatrixHierType,
             typename DenseMatrixFlatType,
             typename OrdinalType>
    KOKKOS_INLINE_FUNCTION
    static void
    createHierMatrix(DenseMatrixHierType &hier,
                     const DenseMatrixFlatType &flat,
                     const OrdinalType mb,
                     const OrdinalType nb) {
      OrdinalType hm, hn;
      getDimensionOfHierMatrix(hm, hn,
                               flat, 
                               mb, nb);
      
      hier.create(hm, hn);
      
      getHierMatrix(hier, flat, mb , nb);
    }
    
  };

}

#endif



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
