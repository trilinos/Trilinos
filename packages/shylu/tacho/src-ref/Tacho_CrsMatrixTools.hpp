#ifndef __TACHO_CRS_MATRIX_TOOLS_HPP__
#define __TACHO_CRS_MATRIX_TOOLS_HPP__

/// \file Tacho_CrsMatrixTools.hpp
/// \brief Generic utility function for crs matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)  

#include "Tacho_Util.hpp"

namespace Tacho { 

  class CrsMatrixTools {
  public:
    /// Elementwise copy
    /// ------------------------------------------------------------------
    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (x), no team interface
    /// - Blocking with fence (o)

    /// \brief elementwise copy 
    template<typename CrsMatBaseTypeA, 
             typename CrsMatBaseTypeB>
    KOKKOS_INLINE_FUNCTION
    static void
    copy(CrsMatBaseTypeA &A,
         const CrsMatBaseTypeB &B) {
      static_assert( Kokkos::Impl::is_same<
                     typename CrsMatBaseTypeA::space_type,
                     typename CrsMatBaseTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename CrsMatBaseTypeA::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;

      space_type::execution_space::fence();      

      // assume that rowpt begin and end arrays are separated.
      Kokkos::parallel_for( range_policy_type(0, B.NumRows()),
                            [&](const ordinal_type i) 
                            {
                              A.RowPtrBegin(i) = B.RowPtrBegin(i);
                              A.RowPtrEnd(i) = B.RowPtrEnd(i);
                            } );

      Kokkos::parallel_for( range_policy_type(0, B.NumNonZeros()), 
                            [&](const ordinal_type k) 
                            {
                              A.Col(k) = B.Col(k);
                              A.Value(k) = B.Value(k);
                            } );
      
      space_type::execution_space::fence();
      A.setNumNonZeros();
    }
    
    /// \brief elementwise copy of lower/upper triangular of matrix
    template<typename CrsMatBaseTypeA, 
             typename CrsMatBaseTypeB>
    KOKKOS_INLINE_FUNCTION
    static void
    copy(CrsMatBaseTypeA &A,
         const int uplo,
         const int offset,
         const CrsMatBaseTypeB &B) {
      static_assert( Kokkos::Impl
                     ::is_same<
                     typename CrsMatBaseTypeA::space_type,
                     typename CrsMatBaseTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      //typedef typename CrsMatBaseTypeA::space_type space_type;
      //typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;
      
      //space_type::execution_space::fence();      

      switch (uplo) {
      case Uplo::Upper: {
        // parallel for  : compute size of each row dimension
        // parallel scan : compute offsets 
        // parallel for  : assignment
        // for now, sequential
        size_type nnz = 0;
        for (ordinal_type i=0;i<B.NumRows();++i) {
          auto cols = B.ColsInRow(i);
          auto vals = B.ValuesInRow(i);
          const int ioffset = i + offset;          
          A.RowPtrBegin(i) = nnz;
          for (auto idx=0;idx<cols.dimension_0();++idx) {
            if (ioffset <= cols[idx]) {
              A.Col(nnz) = cols[idx];
              A.Value(nnz) = vals[idx];
              ++nnz;
            }
          }
          A.RowPtrEnd(i) = nnz;
        }
        A.setNumNonZeros();
        break;
      }
      case Uplo::Lower: {
        // parallel for  : compute size of each row dimension
        // parallel scan : compute offsets 
        // parallel for  : assignment
        // for now, sequential
        size_type nnz = 0;
        for (ordinal_type i=0;i<B.NumRows();++i) {
          auto cols = B.ColsInRow(i);
          auto vals = B.ValuesInRow(i);
          const int ioffset = i - offset;
          A.RowPtrBegin(i) = nnz;
          for (auto idx=0;idx<cols.dimension_0();++idx) {
            if (ioffset >= cols[idx]) {
              A.Col(nnz) = cols[idx];
              A.Value(nnz) = vals[idx];
              ++nnz;
            }
          }
          A.RowPtrEnd(i) = nnz;
        }
        A.setNumNonZeros();
        break;
      }
      }
      //space_type::execution_space::fence();
    }

    /// \brief elementwise copy with permutation
    template<typename CrsMatBaseTypeA, 
             typename CrsMatBaseTypeB,
             typename OrdinalTypeArray>
    KOKKOS_INLINE_FUNCTION
    static void
    copy(CrsMatBaseTypeA &A,
         const OrdinalTypeArray &p,
         const OrdinalTypeArray &ip,
         const CrsMatBaseTypeB &B) {
      static_assert( Kokkos::Impl
                     ::is_same<
                     typename CrsMatBaseTypeA::space_type,
                     typename CrsMatBaseTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      
      typedef typename CrsMatBaseTypeA::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;

      // create work space
      CrsMatBaseTypeA W;
      W.createConfTo(A);

      space_type::execution_space::fence();      

      // column exchange
      if (p.dimension_0()) {
        // structure copy
        Kokkos::parallel_for( range_policy_type(0, B.NumRows()),
                              [&](const ordinal_type i)
                              {
                                W.RowPtrBegin(i) = B.RowPtrBegin(i);
                                W.RowPtrEnd(i) = B.RowPtrEnd(i);
                              } );
        // value copy
        Kokkos::parallel_for( range_policy_type(0, B.NumRows()),
                              [&](const ordinal_type i)
                              {
                                const auto B_cols = B.ColsInRow(i);
                                const auto B_vals = B.ValuesInRow(i);
                                
                                const auto W_cols = W.ColsInRow(i);
                                const auto W_vals = W.ValuesInRow(i);
                                
                                for (size_type j=0;j<B_cols.dimension_0();++j) {
                                  W_cols[j] = p[B_cols[j]];
                                  W_vals[j] = B_vals[j];
                                }
                              } );
      } else {
        copy(W, B);
      }

      // row exchange and sort
      if (ip.dimension_0()) {
        // structure copy
        size_type offset = 0;
        for (ordinal_type i=0;i<W.NumRows();++i) {
          A.RowPtrBegin(i) = offset; 
          offset += W.NumNonZerosInRow(ip[i]);
          A.RowPtrEnd(i) = offset;
        } 
        // value copy
        Kokkos::parallel_for( range_policy_type(0, W.NumRows()),
                              [&](const ordinal_type i)
                              {
                                const auto ii = ip[i];
                                const auto W_cols = W.ColsInRow(i);
                                const auto W_vals = W.ValuesInRow(i);
                                
                                const auto A_cols = A.ColsInRow(i);
                                const auto A_vals = A.ValuesInRow(i);
                                
                                for (size_type j=0;j<W_cols.dimension_0();++j) {
                                  A_cols[j] = W_cols[j];
                                  W_cols[j] = j;  // use W as workspace of indices
                                }
                                
                                Util::sort(A_cols, W_cols);

                                for (size_type j=0;j<W_cols.dimension_0();++j) 
                                  A_vals[W_cols[j]] = W_vals[j];
                              } );
      } else {
        copy(A, W);
      }
      space_type::execution_space::fence();
    }


    // /// Flat to Hier
    // /// ------------------------------------------------------------------

    // /// Properties: 
    // /// - Compile with Device (o), 
    // /// - Callable in KokkosFunctors (o)

    // /// \brief compute dimension of hier matrix when the flat is divided by mb x nb
    // template<typename DenseMatrixFlatType,
    //          typename OrdinalType>
    // KOKKOS_INLINE_FUNCTION
    // static void
    // getDimensionOfHierMatrix(OrdinalType &hm,
    //                          OrdinalType &hn,
    //                          const DenseMatrixFlatType &flat,
    //                          const OrdinalType mb,
    //                          const OrdinalType nb) {
    //   const auto fm = flat.NumRows();
    //   const auto fn = flat.NumCols();

    //   const auto mbb = (mb == 0 ? fm : mb);
    //   const auto nbb = (nb == 0 ? fn : nb);

    //   hm = fm/mbb + (fm%mbb > 0);
    //   hn = fn/nbb + (fn%nbb > 0);
    // }
      
    // /// Properties: 
    // /// - Compile with Device (o), 
    // /// - Callable in KokkosFunctors (o)
    // /// - Blocking with fence (o)

    // /// \brief fill hier matrix 
    // template<typename DenseMatrixHierType,
    //          typename DenseMatrixFlatType,
    //          typename OrdinalType>
    // KOKKOS_INLINE_FUNCTION
    // static void
    // getHierMatrix(DenseMatrixHierType &hier,
    //               const DenseMatrixFlatType &flat,
    //               const OrdinalType mb,
    //               const OrdinalType nb,
    //               const bool create = false) {
    //   static_assert( Kokkos::Impl::is_same<
    //                  typename DenseMatrixHierType::space_type,
    //                  typename DenseMatrixFlatType::space_type
    //                  >::value,
    //                  "Space type of input matrices does not match" );
      
    //   typedef typename DenseMatrixHierType::space_type space_type;
    //   typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;

    //   OrdinalType hm, hn;
    //   getDimensionOfHierMatrix(hm, hn,
    //                            flat, 
    //                            mb, nb);
      
    //   if (create)
    //     hier.create(hm, hn);

    //   const OrdinalType fm = flat.NumRows(), fn = flat.NumCols();

    //   space_type::execution_space::fence();

    //   Kokkos::parallel_for( range_policy_type(0, hn),
    //                         [&](const ordinal_type j)
    //                         {
    //                           const OrdinalType offn = nb*j;
    //                           const OrdinalType ntmp = offn + nb; 
    //                           const OrdinalType n    = ntmp < fn ? nb : (fn - offn); 

    //                           //#pragma unroll
    //                           for (ordinal_type i=0;i<hm;++i) {
    //                             const OrdinalType offm = mb*i;
    //                             const OrdinalType mtmp = offm + mb; 
    //                             const OrdinalType m    = mtmp < fm ? mb : (fm - offm); 
    //                             hier.Value(i, j).setView(flat, offm, m,
    //                                                      /**/  offn, n);
    //                           }
    //                         } );

    //   space_type::execution_space::fence();
    // }

    // /// \brief create hier matrix 
    // template<typename DenseMatrixHierType,
    //          typename DenseMatrixFlatType,
    //          typename OrdinalType>
    // KOKKOS_INLINE_FUNCTION
    // static void
    // createHierMatrix(DenseMatrixHierType &hier,
    //                  const DenseMatrixFlatType &flat,
    //                  const OrdinalType mb,
    //                  const OrdinalType nb) {
    //   OrdinalType hm, hn;
    //   getDimensionOfHierMatrix(hm, hn, 
    //                            flat,
    //                            mb, nb);
    //   getHierMatrix(hier, flat, mb , nb, true);
    // }
    
  };

}

#endif
