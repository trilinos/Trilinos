#ifndef __TACHO_CRS_MATRIX_TOOLS_HPP__
#define __TACHO_CRS_MATRIX_TOOLS_HPP__

/// \file Tacho_CrsMatrixTools.hpp
/// \brief Generic utility function for crs matrices.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

  class CrsMatrixTools {
  public:

    // Tim Davis, Direct Methods for Sparse Linear Systems, Siam, p 42.
    template<typename CrsMatBaseType,
             typename OrdinalTypeArray>
    static void
    createEliminationTree(OrdinalTypeArray &parent,
                          const CrsMatBaseType A) {
      const auto m = A.NumRows();
      auto ancestor = OrdinalTypeArray("CrsMatrixTools::createEliminationTree::ancestor", m);

      if (parent.data() == NULL || parent.dimension(0) < m)
        parent = OrdinalTypeArray("CrsMatrixTools::createEliminationTree::parent", m);

      for (auto i=0;i<m;++i) {
        parent(i) = -1;
        ancestor(i) = -1;
        
        const auto nnz_in_row = A.NumNonZerosInRow(i);
        const auto cols_in_row = A.ColsInRow(i);        

        for (auto idx=0;idx<nnz_in_row;++idx) 
          for (auto j=cols_in_row(idx);j != -1 && j < i;) {
            const auto next = ancestor(j);
            ancestor(j) = i ;
            if (next == -1) parent(j) = i;
            j = next;
          }
      }
    }

    // Tim Davis, Algorithm 849: A Concise Sparse Cholesky Factorization Package
    template<typename CrsMatBaseType,
             typename OrdinalTypeArray,
             typename SizeTypeArray>
    static void
    createSymbolicStructure(SizeTypeArray &ap,  // row pointer of upper triangular
                            OrdinalTypeArray &aj, // col indices
                            const CrsMatBaseType A) {
      const auto m = A.NumRows();
      auto parent = OrdinalTypeArray("CrsMatrixTools::createEliminationTree::parent", m);
      auto flag = OrdinalTypeArray("CrsMatrixTools::createEliminationTree::flag", m);
      auto upper_row_cnt = OrdinalTypeArray("CrsMatrixTools::createEliminationTree::upper_row_cnt", m);

      for (auto i=0;i<m;++i) {
        parent(i) = -1;
        flag(i) = -1;
        upper_row_cnt(i) = 1; // add diagonal entry

        const auto nnz_in_row = A.NumNonZerosInRow(i);
        const auto cols_in_row = A.ColsInRow(i);        
        
        for (auto idx=0;idx<nnz_in_row;++idx) 
          for (auto j=cols_in_row(idx) ; flag(j) != i && j < i; j = parent(j) ) {
            if (parent(j) == -1) parent(j) = i;
            ++upper_row_cnt(j);
            flag(j) = i ;
          }
      }

      // prefix scan
      ap = SizeTypeArray("CrsMatrixTools::createEliminationTree::ap", m+1);
      for (auto i=0;i<m;++i) 
        ap(i+1) = ap(i) + upper_row_cnt(i);

      // detect aj
      aj = OrdinalTypeArray("CrsMatrixTools::createEliminationTree::aj", ap(m));
      for (auto i=0;i<m;++i) {
        parent(i) = -1;
        flag(i) = -1;

        // diagonal entry
        aj(ap(i)) = i; 
        upper_row_cnt(i) = 1;

        const auto nnz_in_row = A.NumNonZerosInRow(i);
        const auto cols_in_row = A.ColsInRow(i);        
        
        for (auto idx=0;idx<nnz_in_row;++idx) 
          for (auto j=cols_in_row(idx) ; flag(j) != i && j < i; j = parent(j) ) {
            if (parent(j) == -1) parent(j) = i;
            aj(ap(j)+upper_row_cnt(j)) = i;
            ++upper_row_cnt(j);
            flag(j) = i ;
          }
      }
    }

    /// Elementwise copy
    /// ------------------------------------------------------------------

    /// \brief sort columns
    template<typename CrsMatBaseType>
    static void
    sortColumnsPerRow(CrsMatBaseType &A) {
      typedef typename CrsMatBaseType::ordinal_type       ordinal_type;
      typedef typename CrsMatBaseType::size_type          size_type;
      typedef typename CrsMatBaseType::ordinal_type_array ordinal_type_array;
      typedef typename CrsMatBaseType::value_type_array   value_type_array;

      typedef typename CrsMatBaseType::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;

      const size_type nnz = A.NumNonZeros();
      ordinal_type_array idx(Kokkos::ViewAllocateWithoutInitializing("CrsMatrixTools::sortColumns::idx"), nnz);
      value_type_array   tmp(Kokkos::ViewAllocateWithoutInitializing("CrsMatrixTools::sortColumns::tmp"), nnz);

      Kokkos::deep_copy(tmp, A.Values());

      space_type::execution_space::fence();

      // assume that rowpt begin and end arrays are separated.
      Kokkos::parallel_for( range_policy_type(0, A.NumRows()),
                            [&](const ordinal_type i)
                            {
                              const size_type
                                row_begin  = A.RowPtrBegin(i),
                                row_end    = A.RowPtrEnd(i),
                                sort_begin = 0,
                                nnz_in_row = (row_end - row_begin);

                              // setup index in row
                              Kokkos::pair<size_type,size_type> range(row_begin, row_end);
                              auto idx_in_row = Kokkos::subview(idx, range);
                              for (auto k=0;k<nnz_in_row;++k)
                                idx_in_row(k) = k;

                              // sort columnes and index
                              auto cols_in_row = A.ColsInRow(i);
                              Util::sort(cols_in_row, idx_in_row, sort_begin, nnz_in_row);

                              // apply index to values
                              auto vals_in_row = A.ValuesInRow(i);
                              auto tmp_in_row = Kokkos::subview(tmp, range);                                
                              for (auto k=0;k<nnz_in_row;++k) 
                                vals_in_row(k) = tmp_in_row(idx_in_row(k));
                            } );
    }

    /// \brief elementwise copy
    template<typename CrsMatBaseTypeA,
             typename CrsMatBaseTypeB>
    static void
    copy(CrsMatBaseTypeA &A,
         const CrsMatBaseTypeB B) {
      static_assert( Kokkos::Impl::is_same<
                     typename CrsMatBaseTypeA::space_type,
                     typename CrsMatBaseTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      typedef typename CrsMatBaseTypeA::ordinal_type ordinal_type;
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
    static void
    copy(CrsMatBaseTypeA &A,
         const int uplo,
         const int offset,
         const CrsMatBaseTypeB B) {
      static_assert( Kokkos::Impl
                     ::is_same<
                     typename CrsMatBaseTypeA::space_type,
                     typename CrsMatBaseTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      typedef typename CrsMatBaseTypeA::ordinal_type ordinal_type;
      typedef typename CrsMatBaseTypeA::size_type size_type;
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
          const auto ioffset = i + offset;
          A.RowPtrBegin(i) = nnz;
          for (ordinal_type idx=0;idx<cols.extent(0);++idx) {
            if (ioffset <= cols(idx)) {
              A.Col(nnz) = cols(idx);
              A.Value(nnz) = vals(idx);
              ++nnz;
            }
          }
          A.RowPtrEnd(i) = nnz;
        }
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
          const auto ioffset = i - offset;
          A.RowPtrBegin(i) = nnz;
          for (ordinal_type idx=0;idx<cols.extent(0);++idx) {
            if (ioffset >= cols(idx)) {
              A.Col(nnz) = cols(idx);
              A.Value(nnz) = vals(idx);
              ++nnz;
            }
          }
          A.RowPtrEnd(i) = nnz;
        }
        break;
      }
      }
      //space_type::execution_space::fence();
      A.setNumNonZeros();
    }

    /// \brief elementwise copy with permutation
    template<typename CrsMatBaseTypeA,
             typename CrsMatBaseTypeB,
             typename OrdinalTypeArray>
    static void
    copy(CrsMatBaseTypeA &A,
         const OrdinalTypeArray p,
         const OrdinalTypeArray ip,
         const CrsMatBaseTypeB B) {
      static_assert( Kokkos::Impl
                     ::is_same<
                     typename CrsMatBaseTypeA::space_type,
                     typename CrsMatBaseTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      typedef typename CrsMatBaseTypeA::ordinal_type ordinal_type;
      typedef typename CrsMatBaseTypeA::size_type size_type;
      typedef typename CrsMatBaseTypeA::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;

      // create work space
      CrsMatBaseTypeA W("CrsMatrixTools::copy::W");
      W.createConfTo(A);

      space_type::execution_space::fence();

      // column exchange
      if (p.extent(0)) {
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

                                for (size_type j=0;j<B_cols.extent(0);++j) {
                                  W_cols(j) = p(B_cols(j));
                                  W_vals(j) = B_vals(j);
                                }
                              } );
      } else {
        copy(W, B);
      }

      // row exchange and sort
      if (ip.extent(0)) {
        // structure copy
        size_type offset = 0;
        for (ordinal_type i=0;i<W.NumRows();++i) {
          A.RowPtrBegin(i) = offset;
          offset += W.NumNonZerosInRow(ip(i));
          A.RowPtrEnd(i) = offset;
        }
        // value copy
        Kokkos::parallel_for( range_policy_type(0, W.NumRows()),
                              [&](const ordinal_type i)
                              {
                                const auto ii = ip(i);
                                const auto W_cols = W.ColsInRow(ii);
                                const auto W_vals = W.ValuesInRow(ii);

                                const auto A_cols = A.ColsInRow(i);
                                const auto A_vals = A.ValuesInRow(i);

                                for (size_type j=0;j<W_cols.extent(0);++j) {
                                  A_cols(j) = W_cols(j);
                                  //A_vals(j) = W_vals(j);
                                  W_cols(j) = j;  // use W as workspace of indices
                                }

                                const ordinal_type begin = 0, end = A_cols.extent(0);
                                Util::sort(A_cols, W_cols, begin, end);

                                for (size_type j=begin;j<end;++j)
                                  A_vals(j) = W_vals(W_cols(j));
                              } );
      } else {
        copy(A, W);
      }

      space_type::execution_space::fence();
      A.setNumNonZeros();
    }


    /// \brief elementwise add
    template<typename CrsMatBaseTypeA,
             typename CrsMatBaseTypeB>
    static void
    add(CrsMatBaseTypeA &A,
        const CrsMatBaseTypeB B) {
      static_assert( Kokkos::Impl
                     ::is_same<
                     typename CrsMatBaseTypeA::space_type,
                     typename CrsMatBaseTypeB::space_type
                     >::value,
                     "Space type of input matrices does not match" );
      typedef typename CrsMatBaseTypeA::ordinal_type ordinal_type;
      typedef typename CrsMatBaseTypeA::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;

      const ordinal_type m = Util::min(A.NumRows(), B.NumRows());

      space_type::execution_space::fence();
      Kokkos::parallel_for( range_policy_type(0, m),
                            [&](const ordinal_type i)
                            {
                              auto ja = A.RowPtrBegin(i);
                              auto jb = B.RowPtrBegin(i);

                              const auto jaend = A.RowPtrEnd(i);
                              const auto jbend = B.RowPtrEnd(i);

                              for ( ;jb<jbend;++jb) {
                                for ( ;(A.Col(ja)<B.Col(jb) && ja<jaend);++ja);
                                A.Value(ja) += (A.Col(ja) == B.Col(jb))*B.Value(jb);
                              }
                            } );
      space_type::execution_space::fence();
    }

    /// \brief elementwise add
    template<typename CrsMatBaseType,
             typename ScaleVectorType>
    static void
    computeRowScale(ScaleVectorType &s,
                    const CrsMatBaseType &A) {
      typedef typename CrsMatBaseType::ordinal_type ordinal_type;
      typedef typename CrsMatBaseType::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;
      
      const ordinal_type m = A.NumRows();
      
      space_type::execution_space::fence();
      Kokkos::parallel_for( range_policy_type(0, m),
                            [&](const ordinal_type i)
                            {
                              auto j = A.RowPtrBegin(i);
                              
                              const auto jend = A.RowPtrEnd(i);

                              // for complex abs should return magnitude type
                              typename ScaleVectorType::value_type sum = 0.0;
                              for ( ;j<jend;++j) 
                                sum += Util::abs(A.Value(j));

                              s(i) = sum;
                            } );
      space_type::execution_space::fence();
    }

    /// \brief elementwise add
    template<typename CrsMatBaseType,
             typename ScaleVectorType>
    static void
    applyScale(CrsMatBaseType &A,
               const ScaleVectorType &s) {
      typedef typename CrsMatBaseType::ordinal_type ordinal_type;
      typedef typename CrsMatBaseType::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;
      
      const ordinal_type m = A.NumRows();
      
      space_type::execution_space::fence();
      Kokkos::parallel_for( range_policy_type(0, m),
                            [&](const ordinal_type i)
                            {
                              auto j = A.RowPtrBegin(i);
                              
                              const auto jend = A.RowPtrEnd(i);

                              for ( ;j<jend;++j) 
                                A.Value(j) *= s(i);
                            } );
      space_type::execution_space::fence();
    }

    /// \brief elementwise add
    template<typename CrsMatBaseType,
             typename ScaleVectorType>
    static void
    applyInverseScale(CrsMatBaseType &A,
                      const ScaleVectorType &s) {
      typedef typename CrsMatBaseType::ordinal_type ordinal_type;
      typedef typename CrsMatBaseType::space_type space_type;
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy_type;
      
      const ordinal_type m = A.NumRows();
      
      space_type::execution_space::fence();
      Kokkos::parallel_for( range_policy_type(0, m),
                            [&](const ordinal_type i)
                            {
                              auto j = A.RowPtrBegin(i);
                              
                              const auto jend = A.RowPtrEnd(i);

                              for ( ;j<jend;++j) 
                                A.Value(j) /= s(i);
                            } );
      space_type::execution_space::fence();
    }


    /// \brief filter zero
    template<typename CrsMatBaseTypeA>
    static void
    filterZeros(CrsMatBaseTypeA &A) {
      typedef typename CrsMatBaseTypeA::size_type size_type;

      size_type nnz = 0, nz = 0;
      const auto mA = A.NumRows();

      for (auto i=0;i<mA;++i) {
        const auto nnz_in_row  = A.NumNonZerosInRow(i);
        const auto cols_in_row = A.ColsInRow(i);
        const auto vals_in_row = A.ValuesInRow(i);

        A.RowPtrBegin(i) = nnz;
        for (size_type j=0;j<nnz_in_row;++j) {
          const auto col = cols_in_row(j);
          const auto val = vals_in_row(j);

          if (Util::abs(val) == 0) {
            ++nz;
          } else {
            A.Col(nnz)   = col;
            A.Value(nnz) = val;
            ++nnz;
          }
        }
        A.RowPtrEnd(i) = nnz;
      }
      A.setNumNonZeros();
    }

    /// Flat to Hier
    /// ------------------------------------------------------------------

    /// \brief create hier matrix
    template<typename CrsMatrixHierType,
             typename CrsMatrixFlatType,
             typename OrdinalType,
             typename OrdinalTypeArray>
    static void
    createHierMatrix(/**/  CrsMatrixHierType &hier,
                     const CrsMatrixFlatType flat,
                     const OrdinalType       nblks,
                     const OrdinalTypeArray  range,
                     const OrdinalTypeArray  tree,
                     const bool full = false) {
      // this strictly requires disjoint tree of children (sometimes scotch does not return it)
      typedef typename CrsMatrixHierType::size_type size_type;
      typedef typename CrsMatrixHierType::ordinal_type ordinal_type;

      OrdinalTypeArray nnz_in_row("nnz_in_row", nblks);
      for (auto i=0;i<nblks;++i)
        for (auto j=i;j != -1;j=tree(j)) {
          ++nnz_in_row(i);
          if (full && i != j)
            ++nnz_in_row(j);
        }
      
      size_type nnz = 0;
      for (auto i=0;i<nblks;++i) 
        nnz += nnz_in_row(i);
      
      hier.create(nblks, nblks, nnz);

      nnz = 0;
      for (ordinal_type i=0;i<nblks;++i) {
        hier.RowPtrBegin(i) = nnz; 
        nnz += nnz_in_row(i);
        hier.RowPtrEnd(i) = nnz;
        nnz_in_row(i) = 0;
      }

      for (ordinal_type i=0;i<nblks;++i) {
        for (ordinal_type j=i;j != -1;j=tree(j)) {
          {
            const auto idx = nnz_in_row(i);
            
            auto cols = hier.ColsInRow(i);
            cols(idx) = j;

            auto vals = hier.ValuesInRow(i);
            vals(idx).setView(flat, range(i), (range(i+1) - range(i)),
                              /**/  range(j), (range(j+1) - range(j)));

            ++nnz_in_row(i);
          }
          if (full && i != j) {
            const auto idx = nnz_in_row(j);
            
            auto cols = hier.ColsInRow(j);
            cols(idx) = i;
            
            auto vals = hier.ValuesInRow(j);
            vals(idx).setView(flat, range(j), (range(j+1) - range(j)),
                              /**/  range(i), (range(i+1) - range(i)));

            ++nnz_in_row(j);
          }
        }
      }
      hier.setNumNonZeros();
    }

    /// \brief filter zero
    template<typename CrsMatBaseTypeA>
    static void
    filterEmptyBlocks(CrsMatBaseTypeA &A) {
      typedef typename CrsMatBaseTypeA::ordinal_type ordinal_type;
      typedef typename CrsMatBaseTypeA::size_type size_type;
      size_type nnz = 0, emptyblocks = 0;

      const auto mA = A.NumRows();
      for (ordinal_type i=0;i<mA;++i) {
        const auto nnz_in_row  = A.NumNonZerosInRow(i);
        const auto cols_in_row = A.ColsInRow(i);
        const auto vals_in_row = A.ValuesInRow(i);

        A.RowPtrBegin(i) = nnz;
        for (ordinal_type j=0;j<nnz_in_row;++j) {
          const auto col   = cols_in_row(j);
          const auto block = vals_in_row(j);

          if (block.isEmpty()) {
            ++emptyblocks;
          } else {
            A.Col(nnz)   = col;
            A.Value(nnz) = block;
            ++nnz;
          }
        }
        A.RowPtrEnd(i) = nnz;
      }
      A.setNumNonZeros();
    }


    // this is probably too expensive, need to find a better robust way
    // /// \brief create hier matrix
    // template<typename CrsMatrixHierType,
    //          typename CrsMatrixFlatType,
    //          typename OrdinalType,
    //          typename OrdinalTypeArray>
    // static void
    // createHierMatrix(/**/  CrsMatrixHierType &hier,
    //                  const CrsMatrixFlatType flat,
    //                  const OrdinalType       nblks,
    //                  const OrdinalTypeArray  range) {

    //   typedef typename CrsMatrixHierType::size_type size_type;
    //   typedef typename CrsMatrixHierType::ordinal_type ordinal_type;

    //   typedef typename CrsMatrixFlatType::ordinal_type_array flat_ordinal_type_array;

    //   // expensive but does not miss anything
    //   flat_ordinal_type_array ftmp("ftmp", flat.NumCols());

    //   size_type nnz = 0;
    //   for (ordinal_type i=0;i<nblks;++i) {
    //     // mark id for the range of columns
    //     const auto id = (i+1);
    //     const auto range_at_i = range(i+1) - range(i);
    //     for (ordinal_type ii=range(i);ii<range(i+1);++ii) {
    //       const auto cols_in_row = flat.ColsInRow(ii);
    //       const auto nnz_in_row  = flat.NumNonZerosInRow(ii);
    //       for (ordinal_type jj=0;jj<nnz_in_row;++jj)
    //         ftmp(cols_in_row(jj)) = id;
    //     }

    //     // count marked id for each range
    //     for (ordinal_type j=i;j<nblks;++j) {
    //       bool marked = false;
    //       for (ordinal_type jj=range(j);jj<range(j+1);++jj) {
    //         if (ftmp(jj) == id) {
    //           marked = true;
    //           break;
    //         }
    //       }
    //       if (marked) {
    //         hier.Col(nnz) = j;
    //         hier.Value(nnz).setView(flat, range(i), (range(i+1) - range(i)),
    //                                 /**/  range(j), (range(j+1) - range(j)));
    //         ++nnz;
    //       }
    //     }
    //   }
    //   hier.setNumNonZeros();
    // }

  };

}

#endif
