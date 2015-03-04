#pragma once
#ifndef __HERK_U_T_HPP__
#define __HERK_U_T_HPP__

/// \file herk_u_t.hpp
/// \brief Sparse hermitian rank one update on given sparse patterns.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;
  
  template<>
  template<typename ScalarType, 
           typename CrsMatViewType,
           typename ParallelForType>
  KOKKOS_INLINE_FUNCTION 
  int
  Herk<Uplo::Upper,Trans::ConjTranspose,
       AlgoHerk::ForRightBlocked>
  ::invoke(const ParallelForType::member_type member,
           const ScalarType alpha,
           const CrsMatViewType A,
           const ScalarType beta,
           const CrsMatViewType C) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;
    typedef DotTraits<value_type> dot_traits;

    scale<ScalarType,CrsMatViewType,ParallelForType>(member, beta, C);
    
    // row_view_type a, c;
    for (ordinal_type k=0;k<A.NumRows();++k) {
      //a.setView(A, k);
      row_view_type &a = A.RowView(k);
      const ordinal_type nnz = a.NumNonZeros();

      ParallelFor(ParallelFor::TeamThreadLoop(member, 0, nnz),
                  [&](const ordinal_type i) {
                    const ordinal_type row_at_i = a.Col(i);
                    const value_type   val_at_i = conj(a.Value(i));
                    
                    //c.setView(C, row_at_i);
                    row_view_type &c = C.RowView(row_at_i);
                    
                    ordinal_type idx = 0;
                    for (ordinal_type j=i;j<nnz && (idx > -2);++j) {
                      ordinal_type col_at_j = a.Col(j);
                      value_type   val_at_j = a.Value(j);
                      
                      idx = c.Index(col_at_j, idx);
                      if (idx >= 0) 
                        c.Value(idx) += alpha*val_at_i*val_at_j;
                    }
                  });
    }

    return 0;
  }

}

#endif
