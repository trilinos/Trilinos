#pragma once
#ifndef __SCALE_HPP__
#define __SCALE_HPP__

/// \file scale.hpp
/// \brief Scaling sparse matrix.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;
  
  template<typename ScalarType, 
           typename CrsMatViewType>
  KOKKOS_INLINE_FUNCTION 
  int
  scale(const ScalarType alpha,
        CrsMatViewType &A) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;

    for (ordinal_type i=0;i<A.NumRows();++i) {
      row_view_type row = A.extractRow(i);
      for (ordinal_type j=0;j<row.NumNonZeros();++j)
        row.Value(j) *= alpha;
    }

    return 0;
  }

}

#endif
