#pragma once
#ifndef __SCALE_HPP__
#define __SCALE_HPP__

/// \file scale.hpp
/// \brief Scaling sparse matrix.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<typename T> struct ScaleTraits {
    typedef T scale_type;
    static T one; 
    static T zero;
  };

  // assume built-in types have appropriate type conversion
  template<typename T> T ScaleTraits<T>::one  = 1;
  template<typename T> T ScaleTraits<T>::zero = 0;
  
  
  template<typename ScalarType, 
           typename CrsMatViewType,
           typename ParallelForType>
  KOKKOS_INLINE_FUNCTION 
  int
  scale(const ParallelForType::member_type member, 
        const ScalarType alpha, 
        CrsMatViewType &A) {
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::row_view_type row_view_type;

    if (alpha == ScaleTraits<value_type>::one) { 
      // do nothing
    } else {
      row_view_type row;
      ParallelFor(ParallelFor::TeamThreadLoop(member, 0, A.NumRows()),
                  [&](const ordinal_type i) {
                    row.setView(A, i);
                    for (ordinal_type j=0;j<row.NumNonZeros();++j)
                      row.Value(j) *= alpha;
                  });
    }

    return 0;
  }

}

#endif
