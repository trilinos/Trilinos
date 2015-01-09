#pragma once
#ifndef __DOT_HPP__
#define __DOT_HPP__

/// \file dot.hpp
/// \brief Sparse dot product.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Example { 

  using namespace std;

  template<typename T> struct DotTraits {
    typedef T dot_type;

    static KOKKOS_FORCEINLINE_FUNCTION 
    dot_type 
    dot(const T &x, const T &y) { return conj<T>(x)*y; }
  }; 

  template<typename CrsRowViewType>
  KOKKOS_INLINE_FUNCTION 
  typename CrsRowViewType::value_type
  dot(const CrsRowViewType x, const CrsRowViewType y) {
    typedef typename CrsRowViewType::ordinal_type ordinal_type;
    typedef typename CrsRowViewType::value_type   value_type;

    typedef DotTraits<value_type> dot_traits;
    
    value_type r_val(0);
    if (x.NumNonZeros() < y.NumNonZeros())  
      for (ordinal_type j=0;j<x.NumNonZeros();++j) 
        r_val += dot_traits::dot(x.Value(j), y.ValueAtColumn(x.Col(j)));
    else 
      for (ordinal_type j=0;j<y.NumNonZeros();++j) 
        r_val += dot_traits::dot(y.Value(j), x.ValueAtColumn(y.Col(j)));

    return r_val;
  }

}

#endif
