#pragma once
#ifndef __DOT_HPP__
#define __DOT_HPP__

namespace Example { 

  using namespace std;
  
  template<typename CrsRowViewType>
  inline typename CrsRowViewType::real_type 
  dot(const CrsRowViewType x, const CrsRowViewType y) {
    typedef typename CrsRowViewType::ordinal_type ordinal_type;
    typedef typename CrsRowViewType::value_type   value_type;
    typedef typename CrsRowViewType::real_type    real_type;

    real_type r_val = 0.0;
    if (x.NumNonZeros() < y.NumNonZeros())  
      for (ordinal_type j=0;j<x.NumNonZeros();++j) 
        r_val += (x.Value(j) * y.get(x.Col(j)));
    else 
      for (ordinal_type j=0;j<y.NumNonZeros();++j) 
        r_val += (y.Value(j) * x.get(y.Col(j)));

    return r_val;
  }

}

#endif
