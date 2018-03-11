#pragma once
#ifndef LOWER_BANDED_MATRIX_HPP
#define LOWER_BANDED_MATRIX_HPP

#include<vector>





template<typename Real>
class LowerBandedMatrix {
  
  using std::vector;
  using rvector = vector<Real>;
  using ivector = vector<int>;

private:

  ivector         index;
  vector<rvector> value;

public:   

  LowerBandedMatrix( ivector band_index, 
                     vector<rvector> band_value ) : 
    index(band_index), value(band_value) {
  }

  void apply( rvector& Ax, const rvector& x ) {

  }

  void solve( rvector& x, const rvector& Ax ) {

  }
  


}; // class LowerBandedMatrix

#endif // LOWER_BANDED_MATRIX_HPP

