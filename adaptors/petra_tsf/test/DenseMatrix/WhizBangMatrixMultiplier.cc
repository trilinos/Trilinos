
#include "WhizBangMatrixMultiplier.h"

// Constuctor takes two matrices. Multiplier multiplies them together, result in M1
template<class scalarType>
WhizBangMatrixMultiplier<scalarType>::WhizBangMatrixMultiplier(TSF::DenseMatrix<scalarType> * M1, 
                                                               TSF::DenseMatrix<scalarType> * M2, 
                                                               TSF::DenseMatrix<scalarType> * M3):
  M1_(M1),
  M2_(M2),
  M3_(M3)
{}
  
WhizBangMatrixMultiplier<scalarType>::~WhizBangMatrixMultiplier() {}
  
void WhizBangMatrixMultiplier<scalarType>::multiply() { 

  M3_->multiply('N', 'T', 1.0, *M1_, *M2_, 0.0);
  return;
}
