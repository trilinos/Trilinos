#include <cstdlib>
#include <cassert>
#include <iostream>
#include <strstream>
#include <cstring>
#include <cmath>
#include <iomanip>
using namespace std;

#include "TSF_DenseMatrix.h"
#include "TPetra_DenseMatrix.h"

// Constuctor takes three matrices. Multiplier multiplies them together, result in M3
template<class scalarType>
class WhizBangMatrixMultiplier {

public:

  WhizBangMatrixMultiplier(TSF::DenseMatrix<scalarType> * M1, 
                           TSF::DenseMatrix<scalarType> * M2, 
                           TSF::DenseMatrix<scalarType> * M3);
  
    ~WhizBangMatrixMultiplier();
  
    void multiply();  // Multiplies M1 and M2, stores result in M3
  
private:
  
  TSF::DenseMatrix<scalarType> * M1_;  
  TSF::DenseMatrix<scalarType> * M2_;
  TSF::DenseMatrix<scalarType> * M3_;

};
