
#include "Petra_RDP_DenseVector.h"
//=============================================================================
Petra_RDP_DenseVector::Petra_RDP_DenseVector(void)
  : Petra_RDP_DenseMatrix()
{
}

//=============================================================================
Petra_RDP_DenseVector::Petra_RDP_DenseVector(Petra_DataAccess CV, double *Values, int Length)
  : Petra_RDP_DenseMatrix(CV, Values, Length, Length, 1)
{}
//=============================================================================
Petra_RDP_DenseVector::Petra_RDP_DenseVector(const Petra_RDP_DenseVector& Source)
  : Petra_RDP_DenseMatrix(Source)
{}
//=============================================================================
Petra_RDP_DenseVector::~Petra_RDP_DenseVector()
{}
//=========================================================================
double& Petra_RDP_DenseVector::operator() (int Index)  {

  if (Index>=M_) {
    cout << "Index = " << Index << " Out of Range 0 - " << M_-1;
    abort();
  }
   return(A_[Index]);
}

//=========================================================================
const double& Petra_RDP_DenseVector::operator() (int Index) const  {

  if (Index>=M_) {
    cout << "Index = " << Index << " Out of Range 0 - " << M_-1;
    abort();
  }
   return(A_[Index]);
}
//=========================================================================
const double& Petra_RDP_DenseVector::operator [] (int Index) const  {

#ifdef PETRA_ARRAY_BOUNDS_CHECK
  if (Index>=M_) {
    cout << "Index = " << Index << " Out of Range 0 - " << M_-1;
    abort();
  }
#endif
   return(A_[Index]);
}
//=========================================================================
double& Petra_RDP_DenseVector::operator [] (int Index)  {

#ifdef PETRA_ARRAY_BOUNDS_CHECK
  if (Index>=M_) {
    cout << "Index = " << Index << " Out of Range 0 - " << M_-1;
    abort();
  }
#endif
   return(A_[Index]);
}


// Non-member functions

ostream& operator << (ostream& os, const Petra_RDP_DenseVector& A)
{
  const int M = A.M();
  const int N = A.N();
  long olda = os.setf(ios::right,ios::adjustfield);
  long oldf = os.setf(ios::scientific,ios::floatfield);
  int oldp = os.precision(12);

  for (int i=0; i<M; i++) {
      os.width(20);
      os << A(i);
  }
    os << endl;

  // Reset os flags

  os.setf(olda,ios::adjustfield);
  os.setf(oldf,ios::floatfield);
  os.precision(oldp);

  return os;
}
