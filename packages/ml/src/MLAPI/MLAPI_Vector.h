#ifndef ML_VECTOR_H
#define ML_VECTOR_H

#include "ml_config.h"
#include <iostream>
#include "MLAPI_Space.h"

namespace MLAPI {

template<class, class> class BaseObjectSum;
template<class, class> class BaseObjectDiff;
template<class, class> class BaseObjectMult;
template<class, class> class BaseObjectDiv;

class Vector {

public:
  Vector(const Space& VectorSpace)
  {
    VectorSpace_ = new Space(VectorSpace);
    MyLength_ = VectorSpace.NumMyElements();
    Values_ = new double[MyLength_];
    DeleteValues_ = true;
  }

  Vector(const Space& VectorSpace, double* Values)
  {
    VectorSpace_ = new Space(VectorSpace);
    MyLength_ = VectorSpace.NumMyElements();
    Values_ = Values;
    DeleteValues_ = false;
  }

  ~Vector() 
  {
    if (DeleteValues_)
      delete[] Values_;
  }

  Vector& operator=(double rhs) 
  {
    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = rhs;
    return(*this);
  }

  Vector& operator=(const Vector& right) 
  {
    assert (MyLength_ == right.VectorSpace().NumMyElements());
    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = right[i];
    return(*this);
  }

  template<class Left, class Right> Vector&
  operator=(const BaseObjectSum<Left,Right>& right)
  {
    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = right[i];
    return(*this);
  }

  template<class Left, class Right> Vector&
  operator=(const BaseObjectDiff<Left,Right>& right)
  {
    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = right[i];
    return(*this);
  }

  template<class Left, class Right> Vector&
  operator=(const BaseObjectMult<Left,Right>& right)
  {
    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = right[i];
    return(*this);
  }

  template<class Left, class Right> Vector&
  operator=(const BaseObjectDiv<Left,Right>& right)
  {
    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = right[i];
    return(*this);
  }

  const double& operator[](int i) const 
  {
    return Values_[i];
  }

  double& operator[] (int i) 
  {
    return(Values_[i]);
  }

  const Space& VectorSpace() const 
  {
    return(*VectorSpace_);
  }

  double DotProduct(const Vector& RHS) const 
  {
    assert (RHS.VectorSpace() == VectorSpace());

    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      MyResult += Values_[i] * RHS[i];

    VectorSpace().Comm().SumAll(&MyResult,&Result,1);
    return(Result);
  }

  double Norm2() const 
  {
    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      MyResult += Values_[i] * Values_[i];

    VectorSpace().Comm().SumAll(&MyResult,&Result,1);
    return(sqrt(Result));
  }

  double NormInf() const 
  {
    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      if (Values_[i] > MyResult)
        MyResult = Values_[i];

    VectorSpace().Comm().MaxAll(&MyResult,&Result,1);
    return(Result);
  }

private:
  int MyLength_;
  double* Values_;
  bool DeleteValues_;
  Space* VectorSpace_;

}; // Vector

} // namespace MLAPI
#endif // if ML_VECTOR_H
