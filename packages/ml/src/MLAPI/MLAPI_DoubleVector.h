#ifndef ML_VECTOR_H
#define ML_VECTOR_H

#include "ml_include.h"
#include "ml_comm.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
//#include "MLAPI_Workspace.h"

namespace MLAPI {

template<class, class> class BaseObjectSum;
template<class, class> class BaseObjectDiff;
template<class, class> class BaseObjectMult;
template<class, class> class BaseObjectDiv;

class DoubleVector : public BaseObject {

public:
  DoubleVector()
  {
    Initialize();
  }

  DoubleVector(const Space& VectorSpace) 
  {
    Initialize();
    Reshape(VectorSpace);
  }

  DoubleVector(const DoubleVector& RHS)
  {
    exit(-12);
    Reshape(RHS.VectorSpace());

    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = RHS(i);
  }

  DoubleVector(const Space& VectorSpace, double* Values)
  {
    VectorSpace_ = VectorSpace;
    MyLength_ = VectorSpace.NumMyElements();
    Values_ = Values;
    DeleteValues_ = false;
  }

  ~DoubleVector() 
  {
    Destroy();
    Initialize();
  }

  void Initialize()
  {
    MyLength_ = 0;
    Values_ = 0;
    DeleteValues_ = false;
  }

  void Destroy() 
  {
    if (DeleteValues_)
      delete[] Values_;
    Initialize();
  }

  void Reshape(const Space& RHS)
  {
    Destroy();

    MyLength_ = RHS.NumMyElements();
    if (MyLength_) {
      VectorSpace_ = RHS;
      Values_ = new double[MyLength_];
      DeleteValues_ = true;
    }
  }

  DoubleVector& operator=(const string& Name)
  {
    SetName(Name);
    return(*this);
  }

  DoubleVector& operator=(double rhs) 
  {
    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = rhs;
    return(*this);
  }

  DoubleVector& operator=(const DoubleVector& RHS) 
  {
    if (this == &RHS)
      return(*this);

    if (MyLength_ != RHS.VectorSpace().NumMyElements()) {
      assert(1 == 0);
      Reshape(RHS.VectorSpace());
    }

    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = RHS(i);
    return(*this);
  }

  template<class Left, class Right> DoubleVector&
  operator=(const BaseObjectSum<Left,Right>& RHS)
  {
    if (MyLength_ != RHS.VectorSpace().NumMyElements()) {
      Destroy();
      Reshape(RHS.VectorSpace());
    }

    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = RHS(i);
    return(*this);
  }

  template<class Left, class Right> DoubleVector&
  operator=(const BaseObjectDiff<Left,Right>& RHS)
  {
    if (MyLength_ != RHS.VectorSpace().NumMyElements()) {
      Destroy();
      Reshape(RHS.VectorSpace());
    }

    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = RHS(i);
    return(*this);
  }

  template<class Left, class Right> DoubleVector&
  operator=(const BaseObjectMult<Left,Right>& RHS)
  {
    if (MyLength_ != RHS.VectorSpace().NumMyElements()) {
      Destroy();
      Reshape(RHS.VectorSpace());
    }

    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = RHS(i);
    return(*this);
  }

  template<class Left, class Right> DoubleVector&
  operator=(const BaseObjectDiv<Left,Right>& RHS)
  {
    if (MyLength_ != RHS.VectorSpace().NumMyElements()) {
      Destroy();
      Reshape(RHS.VectorSpace());
    }

    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = RHS(i);
    return(*this);
  }

  const double& operator()(int i) const 
  {
    return Values_[i];
  }

  double& operator() (int i) 
  {
    return(Values_[i]);
  }

  const Space& VectorSpace() const 
  {
    return(VectorSpace_);
  }

  double DotProduct(const DoubleVector& RHS) const 
  {
    assert (RHS.VectorSpace() == VectorSpace());

    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      MyResult += Values_[i] * RHS(i);

    Result = ML_Comm_GsumDouble(GetMLComm(),MyResult);
    return(Result);
  }

  double Norm2() const 
  {
    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      MyResult += Values_[i] * Values_[i];

    Result = ML_Comm_GsumDouble(GetMLComm(),MyResult);
    return(sqrt(Result));
  }

  double NormInf() const 
  {
    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      if (Values_[i] > MyResult)
        MyResult = Values_[i];

    Result = ML_Comm_GmaxDouble(GetMLComm(),MyResult);
    return(Result);
  }

  void Random() {
    ML_random_vec(Values_,MyLength_,MLAPI::GetMLComm());
    return;
  }

  void Reciprocal() 
  {
    for (int i = 0 ; i < MyLength_ ; ++i) {
      if (Values_[i] != 0.0)
        Values_[i] = 1.0 / Values_[i];
    }
  }

  int MyLength() const
  {
    return(MyLength_);
  }

private:
  int MyLength_;
  double* Values_;
  bool DeleteValues_;
  Space VectorSpace_;

}; // DoubleVector

} // namespace MLAPI

#endif // if ML_VECTOR_H
