#ifndef ML_VECTOR_H
#define ML_VECTOR_H

#include "ml_include.h"
#include "ml_comm.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"

namespace MLAPI {

/*!
\class DoubleVector

\brief Basic class for distributed vectors.

DoubleVector is the basic class to define distributed vectors. An example 
of usage is as follows:
\code
using namespace MLAPI;
...
int NumMyElements = 5;
Space MySpace(NumMyElements);
DoubleVector x(MySpace), y(MySpace), z(MySpace);

x = 1.0;
for (int i = 0 ; i < y.MyLength() ; ++i)
  y(i) = 1.0 * i;

z = x + y;
double NormZ = z.Norm2();
\endcode

\author Marzio Sala, SNL 9214.

\date Last updated on 07-Jan-05
*/

class DoubleVector : public BaseObject {

public:
  //! Default constructor.
  DoubleVector()
  {
    Initialize();
  }

  //! Constructor for a given Space.
  DoubleVector(const Space& VectorSpace) 
  {
    Initialize();
    Reshape(VectorSpace);
  }

  //! Copy constructor.
  DoubleVector(const DoubleVector& RHS)
  {
    exit(-12);
    Reshape(RHS.VectorSpace());

    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = RHS(i);
  }

  //! Constructor with a given Space, and user-provided array of values.
  DoubleVector(const Space& VectorSpace, double* Values)
  {
    VectorSpace_ = VectorSpace;
    MyLength_ = VectorSpace.NumMyElements();
    Values_ = Values;
    DeleteValues_ = false;
  }

  //! Destructor.
  ~DoubleVector() 
  {
    Destroy();
  }

  //! Sets the name of \c this object.
  DoubleVector& operator=(const string& Name)
  {
    SetName(Name);
    return(*this);
  }

  //! Operator=.
  DoubleVector& operator=(double rhs) 
  {
    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = rhs;
    return(*this);
  }

  //! Operator=.
  DoubleVector& operator=(const DoubleVector& RHS) 
  {
    if (this == &RHS)
      return(*this);

    if (MyLength_ != RHS.VectorSpace().NumMyElements()) {
      Reshape(RHS.VectorSpace());
    }

    for (int i = 0 ; i < MyLength_ ; ++i)
      Values_[i] = RHS(i);
    return(*this);
  }

  //! Returns the value of local element \c i.
  const double& operator()(int i) const 
  {
    return Values_[i];
  }

  //! Returns the value of local element \c i (non-const version).
  double& operator() (int i) 
  {
    return(Values_[i]);
  }

  //! Resets the Space of \this vectors.
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

  //! Returns the Space on which \c this vector is defined.
  const Space& VectorSpace() const 
  {
    return(VectorSpace_);
  }

  //! Computes the dot product between \c this vector and \c RHS.
  double DotProduct(const DoubleVector& RHS) const 
  {
    assert (RHS.VectorSpace() == VectorSpace());

    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      MyResult += Values_[i] * RHS(i);

    Result = ML_Comm_GsumDouble(GetMLComm(),MyResult);
    return(Result);
  }

  //! Computes the 2-norm of \c this vector.
  double Norm2() const 
  {
    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      MyResult += Values_[i] * Values_[i];

    Result = ML_Comm_GsumDouble(GetMLComm(),MyResult);
    return(sqrt(Result));
  }

  //! Computes the infinite norm of \c this vector.
  double NormInf() const 
  {
    double MyResult = 0.0, Result = 0.0;
    for (int i = 0 ; i < MyLength_ ; ++i) 
      if (Values_[i] > MyResult)
        MyResult = Values_[i];

    Result = ML_Comm_GmaxDouble(GetMLComm(),MyResult);
    return(Result);
  }

  //! Populates the vector with random elements.
  void Random() {
    ML_random_vec(Values_,MyLength_,MLAPI::GetMLComm());
    return;
  }

  //! Replaces each element of the vector with its reciprocal.
  void Reciprocal() 
  {
    for (int i = 0 ; i < MyLength_ ; ++i) {
      if (Values_[i] != 0.0)
        Values_[i] = 1.0 / Values_[i];
    }
  }

  //! Returns the length of \c this vector.
  int MyLength() const
  {
    return(MyLength_);
  }

private:
  //! Initialize \c this object.
  void Initialize()
  {
    MyLength_ = 0;
    Values_ = 0;
    DeleteValues_ = false;
  }

  //! Destroys data contained in \c this object.
  void Destroy() 
  {
    if (DeleteValues_)
      delete[] Values_;
    Initialize();
  }

  //! Length of the vector.
  int MyLength_;
  //! Pointer to locally own values.
  double* Values_;
  //! If \c true, Destroy() will delete \c Values.
  bool DeleteValues_;
  //! Data layout.
  Space VectorSpace_;

}; // DoubleVector

std::ostream& operator<< (std::ostream& os, const DoubleVector& v) 
{
  for (int i = 0 ; i < v.VectorSpace().NumMyElements() ; ++i)
    os << v(i) << ' ';
  os << std::endl;
  return(os);
}

} // namespace MLAPI

#endif // if ML_VECTOR_H
