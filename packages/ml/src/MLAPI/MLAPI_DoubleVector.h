#ifndef ML_VECTOR_H
#define ML_VECTOR_H

#include "ml_include.h"
#include "ml_lapack.h"
#include "ml_comm.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
#include "Teuchos_RefCountPtr.hpp"
#include <iomanip>

namespace MLAPI {

/*!
\class DoubleVector

\brief Basic class for distributed double-precision vectors.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

class DoubleVector : public BaseObject {

public:

  //@{ Constructors and destructors

  //! Default constructor.
  DoubleVector() { }

  //! Constructor for a given Space.
  DoubleVector(const Space& VectorSpace, bool SetToZero = true)
  {
    VectorSpace_ = VectorSpace;
    if (MyLength()) {
      SetRCPValues(Teuchos::rcp(new double[MyLength()]));
      *this = 0.0;
    }
  }

  //! Constructor with a given Space, and user-provided array of values.
  DoubleVector(const Space& VectorSpace, double* Values)
  {
    VectorSpace_ = VectorSpace;
    SetRCPValues(Teuchos::rcp(Values, false));
  }

  //! Copy constructor.
  DoubleVector(const DoubleVector& rhs)
  {
    VectorSpace_ = rhs.VectorSpace();
    SetRCPValues(rhs.RCPValues());
  }

  //! Destructor.
  ~DoubleVector() 
  {
    SetRCPValues(Teuchos::null);
  }

  // @}
  // @{ Overloaded operators

  //! Sets all elements of this vector to \c rhs.
  DoubleVector& operator=(double rhs) 
  {
    for (int i = 0 ; i < MyLength() ; ++i)
      Values()[i] = rhs;

    return(*this);
  }

  //! Copies the \c rhs into \c this object.
  DoubleVector& operator=(const DoubleVector& rhs) 
  {
    if (this != &rhs) {
      VectorSpace_ = rhs.VectorSpace();
      SetRCPValues(rhs.RCPValues());
      SetLabel(rhs.GetLabel());
    }

    return(*this);
  }

  //! Sets the name of \c this object, does not touch vector elements or space.
  DoubleVector& operator=(const string& Label)
  {
    SetLabel(Label);
    return(*this);
  }

  //! Returns the value of local element \c i.
  inline const double& operator()(int i) const 
  {
#ifdef MLAPI_CHECK
    if (i < 0) || (i >= MyLength())
      ML_THROW("Requested component " + toString(i) +
               ", while MyLength() = " + toString(MyLength()), -1);
#endif
    return(Values()[i]);
  }

  //! Returns the value of local element \c i (non-const version).
  inline double& operator() (int i) 
  {
#ifdef MLAPI_CHECK
    if ((i < 0) || (i >= MyLength()))
      ML_THROW("Requested component " + toString(i) +
               ", while MyLength() = " + toString(MyLength()), -1);
#endif
    return(Values()[i]);
  }

  // @}
  // @{ Mathematical methods
  
  //! Sets the space of this vector.
  void Reshape(const Space& S)
  {
    VectorSpace_ = S;
    if (MyLength())
      SetRCPValues(Teuchos::rcp(new double[MyLength()]));
    else
      SetRCPValues(Teuchos::null);
  }

  //! Sets this = rhs.
  void Update(const DoubleVector& rhs)
  {
    int n = MyLength();
    if (n == 0) return;

    CheckSpaces(rhs);

    int incr = 1;
    // copy rhs into this
    DCOPY_F77(&n, (double*)rhs.Values(), &incr, Values(), &incr);
  }
  
  //! Sets this = alpha * rhs.
  void Update(double alpha, const DoubleVector& rhs)
  {
    int n = MyLength();
    if (n == 0) return;

    CheckSpaces(rhs);

    int incr = 1;
    // copy rhs into this
    DCOPY_F77(&n, (double*)rhs.Values(), &incr, Values(), &incr);
    // scale this
    DSCAL_F77(&n, &alpha, Values(), &incr);
  }

  //! Sets this = alpha * x + beta * y.
  void Update(double alpha, const DoubleVector& x,
              double beta,  const DoubleVector& y)
  {
    int n = MyLength();
    if (n == 0) return;

    CheckSpaces(x);
    CheckSpaces(y);

    int incr = 1;
    // copy rhs into this
    DCOPY_F77(&n, (double*)x.Values(), &incr, Values(), &incr);
    Update(beta,y,alpha);

  }

  //! Sets this = alpha * rhs + beta * this.
  void Update(double alpha, const DoubleVector& rhs, double beta)
  {
    int n = MyLength();
    if (n == 0) return;

    CheckSpaces(rhs);

    int incr = 1;
    // scale this by beta
    DSCAL_F77(&n, &beta, Values(), &incr);
    // computes this = alpha * rhs + this
    DAXPY_F77(&n, &alpha, (double*)rhs.Values(), &incr, Values(), &incr);
  }

  //! Computes the dot product between \c this vector and \c rhs.
  inline double DotProduct(const DoubleVector& rhs) const 
  {
    assert (rhs.VectorSpace() == VectorSpace());

    double MyResult = 0.0, Result = 0.0;
    int n = MyLength();
    int incr = 1;
    MyResult = DDOT_F77(&n, (double*)Values(), &incr, (double*)rhs.Values(), &incr);

    Result = ML_Comm_GsumDouble(GetML_Comm(),MyResult);
    return(Result);
  }

  //! Computes the 2-norm of \c this vector.
  inline double Norm2() const 
  {
    double MyResult = 0.0, Result = 0.0;
    int n = MyLength();
    int incr = 1;
    MyResult = DDOT_F77(&n, (double*)Values(), &incr, (double*)Values(), &incr);

    Result = ML_Comm_GsumDouble(GetML_Comm(),MyResult);
    return(sqrt(Result));
  }

  //! Computes the infinite norm of \c this vector.
  inline double NormInf() const 
  {
    double MyResult = 0.0, Result = 0.0;
    int n = MyLength();
    int incr = 1;
    int i = IDAMAX_F77(&n, (double*)Values(), &incr);
    MyResult = Values()[i - 1];

    Result = ML_Comm_GmaxDouble(GetML_Comm(),MyResult);
    return(Result);
  }

  //! Populates the vector with random elements.
  inline void Random() {
    ML_random_vec(Values(),MyLength(),MLAPI::GetML_Comm());
    return;
  }

  //! Replaces each element of the vector with its reciprocal.
  inline void Reciprocal() 
  {
    for (int i = 0 ; i < MyLength() ; ++i) {
      if (Values()[i] != 0.0)
        Values()[i] = 1.0 / Values()[i];
    }
  }

  //! Scales each element by the specified factor.
  inline void Scale(const double Factor) 
  {
    if (MyLength() == 0)
      return;

    int n = MyLength();
    int incr = 1;
    DSCAL_F77(&n, (double*)&Factor, Values(), &incr);
  }

  // @}
  // @{ Query methods
  
  //! Returns the Space on which \c this vector is defined.
  inline const Space& VectorSpace() const 
  {
    return(VectorSpace_);
  }

  //! Returns the length of \c this vector.
  inline int MyLength() const
  {
    return(VectorSpace_.NumMyElements());
  }

  //! Returns the length of \c this vector.
  inline int GlobalLength() const
  {
    return(VectorSpace_.NumGlobalElements());
  }

  //! Returns a pointer to the double array (non-const version)
  inline double* Values()
  {
    return(Values_.get());
  }

  //! Returns a pointer to the double array (const version)
  inline const double* Values() const
  {
    return(Values_.get());
  }

  //! Prints basic information about \c this object on ostream
  virtual std::ostream& Print(std::ostream& os,
                              const bool verbose = true) const
  {
    os << std::endl;
    if (MyPID() == 0) {
      os << "*** MLAPI::DoubleVector ***" << endl;
      os << "Label          = " << GetLabel() << endl;
      os << "MyLength()     = " << MyLength() << endl;
      os << "GlobalLength() = " << GlobalLength() << endl;
      os << endl << endl;
    }

    if (verbose) {

      if (MyPID() == 0) {
        os.width(10);
        os << "ProcID";
        os.width(20);
        os << "LID";
        os.width(20);
        os << "GID";
        os.width(20);
        os << "value" << endl << endl;
      }
      
      for (int iproc = 0 ; iproc < NumProc() ; ++iproc) {

        if (MyPID() == iproc) {

          for (int i = 0 ; i < MyLength() ; ++i) {
            os.width(10);
            os << MyPID();
            os.width(20);
            os << i;
            os.width(20);
            os << VectorSpace()(i);
            os.width(20);
            os << (*this)(i) << endl;
          }
        }

        GetEpetra_Comm().Barrier();
      }
    }

    return(os);
  }

private:

  // @}
  // @{ Internally used methods
  
  //! Initialize \c this object.
  inline void Initialize()
  {
    Values_ = Teuchos::null;
  }

  //! Verifies that \c rhs is compatible with \c this, and not its alias.
  void CheckSpaces(const DoubleVector rhs) 
  {
    if (rhs.VectorSpace() != VectorSpace()) {
      rhs.Print(cout,false);
      VectorSpace().Print(cout,false);
      ML_THROW("rhs.VectorSpace() is not equal to this->VectorSpace()", -1);
    }

    if (rhs.Values() == Values())
      ML_THROW("updating a vector with its alias...", -1);

  }

  //! Returns a pointer to the double array (non-const version)
  inline Teuchos::RefCountPtr<double> RCPValues() 
  {
    return(Values_);
  }

  //! Returns a pointer to the double array (const version)
  inline const Teuchos::RefCountPtr<double> RCPValues() const
  {
    return(Values_);
  }

  //! Sets the RefCountPtr<Values_>
  inline void SetRCPValues(const Teuchos::RefCountPtr<double> Values)
  {
    Values_ = Values;
  }
  // @}
  // @{ Internal data
  
  //! Pointer to locally own values.
  Teuchos::RefCountPtr<double> Values_;
  //! Data layout.
  Space VectorSpace_;

  // @}
  
}; // DoubleVector

} // namespace MLAPI

#endif // if ML_VECTOR_H
