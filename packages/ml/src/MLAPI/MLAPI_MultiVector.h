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
\class MultiVector

\brief Basic class for distributed double-precision vectors.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

class MultiVector : public BaseObject {

public:

  //@{ Constructors and destructors

  //! Default constructor.
  MultiVector() { }

  //! Constructor for a given Space.
  MultiVector(const Space& VectorSpace, const int NumVectors = 1,
               bool SetToZero = true)
  {
    NumVectors_  = NumVectors;
    VectorSpace_ = VectorSpace;
    if (GetMyTotalLength()) {
      SetRCPValues(Teuchos::rcp(new double[GetMyTotalLength()]));
      *this = 0.0;
    }
  }

  //! Constructor with a given Space, and user-provided array of values.
  MultiVector(const Space& VectorSpace, double* Values,
               const int NumVectors = 1)
  {
    NumVectors_  = NumVectors;
    VectorSpace_ = VectorSpace;
    SetRCPValues(Teuchos::rcp(Values, false));
  }

  //! Copy constructor.
  MultiVector(const MultiVector& rhs)
  {
    NumVectors_  = rhs.GetNumVectors();
    VectorSpace_ = rhs.GetVectorSpace();
    SetRCPValues(rhs.GetRCPValues());
  }

  //! Destructor.
  ~MultiVector() 
  {
    SetRCPValues(Teuchos::null);
  }

  // @}
  // @{ Overloaded operators

  //! Sets all elements of this vector to \c rhs.
  MultiVector& operator=(double rhs) 
  {
    for (int i = 0 ; i < GetMyTotalLength() ; ++i)
      GetValues()[i] = rhs;

    return(*this);
  }

  //! Copies the \c rhs into \c this object.
  MultiVector& operator=(const MultiVector& rhs) 
  {
    if (this != &rhs) {
      NumVectors_  = rhs.GetNumVectors();
      VectorSpace_ = rhs.GetVectorSpace();
      SetRCPValues(rhs.GetRCPValues());
      SetLabel(rhs.GetLabel());
    }

    return(*this);
  }

  //! Sets the name of \c this object, does not touch vector elements or space.
  MultiVector& operator=(const string& Label)
  {
    SetLabel(Label);
    return(*this);
  }

  //! Returns the value of local element \c i (const version).
  inline const double& operator() (const int i) const
  {
#ifdef MLAPI_CHECK
    if ((i < 0) || (i >= GetMyLength()))
      ML_THROW("Requested component " + toString(i) +
               ", while MyLength() = " + toString(GetMyLength()), -1);
#endif
    return(GetValues()[i]);
  }

  //! Returns the value of local element \c i (non-const version).
  inline double& operator() (const int i) 
  {
#ifdef MLAPI_CHECK
    if ((i < 0) || (i >= GetMyLength()))
      ML_THROW("Requested component " + toString(i) +
               ", while MyLength() = " + toString(GetMyLength()), -1);
#endif
    return(GetValues()[i]);
  }

  //! Returns the value of local element \c i.
  inline const double& operator()(const int i, const int v) const 
  {
#ifdef MLAPI_CHECK
    if (i < 0) || (i >= GetMyLength())
      ML_THROW("Requested component " + toString(i) +
               ", while MyLength() = " + toString(GetMyLength()), -1);
    if (v < 0) || (v >= GetNumVectors())
      ML_THROW("Requested vector " + toString(v) +
               ", while NumVectors() = " + toString(GetNumVectors()), -1);
#endif
    return(GetValues()[i + v * GetMyLength()]);
  }

  //! Returns the value of local element \c i (non-const version)
  inline double& operator()(const int i, const int v) 
  {
#ifdef MLAPI_CHECK
    if (i < 0) || (i >= GetMyLength())
      ML_THROW("Requested component " + toString(i) +
               ", while MyLength() = " + toString(GetMyLength()), -1);
    if (v < 0) || (v >= GetNumVectors())
      ML_THROW("Requested vector " + toString(v) +
               ", while NumVectors() = " + toString(GetNumVectors()), -1);
#endif
    return(GetValues()[i + v * GetMyLength()]);
  }

  // @}
  // @{ Mathematical methods
  
  //! Sets the space of this vector.
  void Reshape(const Space& S, const int NumVectors = 1)
  {
    NumVectors_ = NumVectors;
    VectorSpace_ = S;
    if (GetMyTotalLength())
      SetRCPValues(Teuchos::rcp(new double[GetMyTotalLength()]));
    else
      SetRCPValues(Teuchos::null);
  }

  //! Sets this = rhs.
  void Update(const MultiVector& rhs)
  {
    int n = GetMyTotalLength();
    if (n == 0) return;

    CheckSpaces(rhs);

    int incr = 1;
    // copy rhs into this
    DCOPY_F77(&n, (double*)rhs.GetValues(), &incr, GetValues(), &incr);
  }
  
  //! Sets this = alpha * rhs.
  void Update(double alpha, const MultiVector& rhs)
  {
    int n = GetMyTotalLength();
    if (n == 0) return;

    CheckSpaces(rhs);

    int incr = 1;
    // copy rhs into this
    DCOPY_F77(&n, (double*)rhs.GetValues(), &incr, GetValues(), &incr);
    // scale this
    DSCAL_F77(&n, &alpha, GetValues(), &incr);
  }

  //! Sets this = alpha * x + beta * y.
  void Update(double alpha, const MultiVector& x,
              double beta,  const MultiVector& y)
  {
    int n = GetMyTotalLength();
    if (n == 0) return;

    CheckSpaces(x);
    CheckSpaces(y);

    int incr = 1;
    // copy rhs into this
    DCOPY_F77(&n, (double*)x.GetValues(), &incr, GetValues(), &incr);
    Update(beta,y,alpha);

  }

  //! Sets this = alpha * rhs + beta * this.
  void Update(double alpha, const MultiVector& rhs, double beta)
  {
    int n = GetMyTotalLength();
    if (n == 0) return;

    CheckSpaces(rhs);

    int incr = 1;
    // scale this by beta
    DSCAL_F77(&n, &beta, GetValues(), &incr);
    // computes this = alpha * rhs + this
    DAXPY_F77(&n, &alpha, (double*)rhs.GetValues(), &incr, GetValues(), &incr);
  }

  //! Computes the dot product between \c this vector and \c rhs.
  inline double DotProduct(const MultiVector& rhs) const 
  {
    assert (rhs.GetVectorSpace() == GetVectorSpace());

    double MyResult = 0.0, Result = 0.0;
    int n = GetMyTotalLength();
    int incr = 1;
    MyResult = DDOT_F77(&n, (double*)GetValues(), &incr, (double*)rhs.GetValues(), &incr);

    Result = ML_Comm_GsumDouble(GetML_Comm(),MyResult);
    return(Result);
  }

  //! Computes the 2-norm of \c this vector.
  inline double Norm2() const 
  {
    double MyResult = 0.0, Result = 0.0;
    int n = GetMyTotalLength();
    int incr = 1;
    MyResult = DDOT_F77(&n, (double*)GetValues(), &incr, (double*)GetValues(), &incr);

    Result = ML_Comm_GsumDouble(GetML_Comm(),MyResult);
    return(sqrt(Result));
  }

  //! Computes the infinite norm of \c this vector.
  inline double NormInf() const 
  {
    double MyResult = 0.0, Result = 0.0;
    int n = GetMyTotalLength();
    int incr = 1;
    int i = IDAMAX_F77(&n, (double*)GetValues(), &incr);
    MyResult = GetValues()[i - 1];

    Result = ML_Comm_GmaxDouble(GetML_Comm(),MyResult);
    return(Result);
  }

  //! Populates the vector with random elements.
  inline void Random() 
  {
    ML_random_vec(GetValues(),GetMyTotalLength(),MLAPI::GetML_Comm());
    return;
  }

  //! Replaces each element of the vector with its reciprocal.
  inline void Reciprocal() 
  {
    for (int i = 0 ; i < GetMyTotalLength() ; ++i) {
      if (GetValues()[i] != 0.0)
        GetValues()[i] = 1.0 / GetValues()[i];
    }
  }

  //! Scales each element by the specified factor.
  inline void Scale(const double Factor) 
  {
    int n = GetMyTotalLength();
    if (n == 0) return;

    int incr = 1;
    DSCAL_F77(&n, (double*)&Factor, GetValues(), &incr);
  }

  // @}
  // @{ Query methods
  
  //! Returns the Space on which \c this vector is defined.
  inline const Space& GetVectorSpace() const 
  {
    return(VectorSpace_);
  }

  //! Returns the number of vectors.
  inline int GetNumVectors() const
  {
    return(NumVectors_);
  }

  //! Returns the local length of each vector.
  inline int GetMyLength() const
  {
    return(VectorSpace_.GetNumMyElements());
  }

  //! Returns the global length of each vector.
  inline int GetGlobalLength() const
  {
    return(VectorSpace_.GetNumGlobalElements());
  }

  //! Returns the local length of allocated vector, MyLength() * NumVectors()
  inline int GetMyTotalLength() const
  {
    return(VectorSpace_.GetNumMyElements() * GetNumVectors());
  }

  //! Returns the global length of allocated vector, GlobalLength() * NumVectors()
  inline int GetGlobalTotalLength() const
  {
    return(VectorSpace_.GetNumGlobalElements() * GetNumVectors());
  }

  //! Returns a pointer to the double array (non-const version)
  inline double* GetValues()
  {
    return(RCPValues_.get());
  }

  //! Returns a pointer to the double array (const version)
  inline const double* GetValues() const
  {
    return(RCPValues_.get());
  }

  //! Prints basic information about \c this object on ostream
  virtual std::ostream& Print(std::ostream& os,
                              const bool verbose = true) const
  {
    os << std::endl;
    if (GetMyPID() == 0) {
      os << "*** MLAPI::MultiVector ***" << endl;
      os << "Label             = " << GetLabel() << endl;
      os << "Local length      = " << GetMyLength() << endl;
      os << "Global length     = " << GetGlobalLength() << endl;
      os << "Number of vectors = " << GetNumVectors() << endl;
      os << endl << endl;
    }

    if (verbose) {

      if (GetMyPID() == 0) {
        os.width(10);
        os << "ProcID";
        os.width(20);
        os << "LID";
        os.width(20);
        os << "GID";
        for (int v = 0 ; v < GetNumVectors() ; ++v) {
          os.width(20);
          os << "value " << v;
        }
        os << endl << endl;
      }
      
      for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

        if (GetMyPID() == iproc) {

          for (int i = 0 ; i < GetMyLength() ; ++i) {
            os.width(10);
            os << GetMyPID();
            os.width(20);
            os << i;
            os.width(20);
            os << GetVectorSpace()(i);
            for (int v = 0 ; v < GetNumVectors() ; ++v) {
              os.width(20);
              os << (*this)(i,v);
            }
            os << endl;
          }
        }

        Barrier();
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
    RCPValues_ = Teuchos::null;
  }

  //! Verifies that \c rhs is compatible with \c this, and not its alias.
  void CheckSpaces(const MultiVector rhs) 
  {
    if (rhs.GetVectorSpace() != GetVectorSpace()) {
      ML_THROW("rhs.GetVectorSpace() is not equal to this->GetVectorSpace()", -1);
    }

    if (rhs.GetValues() == GetValues())
      ML_THROW("updating a vector with its alias...", -1);

    if (rhs.GetNumVectors() != GetNumVectors())
      ML_THROW("rhs and this have different number of vectors" +
               toString(rhs.GetNumVectors()) + " vs. " +
               toString(GetNumVectors()) + ")", -1);

  }

  //! Returns a pointer to the double array (non-const version)
  inline Teuchos::RefCountPtr<double> GetRCPValues() 
  {
    return(RCPValues_);
  }

  //! Returns a pointer to the double array (const version)
  inline const Teuchos::RefCountPtr<double> GetRCPValues() const
  {
    return(RCPValues_);
  }

  //! Sets the RefCountPtr<Values_>
  inline void SetRCPValues(const Teuchos::RefCountPtr<double> RCPValues)
  {
    RCPValues_ = RCPValues;
  }

  // @}
  // @{ Internal data
  
  //! Pointer to locally own values.
  Teuchos::RefCountPtr<double> RCPValues_;
  //! Data layout.
  Space VectorSpace_;
  //! Number of vectors.
  int NumVectors_;

  // @}
  
}; // MultiVector

} // namespace MLAPI

#endif // if ML_VECTOR_H
