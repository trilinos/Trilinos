#ifndef ML_VECTOR_H
#define ML_VECTOR_H

#include "ml_lapack.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "Teuchos_RefCountPtr.hpp"
#include <iomanip>
#include <vector>
#include <algorithm>

namespace MLAPI {

class DoubleVector {

public:
  DoubleVector(const int size)
  {
    ownership_ = true;
    ptr_ = new double[size];
  }

  DoubleVector(double* ptr)
  {
    ownership_ = false;
    ptr_ = ptr;
  }

  inline double& operator[](const int i)
  {
    return(ptr_[i]);
  }

  inline const double& operator[](const int i) const
  {
    return(ptr_[i]);
  }

  ~DoubleVector()
  {
    if (ownership_)
      delete[] ptr_;
  }

  double* Values()
  {
    return(ptr_);
  }

  const double* Values() const
  {
    return(ptr_);
  }

private:
  DoubleVector(const DoubleVector& rhs)
  {}

  DoubleVector& operator=(const DoubleVector& rhs)
  {
    return(*this);
  }
  
  double* ptr_;
  bool    ownership_;
};

/*!
\class MultiVector

\brief Basic class for distributed double-precision vectors.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

class MultiVector : public BaseObject, public CompObject, public TimeObject {

public:

  //@{ \name Constructors and destructors

  //! Default constructor.
  MultiVector() 
  { 
    NumVectors_ = 0;
  }

  //! Constructor for a given Space.
  MultiVector(const Space& VectorSpace, const int NumVectors = 1,
               bool SetToZero = true)
  {
    NumVectors_  = NumVectors;
    VectorSpace_ = VectorSpace;
    if (GetMyTotalLength()) {
      SetRCPValues(Teuchos::rcp(new DoubleVector(GetMyTotalLength())));
      *this = 0.0;
    }
  }

  //! Constructor with a given Space, and user-provided array of values.
  MultiVector(const Space& VectorSpace, double* Values,
              const int NumVectors = 1)
  {
    NumVectors_  = NumVectors;
    VectorSpace_ = VectorSpace;
    SetRCPValues(Teuchos::rcp(new DoubleVector(Values)));
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
  // @{ \name Reshape methods
  
  //! Resets \c this object.
  void Reshape()
  {
    SetRCPValues(Teuchos::null);
    GetVectorSpace().Reshape();
    NumVectors_ = 0;
  }

  //! Sets the space of this vector.
  void Reshape(const Space& S, const int NumVectors = 1)
  {
    NumVectors_ = NumVectors;
    VectorSpace_ = S;
    if (GetMyTotalLength())
      SetRCPValues(Teuchos::rcp(new DoubleVector(GetMyTotalLength())));
    else
      SetRCPValues(Teuchos::null);
  }

  // @}
  // @{ \name Overloaded operators

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
      ML_THROW("Requested component " + GetString(i) +
               ", while MyLength() = " + GetString(GetMyLength()), -1);
#endif
    return(GetValues()[i]);
  }

  //! Returns the value of local element \c i (non-const version).
  inline double& operator() (const int i) 
  {
#ifdef MLAPI_CHECK
    if ((i < 0) || (i >= GetMyLength()))
      ML_THROW("Requested component " + GetString(i) +
               ", while MyLength() = " + GetString(GetMyLength()), -1);
#endif
    return(GetValues()[i]);
  }

  //! Returns the value of local element \c i.
  inline const double& operator()(const int i, const int v) const 
  {
#ifdef MLAPI_CHECK
    if (i < 0) || (i >= GetMyLength())
      ML_THROW("Requested component " + GetString(i) +
               ", while MyLength() = " + GetString(GetMyLength()), -1);
    if (v < 0) || (v >= GetNumVectors())
      ML_THROW("Requested vector " + GetString(v) +
               ", while NumVectors() = " + GetString(GetNumVectors()), -1);
#endif
    return(GetValues()[i + v * GetMyLength()]);
  }

  //! Returns the value of local element \c i (non-const version)
  inline double& operator()(const int i, const int v) 
  {
#ifdef MLAPI_CHECK
    if (i < 0) || (i >= GetMyLength())
      ML_THROW("Requested component " + GetString(i) +
               ", while MyLength() = " + GetString(GetMyLength()), -1);
    if (v < 0) || (v >= GetNumVectors())
      ML_THROW("Requested vector " + GetString(v) +
               ", while NumVectors() = " + GetString(GetNumVectors()), -1);
#endif
    return(GetValues()[i + v * GetMyLength()]);
  }

  // @}
  // @{ \name Set and Get methods
  
  //! Returns the Space on which \c this vector is defined.
  inline const Space& GetVectorSpace() const 
  {
    return(VectorSpace_);
  }

  //! Returns the Space on which \c this vector is defined (non-const)
  inline Space& GetVectorSpace() 
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
    return(RCPValues_.get()->Values());
  }

  //! Returns a pointer to the double array (const version)
  inline const double* GetValues() const
  {
    return(RCPValues_.get()->Values());
  }
  
  //! Returns a pointer to the double array (non-const version)
  inline Teuchos::RefCountPtr<DoubleVector>& GetRCPValues() 
  {
    return(RCPValues_);
  }

  //! Returns a pointer to the double array (const version)
  inline const Teuchos::RefCountPtr<DoubleVector>& GetRCPValues() const
  {
    return(RCPValues_);
  }

  //! Sets the RefCountPtr<Values_>
  inline void SetRCPValues(const Teuchos::RefCountPtr<DoubleVector>& RCPValues)
  {
    RCPValues_ = RCPValues;
  }

  // @}
  // @{ \name Mathematical methods
  
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
    ResetTimer();

    int n = GetMyTotalLength();
    if (n == 0) return;

    CheckSpaces(rhs);

    int incr = 1;
    // copy rhs into this
    DCOPY_F77(&n, (double*)rhs.GetValues(), &incr, GetValues(), &incr);
    // scale this
    DSCAL_F77(&n, &alpha, GetValues(), &incr);

    UpdateFlops(1.0 * GetGlobalTotalLength());
    UpdateTime();

  }

  //! Sets this = alpha * x + beta * y.
  void Update(double alpha, const MultiVector& x,
              double beta,  const MultiVector& y)
  {
    ResetTimer();

    int n = GetMyTotalLength();
    if (n == 0) return;

    CheckSpaces(x);
    CheckSpaces(y);

    int incr = 1;
    // copy rhs into this
    DCOPY_F77(&n, (double*)x.GetValues(), &incr, GetValues(), &incr);

    Update(beta,y,alpha);
    UpdateTime();

  }

  //! Sets this = alpha * rhs + beta * this.
  void Update(double alpha, const MultiVector& rhs, double beta)
  {
    ResetTimer();

    int n = GetMyTotalLength();
    if (n == 0) return;

    CheckSpaces(rhs);

    int incr = 1;
    // scale this by beta
    DSCAL_F77(&n, &beta, GetValues(), &incr);
    // computes this = alpha * rhs + this
    DAXPY_F77(&n, &alpha, (double*)rhs.GetValues(), &incr, GetValues(), &incr);

    UpdateFlops(1.0 * GetGlobalTotalLength());     // DSCAL
    UpdateFlops(2.0 * GetGlobalTotalLength()); // DAXPY
    UpdateTime();
  }

  //! Computes the dot product between \c this vector and \c rhs.
  inline double DotProduct(const MultiVector& rhs, const int vector = 0) const 
  {
    assert (rhs.GetVectorSpace() == GetVectorSpace());
    assert (rhs.GetNumVectors() == GetNumVectors());

    ResetTimer();

    double MyResult = 0.0;
    double Result   = 0.0;
    int n           = GetMyLength();
    int incr        = 1;
    double* ptr     = (double*)GetValues() + vector * n;
    double* rhs_ptr = (double*)rhs.GetValues() + vector * n;
    MyResult        = DDOT_F77(&n, ptr, &incr, rhs_ptr, &incr);
    Result          = ML_Comm_GsumDouble(GetML_Comm(),MyResult);

    UpdateFlops(2.0 * GetGlobalTotalLength()); // DDOT
    UpdateTime();

    return(Result);
  }

  //! Computes the 2-norm of \c this vector.
  inline double Norm2(const int vector = 0) const 
  {
    ResetTimer();

    double MyResult = 0.0;
    double Result   = 0.0;
    int n           = GetMyLength();
    int incr        = 1;
    double* ptr     = (double*)GetValues() + vector * n;
    MyResult        = DDOT_F77(&n, ptr, &incr, ptr, &incr);
    Result          = ML_Comm_GsumDouble(GetML_Comm(),MyResult);
    
    UpdateFlops(2.0 * GetGlobalTotalLength()); // DDOT
    ResetTimer();


    return(sqrt(Result));
  }

  //! Computes the infinite norm of \c this vector.
  inline double NormInf(const int vector = 0) const 
  {
    double MyResult = 0.0;
    double Result   = 0.0;
    int n           = GetMyLength();
    int incr        = 1;
    double* ptr     = (double*)GetValues() + vector * n;
    int i           = IDAMAX_F77(&n, ptr, &incr);
    MyResult        = ptr[i - 1];
    Result          = ML_Comm_GmaxDouble(GetML_Comm(),MyResult);

    return(Result);
  }

  //! Replaces each element of the vector with its reciprocal.
  inline void Reciprocal() 
  {
    ResetTimer();

    for (int i = 0 ; i < GetMyTotalLength() ; ++i) {
      if (GetValues()[i] != 0.0)
        GetValues()[i] = 1.0 / GetValues()[i];
    }

    UpdateFlops(1.0 * GetGlobalTotalLength());
    UpdateTime();
  }

  //! Scales each element by the specified factor.
  inline void Scale(const double Factor) 
  {
    ResetTimer();

    int n = GetMyTotalLength();
    if (n == 0) return;

    int incr = 1;
    DSCAL_F77(&n, (double*)&Factor, GetValues(), &incr);

    UpdateFlops(1.0 * GetGlobalTotalLength()); 
    UpdateTime();
  }

  // @}
  // @{ \name Miscellanous methods

  //! Populates the vector with random elements.
  inline void Random() 
  {
    ML_random_vec(GetValues(),GetMyTotalLength(),MLAPI::GetML_Comm());
    return;
  }

  //! Sorts the component of the vector.
  void Sort(const int v = 0, const bool IsIncreasing = false)
  {
#ifdef MLAPI_CHECK
    if (v < 0) || (v >= GetNumVectors())
      ML_THROW("Requested vector " + GetString(v) +
               ", while NumVectors() = " + GetString(GetNumVectors()), -1);
#endif

    vector<double> tmp(GetMyLength());
    for (int i = 0 ; i < GetMyLength() ; ++i)
      tmp[i] = (*this)(i, v);

    if (IsIncreasing)
      sort(tmp.begin(), tmp.end(), greater<double>());
    else
      sort(tmp.begin(), tmp.end());

    for (int i = 0 ; i < GetMyLength() ; ++i)
      (*this)(i,v) = tmp[i];

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
      os << "Flop count        = " << GetFlops() << endl;
      os << "Cumulative time   = " << GetTime() << endl;
      os << "MFlops rate       = " << 1.0e-6 * GetFlops() / GetTime() << endl;
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
  // @}

private:

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
               GetString(rhs.GetNumVectors()) + " vs. " +
               GetString(GetNumVectors()) + ")", -1);

  }

  //! Pointer to locally own values.
  Teuchos::RefCountPtr<DoubleVector> RCPValues_;
  //! Data layout.
  Space VectorSpace_;
  //! Number of vectors.
  int NumVectors_;

}; // MultiVector

} // namespace MLAPI

#endif // if ML_VECTOR_H
