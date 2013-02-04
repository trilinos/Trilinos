#ifndef ML_MULTIVECTOR_H
#define ML_MULTIVECTOR_H

/*!
\file MLAPI_MultiVector.h

\brief MLAPI wrapper for double vectors.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

//#include "ml_lapack.h"
#include "MLAPI_Error.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_BaseLinearCombination.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_BLAS_wrappers.hpp"
#include "Teuchos_LAPACK_wrappers.hpp"
#include <iomanip>
#include <vector>
#include <algorithm>
#include <functional>
#include "ml_utils.h"

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

class MultiVector;

    
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
    StackPush();

    NumVectors_  = NumVectors;
    VectorSpace_ = VectorSpace;
    SetRCPLength(NumVectors);
    //if (GetMyLength()) {
      for (int v = 0 ; v < NumVectors ; ++v)
        SetRCPValues(Teuchos::rcp(new DoubleVector(GetMyLength())), v);

      if (SetToZero)
        *this = 0.0;
    //}

    StackPop();
  }

  //! Constructor with a given Space, and user-provided array of values.
  MultiVector(const Space& VectorSpace, double** Values,
              const int NumVectors = 1)
  {
    StackPush();

    NumVectors_  = NumVectors;
    VectorSpace_ = VectorSpace;
    SetRCPLength(GetNumVectors());
    for (int v = 0 ; v < GetNumVectors() ; ++v)
      SetRCPValues(Teuchos::rcp(new DoubleVector(Values[v])), v);

    StackPop();
  }

  //! Constructor with a given Space, and user-provided RefCountPtr
  MultiVector(const Space& VectorSpace, Teuchos::RefCountPtr<DoubleVector> RCPValues)
  {
    StackPush();

    NumVectors_  = 1;
    VectorSpace_ = VectorSpace;
    SetRCPLength(GetNumVectors());
    SetRCPValues(RCPValues, 0);

    StackPop();
  }

  //! Constructor with a given Space, and user-provided array of values.
  MultiVector(const Space& VectorSpace, 
              std::vector<Teuchos::RefCountPtr<DoubleVector> > RCPValues)
  {
    StackPush();

    NumVectors_  = (int)RCPValues.size();
    VectorSpace_ = VectorSpace;
    SetRCPLength(GetNumVectors());
    for (int v = 0 ; v < GetNumVectors() ; ++v)
      SetRCPValues(RCPValues[v], v);

    StackPop();
  }

  //! Copy constructor.
  MultiVector(const MultiVector& rhs)
  {
    NumVectors_  = rhs.GetNumVectors();
    VectorSpace_ = rhs.GetVectorSpace();
    SetRCPLength(GetNumVectors());
    for (int v = 0 ; v < GetNumVectors() ; ++v)
      SetRCPValues(rhs.GetRCPValues(v), v);
  }

  //! Destructor.
  ~MultiVector() 
  {
    for (int v = 0 ; v < GetNumVectors() ; ++v)
      SetRCPValues(Teuchos::null, v);
  }

  // @}
  // @{ \name Reshape methods
  
  //! Resets \c this object.
  void Reshape()
  {
    StackPush();

    for (int v = 0 ; v < GetNumVectors() ; ++v)
      SetRCPValues(Teuchos::null, v);
    SetRCPLength(0);
    GetVectorSpace().Reshape();
    NumVectors_ = 0;

    StackPop();
  }

  //! Sets the space of this vector.
  void Reshape(const Space& S, const int NumVectors = 1,
               const bool SetToZero = true)
  {
    StackPush();

    NumVectors_ = NumVectors;
    VectorSpace_ = S;
    SetRCPLength(GetNumVectors());
    for (int v = 0 ; v < GetNumVectors() ; ++v) {
      //if (GetMyLength())
        SetRCPValues(Teuchos::rcp(new DoubleVector(GetMyLength())), v);
      //else
      //  SetRCPValues(Teuchos::null, v);
    }

    if (SetToZero) *this = 0.0;

    StackPop();
  }

  //! Appends a new vector.
  void Append(const int NumVectors = 1, const bool SetToZero = true)
  {
    int n = GetMyLength();

    for (int v = 0 ; v < NumVectors ; ++v) {
      //if (GetMyLength())
        RCPValues_.push_back(Teuchos::rcp(new DoubleVector(n)));
      //else
      //  RCPValues_.push_back(Teuchos::null);

      ++NumVectors_;

      if (SetToZero) {
        Update(0.0, GetNumVectors() - 1);
      }
    }
  }
   
  //! Appends a new vector.
  void Append(MultiVector rhs)
  {
    StackPush();

    CheckSpaces(rhs);

    for (int v = 0 ; v < rhs.GetNumVectors() ; ++v) {
      RCPValues_.push_back(rhs.GetRCPValues(v));
      ++NumVectors_;
    }

    StackPop();
  }

  //! Deletes the last vector.
  void Delete(const int v)
  {
    StackPush();

    CheckVector(v);

    std::vector<Teuchos::RefCountPtr<DoubleVector> > NewList;

    for (int i = 0 ; i < GetNumVectors() ; ++i) {
      if (i != v)
        NewList.push_back(GetRCPValues(i));
    }

    RCPValues_ = NewList;

    NumVectors_--;

    StackPop();
  }

  // @}
  // @{ \name Overloaded operators

  //! Sets all elements of this vector to \c rhs.
  MultiVector& operator=(double rhs) 
  {
    StackPush();

    for (int v = 0 ; v < GetNumVectors() ; ++v)
      for (int i = 0 ; i < GetMyLength() ; ++i)
        GetValues(v)[i] = rhs;

    StackPop();

    return(*this);
  }

  //! Copies the \c rhs into \c this object.
  MultiVector& operator=(const MultiVector& rhs) 
  {
    StackPush();

    if (this != &rhs) {
      NumVectors_  = rhs.GetNumVectors();
      VectorSpace_ = rhs.GetVectorSpace();
      SetRCPLength(GetNumVectors());
      for (int v = 0 ; v < GetNumVectors() ; ++v)
        SetRCPValues(rhs.GetRCPValues(v), v);
      SetLabel(rhs.GetLabel());
    }

    StackPop();

    return(*this);
  }

  //! Sets the elements from the input BaseLinearCombination
  MultiVector& operator=(const BaseLinearCombination& rhs)
  {
    StackPush();

    rhs.Set(*this);

    StackPop();

    return(*this);
  }

  bool operator==(const MultiVector& rhs) const
  {
    return(false);
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
    CheckSingleVector();
    CheckEntry(i);
    CheckVector(0);

    return(GetValues(0)[i]);
  }

  //! Returns the value of local element \c i (non-const version).
  inline double& operator() (const int i) 
  {
    CheckSingleVector();
    CheckEntry(i);
    CheckVector(0);

    return(GetValues(0)[i]);
  }

  //! Returns the value of local element \c i.
  inline const double& operator()(const int i, const int v) const 
  {
    CheckEntry(i);
    CheckVector(v);

    return(GetValues(v)[i]);
  }

  //! Returns the value of local element \c i (non-const version)
  inline double& operator()(const int i, const int v) 
  {
    CheckEntry(i);
    CheckVector(v);

    return(GetValues(v)[i]);
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

  //! Returns a pointer to the double array (non-const version)
  inline double* GetValues(const int v)
  {
    return(RCPValues_[v].get()->Values());
  }

  //! Returns a pointer to the double array (const version)
  inline const double* GetValues(const int v) const
  {
    return(RCPValues_[v].get()->Values());
  }
  
  //! Returns a pointer to the double array (non-const version)
  inline Teuchos::RefCountPtr<DoubleVector>& GetRCPValues(const int v) 
  {
    CheckVector(v);

    return(RCPValues_[v]);
  }

  //! Returns a pointer to the double array (const version)
  inline const Teuchos::RefCountPtr<DoubleVector>& GetRCPValues(const int v) const
  {
    CheckVector(v);

    return(RCPValues_[v]);
  }

  //! Sets the RefCountPtr<Values_>
  inline void SetRCPValues(const Teuchos::RefCountPtr<DoubleVector>& RCPValues,
                           const int v)
  {
    CheckVector(v);

    RCPValues_[v] = RCPValues;
  }

  // @}
  // @{ \name Mathematical methods
  
  //! Sets this(v) = rhs
  void Update(const double alpha, int v = -1)
  {
    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    if (v >= GetNumVectors())
      ML_THROW("Requested wrong vector, " + GetString(v) +
               " while NumVectors = " + GetString(GetNumVectors()), -1);

    for (int i = 0 ; i < GetMyLength() ; ++i)
      GetValues(v)[i] = alpha;
  }

  //! Sets this = rhs.
  void Update(const MultiVector& rhs)
  {
    ResetTimer();
    StackPush();

    int n = GetMyLength();
    if (n == 0) return;

    CheckSpaces(rhs);
    CheckNumVectors(rhs.GetNumVectors());

    int incr = 1;
    for (int v = 0 ; v < GetNumVectors() ; ++v)
      // copy rhs into this
      DCOPY_F77(&n, (double*)rhs.GetValues(v), &incr, GetValues(v), &incr);

    StackPop();
    UpdateTime();
  }
  
  //! Sets this = alpha * rhs.
  void Update(double alpha, const MultiVector& rhs)
  {
    ResetTimer();
    StackPush();

    int n = GetMyLength();
    if (n == 0) return;

    CheckSpaces(rhs);
    CheckNumVectors(rhs.GetNumVectors());

    int incr = 1;
    for (int v = 0 ; v < GetNumVectors() ; ++v) {
      // copy rhs into this
      DCOPY_F77(&n, (double*)rhs.GetValues(v), &incr, GetValues(v), &incr);
      // scale this
      DSCAL_F77(&n, &alpha, GetValues(v), &incr);
    }

    StackPop();
    UpdateFlops(1.0 * GetNumVectors() * GetGlobalLength());
    UpdateTime();

  }

  //! Sets this = alpha * x + beta * y.
  void Update(double alpha, const MultiVector& x,
              double beta,  const MultiVector& y)
  {
    ResetTimer();
    StackPush();

    int n = GetMyLength();
    if (n == 0) return;

    CheckSpaces(x);
    CheckSpaces(y);
    CheckNumVectors(x.GetNumVectors());
    CheckNumVectors(y.GetNumVectors());

    int incr = 1;
    for (int v = 0 ; v < GetNumVectors() ; ++v)
      // copy rhs into this
      DCOPY_F77(&n, (double*)x.GetValues(v), &incr, GetValues(v), &incr);

    StackPop();
    Update(beta,y,alpha);
    UpdateTime();

  }

  //! Sets this = alpha * rhs + beta * this.
  void Update(double alpha, const MultiVector& rhs, double beta)
  {
    ResetTimer();
    StackPush();

    int n = GetMyLength();
    if (n == 0) return;

    CheckSpaces(rhs);
    CheckNumVectors(rhs.GetNumVectors());

    for (int v = 0 ; v < GetNumVectors() ; ++v) 
    {
      double* ptr_this = GetValues(v);
      double* ptr_rhs  = (double*)(rhs.GetValues(v));

      if (alpha == 1.0 && beta == 1.0)
      {
        for (int i = 0 ; i < GetMyLength() ; ++i)
          ptr_this[i] += ptr_rhs[i];
        UpdateFlops(GetGlobalLength());
      }
      else if (alpha == 1.0 && beta == 0.0)
      {
        for (int i = 0 ; i < GetMyLength() ; ++i)
          ptr_this[i] = ptr_rhs[i];
      }
      else if (alpha == 0.0 && beta == 1.0)
      {
        // do nothing here
        if (false)
          cout << "blablablaaaaa" << endl;
      }
      else if (alpha == 1.0 && beta == -1.0)
      {
        for (int i = 0 ; i < GetMyLength() ; ++i)
          ptr_this[i] = ptr_rhs[i] - ptr_this[i];
        UpdateFlops(GetGlobalLength()); 
      }
      else if (alpha == -1.0 && beta == 1.0)
      {
        for (int i = 0 ; i < GetMyLength() ; ++i)
          ptr_this[i] -= ptr_rhs[i];
        UpdateFlops(GetGlobalLength()); 
      }
      else
      {
        for (int i = 0 ; i < GetMyLength() ; ++i)
          ptr_this[i] = ptr_rhs[i] * alpha + ptr_this[i] * beta;
        UpdateFlops(3.0 * GetGlobalLength()); 
      }
    }

    StackPop();
    UpdateTime();
  }

#if 0
  void Add(double alpha)
  {
    for (int v = 0 ; v < GetNumVectors() ; ++v)
      for (int i = 0 ; i < GetMyLength() ; ++i)
        GetValues(v)[i] += alpha;
  }
#endif

  //! Computes the dot product between \c this vector and \c rhs.
  inline double DotProduct(const MultiVector& rhs, int v = -1) const 
  {
    ResetTimer();
    StackPush();

    if (rhs.GetVectorSpace() != GetVectorSpace()) {
      ML_THROW("rhs.GetVectorSpace() is not equal to this->GetVectorSpace()", -1);
    }

    CheckNumVectors(rhs.GetNumVectors());
    
    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    double MyResult = 0.0;
    double Result   = 0.0;
    int n           = GetMyLength();
    int incr        = 1;
    double* ptr     = (double*)GetValues(v);
    double* rhs_ptr = (double*)rhs.GetValues(v);
    MyResult        = DDOT_F77(&n, ptr, &incr, rhs_ptr, &incr);
    Result          = ML_Comm_GsumDouble(GetML_Comm(),MyResult);

    StackPop();
    UpdateFlops(2.0 * GetGlobalLength()); // DDOT
    UpdateTime();

    return(Result);
  }

  //! Computes the 2-norm of \c this vector.
  inline double Norm2(int v = -1) const 
  {
    ResetTimer();
    StackPush();

    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    double MyResult = 0.0;
    double Result   = 0.0;
    int n           = GetMyLength();
    int incr        = 1;
    double* ptr     = (double*)GetValues(v);
    MyResult        = DDOT_F77(&n, ptr, &incr, ptr, &incr);
    Result          = ML_Comm_GsumDouble(GetML_Comm(),MyResult);
    
    StackPop();
    UpdateFlops(2.0 * GetGlobalLength()); // DDOT
    UpdateTime();

    return(sqrt(Result));
  }

  //! Computes the infinite norm of \c this vector.
  inline double NormInf(int v = -1) const 
  {
    ResetTimer();
    StackPush();

    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    double MyResult = 0.0;
    double Result   = 0.0;
    int n           = GetMyLength();
    double* ptr     = (double*)GetValues(v);
    int incr        = 1;
    int i           = IDAMAX_F77(&n, ptr, &incr);
    MyResult        = fabs(ptr[i - 1]);
    // FIXME: delete below
    /*
    for (int i = 0 ; i < n ; ++i)
      if (MyResult < fabs(ptr[i])) 
        MyResult = fabs(ptr[i]);
        */

    Result          = ML_Comm_GmaxDouble(GetML_Comm(),MyResult);

    StackPop();
    UpdateTime();
    return(Result);
  }

  //! Computes the one norm of \c this vector.
  inline double NormOne(int v = -1) const 
  {
    ResetTimer();
    StackPush();

    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    double MyResult = 0.0;
    double Result   = 0.0;
    double* ptr     = (double*)GetValues(v);
    for (int i = 0 ; i < GetMyLength() ; ++i)
      MyResult += fabs(ptr[i]);

    Result          = ML_Comm_GsumDouble(GetML_Comm(),MyResult);

    StackPop();
    UpdateTime();
    return(Result);
  }

  //! Replaces each element of the vector with its reciprocal.
  inline void Reciprocal(int v = -1) 
  {
    ResetTimer();
    StackPush();
    
    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    for (int i = 0 ; i < GetMyLength() ; ++i) {
      if (GetValues(v)[i] != 0.0)
        GetValues(v)[i] = 1.0 / GetValues(v)[i];
    }

    StackPop();

    UpdateFlops(1.0 * GetGlobalLength());
    UpdateTime();
  }

  //! Scales each element by the specified factor.
  inline void Scale(const double Factor, int v = -1) 
  {
    ResetTimer();
    StackPush();

    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    int n = GetMyLength();
    if (n == 0) return;

    int incr = 1;
    DSCAL_F77(&n, (double*)&Factor, GetValues(v), &incr);

    StackPop();

    UpdateFlops(1.0 * GetGlobalLength()); 
    UpdateTime();
  }

  // @}
  // @{ \name Miscellanous methods

  //! Populates the vector with random elements.
  inline void Random(int v = -1) 
  {
    ResetTimer();
    StackPush();

    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    ML_random_vec(GetValues(v),GetMyLength(),MLAPI::GetML_Comm());

    StackPop();
    UpdateTime();
    return;
  }

  //! Sorts the component of the vector.
  void Sort(int v = -1, const bool IsIncreasing = false)
  {
    ResetTimer();
    StackPush();

    if (v == -1) {
      CheckSingleVector();
      v = 0;
    }

    CheckVector(v);

    std::vector<double> tmp(GetMyLength());
    for (int i = 0 ; i < GetMyLength() ; ++i)
      tmp[i] = (*this)(i, v);

    if (IsIncreasing)
      std::sort(tmp.begin(), tmp.end(), std::greater<double>());
    else
      std::sort(tmp.begin(), tmp.end());

    for (int i = 0 ; i < GetMyLength() ; ++i)
      (*this)(i,v) = tmp[i];

    StackPop();
    UpdateTime();
  }

  //! Prints basic information about \c this object on ostream
  virtual std::ostream& Print(std::ostream& os,
                              const bool verbose = true) const
  {
    ResetTimer();
    StackPush();

    if (GetMyPID() == 0) {
      os << endl;
      os << "*** MLAPI::MultiVector ***" << endl;
      os << "Label             = " << GetLabel() << endl;
      os << "Local length      = " << GetMyLength() << endl;
      os << "Global length     = " << GetGlobalLength() << endl;
      os << "Number of vectors = " << GetNumVectors() << endl;
      os << "Flop count        = " << GetFlops() << endl;
      os << "Cumulative time   = " << GetTime() << endl;
      if (GetTime() != 0.0)
        os << "MFlops rate       = " << 1.0e-6 * GetFlops() / GetTime() << endl;
      else
        os << "MFlops rate       = 0.0" << endl;
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
      if (GetMyPID() == 0)
        os << endl;
    }

    StackPop();
    UpdateTime();

    return(os);
  }
  // @}

  friend std::ostream& operator<<(std::ostream& os, MultiVector const &mv)
  {
    mv.Print(os);
    return os;
  }

  bool IsAlias(const MultiVector& rhs) const
  {
    if (rhs.GetValues(0) == GetValues(0))
      return(true);
    else
      return(false);
  }

private:

  //! Sets the length of RCPValues_ array
  void SetRCPLength(const int NumVectors)
  {
    RCPValues_.resize(NumVectors);
  }

  //! Initialize \c this object.
  inline void Initialize()
  {
    for (int v = 0 ; v < GetNumVectors() ; ++v)
      RCPValues_[v] = Teuchos::null;
  }

  //! Verifies that \c rhs is compatible with \c this, and not its alias.
  void CheckSpaces(const MultiVector rhs)  const
  {
    if (rhs.GetVectorSpace() != GetVectorSpace()) {
      ML_THROW("rhs.GetVectorSpace() is not equal to this->GetVectorSpace()", -1);
    }

    if (IsAlias(rhs))
      ML_THROW("updating a vector with its alias...", -1);
  }

  //! Verifies that the requested component is actually stored
  inline void CheckEntry(const int i) const
  {
#ifdef MLAPI_CHECK
    if ((i < 0) || (i >= GetMyLength()))
      ML_THROW("Requested component " + GetString(i) +
               ", while MyLength() = " + GetString(GetMyLength()), -1);
#endif
  }

  //! Verifies that the requested component is actually stored
  inline void CheckVector(const int v) const
  {
#ifdef MLAPI_CHECK
    if (v < 0 || v >= GetNumVectors())
      ML_THROW("Requested vector " + GetString(v) +
               ", while NumVectors() = " + GetString(GetNumVectors()), -1);
#endif
  }

  //! Verifies the number of vectors.
  inline void CheckNumVectors(const int NumVectors) const
  {
    if (GetNumVectors() != NumVectors)
      ML_THROW("Incompatible number of vectors, " +
               GetString(GetNumVectors()) + " vs. " +
               GetString(NumVectors), -1);
  }

  //! Verifies that only one vector is stored
  inline void CheckSingleVector() const
  {
    if (GetNumVectors() != 1)
      ML_THROW("Implicitly requested vector 0, while NumVectors = "
               + GetString(GetNumVectors()), -1);
  }

  //! Pointer to locally own values.
  std::vector<Teuchos::RefCountPtr<DoubleVector> > RCPValues_;
  //! Data layout.
  Space VectorSpace_;
  //! Number of vectors.
  int NumVectors_;

}; // MultiVector

} // namespace MLAPI

#endif // if ML_MULTIVECTOR_H
