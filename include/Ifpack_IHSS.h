/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_IHSS_H
#define IFPACK_IHSS_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"

#ifdef HAVE_IFPACK_EPETRAEXT
#include "Ifpack_Condest.h"
#include "Ifpack_ScalingType.h"
#include "Epetra_CompObject.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"


namespace Teuchos {
  class ParameterList;
}

//! Ifpack_IHSS: A class for constructing and using the Inexact Hermetian Skew-Hermetian splitting method.

/*! Ifpack_IHSS: A class for constructing and using the Inexact Hermetian Skew-Hermetian splitting method, based on Bai, Golub and Ng (2003).  After constructing the appropriate matrices, Ifpack_IHSS then builds Ifpack preconditioners for the resulting matrices.

    \author Chris Siefert, SNL 1431.

    \date Last modified on June 20th, 2010.
*/

class Ifpack_IHSS: public Ifpack_Preconditioner {

public:
  // @{ Constructors and destructors.
  //! Constructor
  /*! Though this says Epetra_RowMatrix it must be an Epetra_CrsMatrix.  Errrors will ensue otherwise.
   */
  Ifpack_IHSS(Epetra_RowMatrix* A);

  //! Destructor
  ~Ifpack_IHSS()
  {
    Destroy();
  }

  // @}
  // @{ Construction methods

  //! Initialize the preconditioner, does not touch matrix values.
  int Initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Computes the Hermetian / Skew-Hermetian matrices & builds their preconditioners.
  /*! This function computes the Hermetian / Skew-Hermetian matrices & builds their preconditioners.
   */
  int Compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes four parameter names: relax_value,
     absolute_threshold, relative_threshold and overlap_mode. These names are
     case insensitive, and in each case except overlap_mode, the ParameterEntry
     must have type double. For overlap_mode, the ParameterEntry must have
     type Epetra_CombineMode.
  */
  int SetParameters(Teuchos::ParameterList& parameterlist);

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
      affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
      does not support transpose use, this method should return a value of -1.

      \param
       UseTranspose_in - (In) If true, multiply by the transpose of operator, otherwise just use operator.

      \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);};
  // @}

  // @{ Mathematical functions.
  // Applies the matrix to X, returns the result in Y.
  int Apply(const Epetra_MultiVector& X,
               Epetra_MultiVector& Y) const
  {
    return(Multiply(false,X,Y));
  }

  int Multiply(bool Trans, const Epetra_MultiVector& X,
               Epetra_MultiVector& Y) const{return A_->Multiply(Trans,X,Y);}

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! In this implementation, we use several existing attributes to determine how virtual
      method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(),
      the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
      if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
      be accounted for when doing a triangular solve.

    \param
           X - (In) A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
           Y - (Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Computes the estimated condition number and returns the value.
  double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                 const int MaxIters = 1550,
                 const double Tol = 1e-9,
                 Epetra_RowMatrix* Matrix_in = 0);

  //! Returns the computed estimated condition number, or -1.0 if not computed.
  double Condest() const
  {
    return(Condest_);
  }

  // @}
  // @{ Query methods

  //! Returns a character string describing the operator
  const char* Label() const {return(Label_);}

  //! Sets label for \c this object.
  int SetLabel(const char* Label_in)
  {
    strcpy(Label_,Label_in);
    return(0);
  }

  //! Returns 0.0 because this class cannot compute Inf-norm.
  double NormInf() const {return(0.0);};

  //! Returns false because this class cannot compute an Inf-norm.
  bool HasNormInf() const {return(false);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(A_->OperatorDomainMap());};

  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const{return(A_->OperatorRangeMap());};

  //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
  const Epetra_Comm & Comm() const{return(A_->Comm());};

  //! Returns a reference to the matrix to be preconditioned.
  const Epetra_RowMatrix& Matrix() const
  {
    return(*A_);
  }

  //! Prints on stream basic information about \c this object.
  virtual std::ostream& Print(std::ostream& os) const;

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const
  {
    return(NumInitialize_);
  }

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const
  {
    return(NumCompute_);
  }

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const
  {
    return(NumApplyInverse_);
  }

  //! Returns the time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the time spent in Compute().
  virtual double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const
  {
    return(0.0);
  }

  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

private:

  // @}
  // @{ Private methods

  //! Copy constructor (should never be used)
  Ifpack_IHSS(const Ifpack_IHSS& RHS) :
    Time_(RHS.Comm())
  {}

  //! operator= (should never be used)
  Ifpack_IHSS& operator=(const Ifpack_IHSS& /* RHS */)
  {
    return(*this);
  }

  //! Destroys all internal data
  void Destroy();

  //! Returns the result of a Ifpack_ILU forward/back solve on a Epetra_MultiVector X in Y.
  /*!
    \param In
    Trans -If true, solve transpose problem.
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y - (Out) A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;


#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  //! Returns the number of global matrix rows.
  int NumGlobalRows() const {return(A_->NumGlobalRows());};

  //! Returns the number of global matrix columns.
  int NumGlobalCols() const {return(A_->NumGlobalCols());};
#endif
  long long NumGlobalRows64() const {return(A_->NumGlobalRows64());};
  long long NumGlobalCols64() const {return(A_->NumGlobalCols64());};

  //! Returns the number of local matrix rows.
  int NumMyRows() const {return(A_->NumMyRows());};

  //! Returns the number of local matrix columns.
  int NumMyCols() const {return(A_->NumMyCols());};

  //! Power method for maximum eigenvalue estimation
  int PowerMethod(Epetra_Operator * Op,const int MaximumIterations,  double& lambda_max);

  //! Returns a reference to the matrix.
  /*  Epetra_RowMatrix& Matrix()
  {
    return(*A_);
    }*/

  // @}
  // @{ Internal data

  //! Pointer to the Epetra_RowMatrix to factorize
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A_;
  Teuchos::RefCountPtr<Epetra_RowMatrix> Aherm_;
  Teuchos::RefCountPtr<Epetra_RowMatrix> Askew_;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Pherm_;
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Pskew_;

  Teuchos::ParameterList List_;

  bool UseTranspose_;
  double Condest_;


  //! If \c true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Label of \c this object.
  char Label_[160];
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;

  //! Number of iterations in computing LambdaMax_
  int EigMaxIters_;
  //! Maximum-to-minimum eig ratio to use for Aherm.
  double EigRatio_;
  //! Contains an approximation to the largest eigenvalue of Aherm
  double LambdaMax_;
  //! Contains the shift parameter, alpha.
  double Alpha_;
  //! Number of IHSS iterations
  int NumSweeps_;



  //! Contains the number of successful call to ApplyInverse().
  mutable int NumApplyInverse_;
  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;
  //! Used for timing issues
  mutable Epetra_Time Time_;

};
#else

#ifdef _MSC_VER
#pragma message("No EpetraEXT")
#else
#warning "No EpetraEXT"
#endif

#endif
#endif /* IFPACK_IHSS_H */
