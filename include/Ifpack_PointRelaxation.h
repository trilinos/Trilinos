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

#ifndef IFPACK_POINTRELAXATION_H
#define IFPACK_POINTRELAXATION_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"

#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Import.h"

#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos {
  class ParameterList;
}
class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_CrsMatrix;

//! Ifpack_PointRelaxation: a class to define point relaxation preconditioners of for Epetra_RowMatrix's.

/*!
  The Ifpack_PointRelaxation class enables the construction of point
  relaxation
  preconditioners of an Epetra_RowMatrix. Ifpack_PointRelaxation
  is derived from
  the Ifpack_Preconditioner class, which is itself derived from Epetra_Operator.
  Therefore this object can be used as preconditioner everywhere an
  ApplyInverse() method is required in the preconditioning step.

This class enables the construction of the following simple preconditioners:
- Jacobi;
- Gauss-Seidel;
- symmetric Gauss-Seidel.

<P>We now briefly describe the main features of the above preconditioners.
Consider a linear system of type
\f[
A x = b,
\f]
where \f$ A\f$ is a square, real matrix, and \f$ x, b\f$ are two real
vectors. We begin with the decomposition
\f[
A = D - E - F
\f]
where \f$ D\f$ is the diagonal of A, \f$ -E\f$ is the strict lower part, and
\f$ -F\f$ is the strict upper part. It is assumed that the diagonal entries
of \f$ A\f$ are different from zero.

<P>Given an starting solution \f$ x_0\f$, an iteration of the (damped) Jacobi
method can be written in matrix form as follows:
\f[
x_{k+1} = \omega D^{-1}(E + F) x_k + D_{-1}b,
\f]
for \f$ k < k_{max}\f$, and \f$\omega \f$ a damping parameter.

Using Ifpack_Jacobi, the user can apply the specified number of sweeps
(\f$ k_{max}\f$), and the damping parameter. If only one sweep is used, then
the class simply applies the inverse of the diagonal of A to the input
vector.

<P>Given an starting solution \f$ x_0\f$, an iteration of the (damped) GaussSeidel
method can be written in matrix form as follows:
\f[
(D - E) x_{k+1} = \omega F x_k + b,
\f]
for \f$ k < k_{max}\f$, and \f$\omega \f$ a damping parameter. Equivalently,
the Gauss-Seidel preconditioner can be defined as
\f[
P_{GS}^{-1} = (D - E)^{-1}.
\f]
Clearly, the role of E and F can be interchanged. However,
Ifpack_GaussSeidel does not consider backward Gauss-Seidel methods.

<P>For a list of supported parameters, please refer to page \ref ifp_params.

<P>The complete list of supported parameters is reported in page \ref ifp_params. For a presentation of basic relaxation schemes, please refer to page
\ref Ifpack_PointRelaxation.

\author Marzio Sala, SNL 9214.

\date Last modified on 22-Jan-05.

*/
class Ifpack_PointRelaxation : public Ifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack_PointRelaxation constructor with given Epetra_RowMatrix.
  /*! Creates an instance of Ifpack_PointRelaxation class.
   *
   * \param Matrix - (In) Pointer to matrix to precondition.
   */
  Ifpack_PointRelaxation(const Epetra_RowMatrix* Matrix);

  //! Destructor.
  virtual ~Ifpack_PointRelaxation() {}

  //@}

  /*! This flag can be used to apply the preconditioner to the transpose of
   * the input operator.
   *
   * \return Integer error code, set to 0 if successful.
   * Set to -1 if this implementation does not support transpose.
    */
  virtual inline int SetUseTranspose(bool UseTranspose_in)
  {
    UseTranspose_ = UseTranspose_in;
    return(0);
  }

  //@}

  //@{ \name Mathematical functions.

  //! Applies the matrix to an Epetra_MultiVector.
  /*!
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param
    Y - (Out) A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
  virtual inline int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  {
    if (IsComputed() == false)
      IFPACK_CHK_ERR(-3);

    if (X.NumVectors() != Y.NumVectors())
      IFPACK_CHK_ERR(-2);

    IFPACK_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
    return(0);
  }

  //! Applies the preconditioner to X, returns the result in Y.
  /*!
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning This routine is NOT AztecOO complaint.
    */
  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix (not implemented)
  virtual double NormInf() const
  {
    return(-1.0);
  }
  //@}

  //@{ \name Attribute access functions

  virtual const char * Label() const
  {
    return(Label_.c_str());
  }

  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const
  {
    return(UseTranspose_);
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const;

  virtual int Initialize();

  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual inline bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Computes the preconditioners.
  virtual int Compute();

  //@}

  //@{ \name Miscellaneous

  virtual const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  //! Computes the condition number estimates and returns the value.
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
                         Epetra_RowMatrix* Matrix = 0);

  //! Returns the condition number estimate, or -1.0 if not computed.
  virtual double Condest() const
  {
    return(Condest_);
  }

  //! Sets all the parameters for the preconditioner
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Prints object to an output stream
  virtual std::ostream& Print(std::ostream & os) const;

  //@}

  //@{ \name Timing and flop count

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

  //! Returns the number of flops in the computation phase.
  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  //! Returns the number of flops for the application of the preconditioner.
  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  // @}

private:

  // @{ Application of the preconditioner

  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  virtual int ApplyInverseJacobi(const Epetra_MultiVector& X,
                                 Epetra_MultiVector& Y) const;

  //! Applies the Gauss-Seidel preconditioner to X, returns the result in Y.
  virtual int ApplyInverseGS(const Epetra_MultiVector& X,
                              Epetra_MultiVector& Y) const;

  virtual int ApplyInverseGS_RowMatrix(const Epetra_MultiVector& X,
                                        Epetra_MultiVector& Y) const;

  virtual int ApplyInverseGS_CrsMatrix(const Epetra_CrsMatrix* A,
                                        const Epetra_MultiVector& X,
                                        Epetra_MultiVector& Y) const;

  virtual int ApplyInverseGS_FastCrsMatrix(const Epetra_CrsMatrix* A,
                                            const Epetra_MultiVector& X,
                                            Epetra_MultiVector& Y) const;
  virtual int ApplyInverseGS_LocalFastCrsMatrix(const Epetra_CrsMatrix* A,
                                            const Epetra_MultiVector& X,
                                            Epetra_MultiVector& Y) const;

  //! Applies the symmetric Gauss-Seidel preconditioner to X, returns the result in Y.
  virtual int ApplyInverseSGS(const Epetra_MultiVector& X,
                              Epetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS_RowMatrix(const Epetra_MultiVector& X,
                                        Epetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS_CrsMatrix(const Epetra_CrsMatrix* A,
                                        const Epetra_MultiVector& X,
                                        Epetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS_FastCrsMatrix(const Epetra_CrsMatrix* A,
                                            const Epetra_MultiVector& X,
                                            Epetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS_LocalFastCrsMatrix(const Epetra_CrsMatrix* A,
                                            const Epetra_MultiVector& X,
                                            Epetra_MultiVector& Y) const;


  //@}

private:

  //! Sets the label.
  virtual void SetLabel();

  //! Copy constructor (PRIVATE, should not be used)
  Ifpack_PointRelaxation(const Ifpack_PointRelaxation& /* rhs */)
  {}

  //! operator = (PRIVATE, should not be used)
  Ifpack_PointRelaxation& operator=(const Ifpack_PointRelaxation& /* rhs */)
  {
    return(*this);
  }

  // @{ Initializations, timing and flops
  //! If \c true, the preconditioner has been computed successfully.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_;
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
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
  // @}

  // @{ Settings
  //! Number of application of the preconditioner (should be greater than 0).
  int NumSweeps_;
  //! Damping factor.
  double DampingFactor_;
  //! If true, use the tranpose of \c Matrix_.
  bool UseTranspose_;
  //! Contains the estimated condition number
  double Condest_;
#if 0
  // Unused; commented out to avoid build warnings

  //! If true, Compute() also computes the condition number estimate.
  bool ComputeCondest_;
#endif // 0
  //! Contains the label of this object.
  std::string Label_;
  int PrecType_;
  double MinDiagonalValue_;
  // @}

  // @{ Other data
  //! Number of local rows.
  int NumMyRows_;
  //! Number of local nonzeros.
  int NumMyNonzeros_;
  //! Number of global rows.
  long long NumGlobalRows_;
  //! Number of global nonzeros.
  long long NumGlobalNonzeros_;
  //! Pointers to the matrix to be preconditioned.
  Teuchos::RefCountPtr<const Epetra_RowMatrix> Matrix_;
  //! Importer for parallel GS and SGS
  Teuchos::RefCountPtr<Epetra_Import> Importer_;
  //! Contains the diagonal elements of \c Matrix.
  mutable Teuchos::RefCountPtr<Epetra_Vector> Diagonal_;
  //! Time object to track timing.
  Teuchos::RefCountPtr<Epetra_Time> Time_;
  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  //! Backward-Mode Gauss Seidel
  bool DoBackwardGS_;
  //! Do L1 Jacobi/GS/SGS
  bool DoL1Method_;
  //! Eta parameter for modified L1 method
  double L1Eta_;

  //! Number of (local) unknowns for local smoothing
  int NumLocalSmoothingIndices_;
  //! List of (local) unknowns for local smoothing (if any)
  int * LocalSmoothingIndices_;

  // @}



};

#endif // IFPACK_POINTRELAXATION_H
