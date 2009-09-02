/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef TIFPACK_POINTRELAXATION_HPP
#define TIFPACK_POINTRELAXATION_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Preconditioner.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_Time.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Import.hpp"

#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos {
  class ParameterList;
}
class Tpetra_MultiVector;
class Tpetra_Vector;
class Tpetra_Map;
class Tpetra_Comm;
class Tpetra_CrsMatrix;

//! Tifpack_PointRelaxation: a class to define point relaxation preconditioners of for Tpetra_RowMatrix's.

/*! 
  The Tifpack_PointRelaxation class enables the construction of point
  relaxation
  preconditioners of an Tpetra_RowMatrix. Tifpack_PointRelaxation 
  is derived from 
  the Tifpack_Preconditioner class, which is itself derived from Tpetra_Operator.
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
where \f$A\f$ is a square, real matrix, and \f$x, b\f$ are two real
vectors. We begin with the decomposition
\f[
A = D - E - F
\f]
where \f$D\f$ is the diagonal of A, \f$-E\f$ is the strict lower part, and
\f$-F\f$ is the strict upper part. It is assumed that the diagonal entries
of \f$A\f$ are different from zero.

<P>Given an starting solution \f$x_0\f$, an iteration of the (damped) Jacobi
method can be written in matrix form as follows:
\f[
x_{k+1} = \omega D^{-1}(E + F) x_k + D_{-1}b,
\f]
for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter.

Using Tifpack_Jacobi, the user can apply the specified number of sweeps
(\f$k_{max}\f$), and the damping parameter. If only one sweep is used, then
the class simply applies the inverse of the diagonal of A to the input
vector.

<P>Given an starting solution \f$x_0\f$, an iteration of the (damped) GaussSeidel
method can be written in matrix form as follows:
\f[
(D - E) x_{k+1} = \omega F x_k + b,
\f]
for \f$k < k_{max}\f$, and \f$\omega \f$ a damping parameter. Equivalently,
the Gauss-Seidel preconditioner can be defined as
\f[
P_{GS}^{-1} = (D - E)^{-1}.
\f]
Clearly, the role of E and F can be interchanged. However,
Tifpack_GaussSeidel does not consider backward Gauss-Seidel methods.

<P>For a list of supported parameters, please refer to page \ref ifp_params.

<P>The complete list of supported parameters is reported in page \ref ifp_params. For a presentation of basic relaxation schemes, please refer to page
\ref Tifpack_PointRelaxation.

\author Michael Heroux, SNL 9214.

\date Last modified on 22-Jan-05.
  
*/
class Tifpack_PointRelaxation : public Tifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Tifpack_PointRelaxation constructor with given Tpetra_RowMatrix.
  /*! Creates an instance of Tifpack_PointRelaxation class.
   *
   * \param
   * Matrix - (In) Pointer to matrix to precondition.
   */
  Tifpack_PointRelaxation(const Tpetra_RowMatrix* Matrix);

  //! Destructor.
  virtual ~Tifpack_PointRelaxation() {}

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

  //! Applies the matrix to an Tpetra_MultiVector.
  /*! 
    \param 
    X - (In) A Tpetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y - (Out) A Tpetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
  virtual inline int Apply(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const;

  //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param
    X - (In) A Tpetra_MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Tpetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning This routine is NOT AztecOO complaint.
    */
  virtual int ApplyInverse(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix (not implemented)
  virtual double NormInf() const
  {
    return(-1.0);
  }
  //@}

  //@{ \name Atribute access functions

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

  //! Returns a pointer to the Tpetra_Comm communicator associated with this operator.
  virtual const Tpetra_Comm & Comm() const;

  //! Returns the Tpetra_Map object associated with the domain of this operator.
  virtual const Tpetra_Map & OperatorDomainMap() const;

  //! Returns the Tpetra_Map object associated with the range of this operator.
  virtual const Tpetra_Map & OperatorRangeMap() const;

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

  virtual const Tpetra_RowMatrix& Matrix() const 
  {
    return(*Matrix_);
  }

  //! Computes the condition number estimates and returns the value.
  virtual double Condest(const Tifpack_CondestType CT = Tifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
			 Tpetra_RowMatrix* Matrix = 0);

  //! Returns the condition number estimate, or -1.0 if not computed.
  virtual double Condest() const
  {
    return(Condest_);
  }

  //! Sets all the parameters for the preconditioner
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Prints object to an output stream
  virtual ostream& Print(ostream & os) const;

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
  virtual int ApplyInverseJacobi(const Tpetra_MultiVector& X, 
                                 Tpetra_MultiVector& Y) const;

  //! Applies the Gauss-Seidel preconditioner to X, returns the result in Y.
  virtual int ApplyInverseGS(const Tpetra_MultiVector& X, 
                              Tpetra_MultiVector& Y) const;

  virtual int ApplyInverseGS_RowMatrix(const Tpetra_MultiVector& X, 
                                        Tpetra_MultiVector& Y) const;

  virtual int ApplyInverseGS_CrsMatrix(const Tpetra_CrsMatrix* A,
                                        const Tpetra_MultiVector& X, 
                                        Tpetra_MultiVector& Y) const;

  virtual int ApplyInverseGS_FastCrsMatrix(const Tpetra_CrsMatrix* A,
                                            const Tpetra_MultiVector& X, 
                                            Tpetra_MultiVector& Y) const;

  //! Applies the symmetric Gauss-Seidel preconditioner to X, returns the result in Y.
  virtual int ApplyInverseSGS(const Tpetra_MultiVector& X, 
                              Tpetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS_RowMatrix(const Tpetra_MultiVector& X, 
                                        Tpetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS_CrsMatrix(const Tpetra_CrsMatrix* A,
                                        const Tpetra_MultiVector& X, 
                                        Tpetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS_FastCrsMatrix(const Tpetra_CrsMatrix* A,
                                            const Tpetra_MultiVector& X, 
                                            Tpetra_MultiVector& Y) const;
  //@}

private:
  
  //! Sets the label.
  virtual void SetLabel();

  //! Copy constructor (PRIVATE, should not be used)
  Tifpack_PointRelaxation(const Tifpack_PointRelaxation& rhs)
  {}
  
  //! operator = (PRIVATE, should not be used)
  Tifpack_PointRelaxation& operator=(const Tifpack_PointRelaxation& rhs)
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
  //! If true, Compute() also computes the condition number estimate.
  bool ComputeCondest_;
  //! Contains the label of this object.
  string Label_;
  int PrecType_;
  double MinDiagonalValue_;
  // @}

  // @{ Other data
  //! Number of local rows.
  int NumMyRows_;
  //! Number of local nonzeros.
  int NumMyNonzeros_;
  //! Number of global rows.
  int NumGlobalRows_;
  //! Number of global nonzeros.
  int NumGlobalNonzeros_;
  //! Pointers to the matrix to be preconditioned.
  Teuchos::RefCountPtr<const Tpetra_RowMatrix> Matrix_;
  //! Importer for parallel GS and SGS
  Teuchos::RefCountPtr<Tpetra_Import> Importer_;
  //! Contains the diagonal elements of \c Matrix.
  mutable Teuchos::RefCountPtr<Tpetra_Vector> Diagonal_;
  //! Time object to track timing.
  Teuchos::RefCountPtr<Tpetra_Time> Time_;
  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  //! Backward-Mode Gauss Seidel 
  bool DoBackwardGS_;
  // @}

  

};

#endif // TIFPACK_POINTRELAXATION_HPP
