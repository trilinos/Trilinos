/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef IFPACK_POLYNOMIAL_H
#define IFPACK_POLYNOMIAL_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Teuchos {
  class ParameterList;
}

class Epetra_MultiVector;
class Epetra_Vector;
class Epetra_Map;
class Epetra_Comm;
class Epetra_Time;
class Epetra_Vector;
class Epetra_Operator;
class Epetra_RowMatrix;

#ifdef HAVE_IFPACK_EPETRAEXT
class EpetraExt_PointToBlockDiagPermute;
#endif

//! Ifpack_Polynomial: class for preconditioning with least squares polynomials in Ifpack

/*!
  The Ifpack_Polynomial class enables the construction of preconditioners
  based on least squares polynomials for an Epetra_RowMatrix. Similar to
  Ifpack_Chebyshev, Ifpack_Polynomial is designed for indefinite linear systems.

The list of parameters is
- RealEigRatio_             = List.get("polynomial: real eigenvalue ratio", RealEigRatio_);
- ImagEigRatio_             = List.get("polynomial: imag eigenvalue ratio", ImagEigRatio_);
these are the ratios to estimate the bounds on the spectrum (if the eigenvalue bounds aren't already given)
- LambdaRealMin_            = List.get("polynomial: min real part", LambdaRealMin_);
- LambdaRealMax_            = List.get("polynomial: max real part", LambdaRealMax_);
- LambdaImagMin_            = List.get("polynomial: min imag part", LambdaImagMin_);
- LambdaImagMax_            = List.get("polynomial: max imag part", LambdaImagMax_);
these values define the rectangular region on which the polynomial will be minimized.
- PolyDegree_           = List.get("polynomial: degree",PolyDegree_);
this is the polynomial degree.
- LSPointsReal_             = List.get("polynomial: real interp points",LSPointsReal_);
- LSPointsImag_             = List.get("polynomial: imag interp points",LSPointsImag_);
these are the number of points used in real and imaginary directions which will be used in the least squares problem.
- MinDiagonalValue_     = List.get("polynomial: min diagonal value", MinDiagonalValue_);
this defines the threshold for diagonal values under which they are not inverted
- ZeroStartingSolution_ = List.get("polynomial: zero starting solution", ZeroStartingSolution_);
this flag allows to set a non-zero initial guess.

\author Paul Tsuji, May 2013.

*/

class Ifpack_Polynomial : public Ifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack_Polynomial constructor with given Epetra_Operator/Epetra_RowMatrix.
  /*! Creates an instance of Ifpack_Polynomial class.
   *
   * \param
   * Matrix - (In) Pointer to the operator to precondition.
   */
  Ifpack_Polynomial(const Epetra_Operator* Matrix);

  //! Ifpack_Polynomial constructor with given Epetra_Operator/Epetra_RowMatrix.
  /*! Creates an instance of Ifpack_Polynomial class.
   *
   * \param
   * Matrix - (In) Pointer to the matrix to precondition.
   */
  Ifpack_Polynomial(const Epetra_RowMatrix* Matrix);

  //! Destructor.
  virtual ~Ifpack_Polynomial() {};

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
  virtual inline int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

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
 
  //@{ \name Miscellaneous

  virtual const Epetra_RowMatrix& Matrix() const 
  {
    return(*Matrix_);
  }

  //! Computes the condition number estimates and returns the value.
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
                         Epetra_RowMatrix* Matrix_in = 0);

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
  // @{ \name Utility methods

  //! Simple power method to compute lambda_max.
  static int PowerMethod(const Epetra_Operator& Operator,
                         const Epetra_Vector& InvPointDiagonal,
                         const int MaximumIterations, 
                         double& LambdaMax);

  //! Uses AztecOO's CG to estimate lambda_min and lambda_max.
  static int CG(const Epetra_Operator& Operator, 
                const Epetra_Vector& InvPointDiagonal, 
                const int MaximumIterations, 
                double& lambda_min, double& lambda_max);
  
#ifdef HAVE_IFPACK_EPETRAEXT
  //! Uses AztecOO's CG to estimate lambda_min and lambda_max.
  // WARNING: This only works in Block Mode.
  int CG(const int MaximumIterations, 
         double& lambda_min, double& lambda_max);
  //! Simple power method to compute lambda_max.
  // WARNING: This only works in Block Mode.
  int PowerMethod(const int MaximumIterations,double& lambda_max);
#endif
private:
  
  // @}
  // @{ \name Private methods
  
  //! Sets the label.
  virtual void SetLabel();

  //! Copy constructor (PRIVATE, should not be used)
  Ifpack_Polynomial(const Ifpack_Polynomial& rhs)
  {}
  
  //! operator = (PRIVATE, should not be used)
  Ifpack_Polynomial& operator=(const Ifpack_Polynomial& rhs)
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
  //! Contains the degree of the least squares polynomial.
  int PolyDegree_;
  //! Contains the number of discretization points of the least squares problem.
  int LSPointsReal_, LSPointsImag_;
  //! If true, use the tranpose of \c Matrix_.
  bool UseTranspose_;
  //! Contains the estimated condition number
  double Condest_;
  //! If true, Compute() also computes the condition number estimate.
  bool ComputeCondest_;
  //! Contains the ratio such that the rectangular domain in the complex plane
  //! is [-LambdaRealMax_ / EigRatio_, LambdaRealMax_] x [-LambdaRealMax_ / ImagEigRatio_, LambdaRealMax_ / ImagEigRatio_]
  double RealEigRatio_, ImagEigRatio_;
  //! Max number of iterations to use in eigenvalue estimation (if automatic).
  int EigMaxIters_;  
  //! Contains the label of this object.
  string Label_;
  //! Bounds on the spectrum
  double LambdaRealMin_, LambdaRealMax_, LambdaImagMin_, LambdaImagMax_;
  //! Contains the minimum value on the diagonal.
  double MinDiagonalValue_;
  //! coefficients of the polynomial
  std::vector<double> coeff_;

  // @{ Other data
  //! Number of local rows.
  int NumMyRows_;
  //! Number of local nonzeros.
  int NumMyNonzeros_;
  //! Number of global rows.
  int NumGlobalRows_;
  //! Number of global nonzeros.
  int NumGlobalNonzeros_;
  //! Pointers to the matrix to be preconditioned as an Epetra_Operator.
  Teuchos::RefCountPtr<const Epetra_Operator> Operator_;
  //! Pointers to the matrix to be preconditioned as an Epetra_RowMatrix.
  Teuchos::RefCountPtr<const Epetra_RowMatrix> Matrix_;
  //! Contains the inverse of diagonal elements of \c Matrix.
  mutable Teuchos::RefCountPtr<Epetra_Vector> InvDiagonal_;
  //! Use Block Preconditioning
  bool UseBlockMode_;
#ifdef HAVE_IFPACK_EPETRAEXT
  //! Max/Min Ration for autocomputing eigenvalues
  Teuchos::ParameterList BlockList_;
  Teuchos::RefCountPtr<EpetraExt_PointToBlockDiagPermute> InvBlockDiagonal_;
#endif


  //! Run on the normal equations
  bool SolveNormalEquations_;

  //! If \c true, the Operator_ is an Epetra_RowMatrix.
  bool IsRowMatrix_;
  //! Time object to track timing.
  Teuchos::RefCountPtr<Epetra_Time> Time_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;

  // @}

};


#endif // IFPACK_POLYNOMIAL_H
