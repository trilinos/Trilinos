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

#ifndef TIFPACK_CHEBYSHEV_DECL_HPP
#define TIFPACK_CHEBYSHEV_DECL_HPP


#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_Condest.hpp"
#include "Tifpack_Parameters.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <Teuchos_TestForException.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <iostream>
#include <string>
#include <sstream>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Tifpack {

//! Tifpack::Chebyshev: class for preconditioning with Chebyshev polynomials

/*!
  The Tifpack::Chebyshev class enables the construction of preconditioners
  based on Chebyshev polynomials for a Tpetra::RowMatrix.
  Tifpack::Chebyshev is derived from the Tifpack::Preconditioner class, 
  which is itself derived from Tpetra::Operator.
  Therefore this object can be used as preconditioner everywhere an
  apply() method is required in the preconditioning step.

  The class is an adaptation of the routine ML_Cheby in Smoother/ml_smoother.hpp

<P> (04/04/06) Flops are not counted in the routine apply()

The list of parameters is
- EigRatio: "chebyshev: ratio eigenvalue"
this is the ratio to define the lower bound on the spectrum; lambda^* = LambdaMax / EigRatio;
a typical value used in ML is 30.0 (30.0 is the default value).
- LambdaMin: "chebyshev: min eigenvalue"
this is the smallest eigenvalue; this parameter is optional and is only
accessed to check whether the input matrix is equal to identity.
- LambdaMax: "chebyshev: max eigenvalue"
this is the largest eigenvalue of the matrix.
- PolyDegree:  "chebyshev: degree"
this is the polynomial degree.
- MinDiagonalValue:  "chebyshev: min diagonal value"
this defines the threshold for diagonal values under which they are not inverted
- ZeroStartingSolution:  "chebyshev: zero starting solution"
this flag allows to set a non-zero initial guess.

\author Ulrich Hetmaniuk. SNL 1416.

\date Last modified on 04-Apr-06.
*/
template<class MatrixType>
class Chebyshev : virtual public Tifpack::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> {

public:
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  // \name Constructors and Destructors
  //@{

  //! Chebyshev constructor with given Tpetra::RowMatrix input.
  explicit Chebyshev(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A);

  //! Chebyshev destructor.
  virtual ~Chebyshev();

  //@}

  //@{ \name Preconditioner computation methods

  //! Sets all the parameters for the preconditioner
  void setParameters(const Teuchos::ParameterList& params);

  //! Initialize
  void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(IsInitialized_);
  }

  //! Computes the preconditioner.
  void compute();

  //! If preconditioner is computed, this query returns true, otherwise it returns false.
  inline bool isComputed() const {
    return(IsComputed_);
  }

  //@}

  //! @name Methods implementing a Tpetra::Operator interface.
  //@{ 

  //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Tpetra::_MultiVector of dimension NumVectors containing result.
    */
  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                   Teuchos::ETransp mode = Teuchos::NO_TRANS,
                 Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
                 Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getRangeMap() const;

  bool hasTransposeApply() const;

  //! Applies the matrix to a Tpetra::MultiVector.
  /*! 
    \param 
    X - (In) A Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing the result.
    */
  void applyMat(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                      Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}

  //@{
  //! \name Mathematical functions.

  //! Computes the estimated condition number and returns the value.
  magnitudeType computeCondEst(CondestType CT = Cheap, 
                               LocalOrdinal MaxIters = 1550,
                               magnitudeType Tol = 1e-9,
                               const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix = Teuchos::null);

  //@}

  //@{ 
  //! \name Attribute accessor methods

  //! Returns the computed estimated condition number, or -1.0 if no computed.
  magnitudeType getCondEst() const;

  //! Returns the Tpetra::BlockMap object associated with the range of this matrix operator.
  const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! Returns a reference to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getMatrix() const;

  //! Returns the number of flops in the computation phase.
  double getComputeFlops() const;

  //! Returns the number of flops for the application of the preconditioner.
  double getApplyFlops() const;

  //! Returns the number of calls to initialize().
  int getNumInitialize() const;

  //! Returns the number of calls to compute().
  int getNumCompute() const;

  //! Returns the number of calls to apply().
  int getNumApply() const;

  //! Returns the time spent in initialize().
  double getInitializeTime() const;

  //! Returns the time spent in compute().
  double getComputeTime() const;

  //! Returns the time spent in apply().
  double getApplyTime() const;

  //@}

  //! @name Overridden from Teuchos::Describable 
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

  // @{ \name Utility methods

  //! Simple power method to compute lambda_max.
  static void PowerMethod(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator,
                         const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal,
                         const int MaximumIterations, 
                         Scalar& LambdaMax);

  //! Not currently implemented: Use CG to estimate lambda_min and lambda_max.
  static void CG(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator, 
                const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal, 
                const int MaximumIterations, 
                Scalar& lambda_min, Scalar& lambda_max);

  //@}

private:
  
  //! Copy constructor (should never be used)
  Chebyshev(const Chebyshev<MatrixType>& src);

  //! operator= (should never be used)
  Chebyshev<MatrixType>& operator=(const Chebyshev<MatrixType>& src);

  // @{ Internal data and parameters

  //! reference to the matrix to be preconditioned
  const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! Reference to the communicator object
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm_;
  //! Contains the inverse of diagonal elements of \c Matrix.
  mutable Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > InvDiagonal_;
  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! Contains the degree of Chebyshev polynomial.
  int PolyDegree_;
  //! The ratio such that [LambdaMax_ / EigRatio_, LambdaMax_], the interval of interest for the Chebyshev polynomial.
  Scalar EigRatio_;
  //! An approximation to the smallest eigenvalue.
  Scalar LambdaMin_;
  //! An approximation to the largest eigenvalue.
  Scalar LambdaMax_;
  //! The minimum value on the diagonal.
  Scalar MinDiagonalValue_;
  //! The estimated condition number
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  magnitudeType Condest_;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_;
  //! Contains the number of successful calls to initialize().
  int NumInitialize_;
  //! Contains the number of successful call to compute().
  int NumCompute_;
  //! Contains the number of successful call to apply().
  mutable int NumApply_;
  //! Contains the time for all successful calls to initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to apply().
  mutable double ApplyTime_;
  //! Contains the number of flops for compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for apply().
  mutable double ApplyFlops_;
  //! Number of local rows.
  size_t NumMyRows_;
  //! Number of global rows.
  global_size_t NumGlobalRows_;
  //! Number of global nonzeros.
  global_size_t NumGlobalNonzeros_;

  //@}

}; // class Chebyshev

}//namespace Tifpack

#endif // TIFPACK_CHEBYSHEV_DECL_HPP

