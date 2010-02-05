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

#ifndef TIFPACK_CHEBYSHEV_HPP
#define TIFPACK_CHEBYSHEV_HPP


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
  based on Chebyshev polynomials for a Tpetra::CrsMatrix.
  Tifpack::Chebyshev is derived from the Tifpack::Preconditioner class, 
  which is itself derived from Tpetra::Operator.
  Therefore this object can be used as preconditioner everywhere an
  apply() method is required in the preconditioning step.

  The class is an adaptation of the routine ML_Cheby in Smoother/ml_smoother.hpp

<P> (04/04/06) Flops are not counted in the routine apply()

The list of parameters is
- EigRatio_             = List.get("chebyshev: ratio eigenvalue", EigRatio_);
this is the ratio to define the lower bound on the spectrum; lambda^* = LambdaMax_ / EigRatio_;
a typical value used in ML is 30.0 (30.0 is the default value).
- LambdaMin_            = List.get("chebyshev: min eigenvalue", LambdaMin_);
this is the smallest eigenvalue; this parameter is optional and is only
accessed to check whether the input matrix is equal to identity.
- LambdaMax_            = List.get("chebyshev: max eigenvalue", LambdaMax_);
this is the largest eigenvalue of the matrix.
- PolyDegree_           = List.get("chebyshev: degree",PolyDegree_);
this is the polynomial degree.
- MinDiagonalValue_     = List.get("chebyshev: min diagonal value", MinDiagonalValue_);
this defines the threshold for diagonal values under which they are not inverted
- ZeroStartingSolution_ = List.get("chebyshev: zero starting solution", ZeroStartingSolution_);
this flag allows to set a non-zero initial guess.

\author Ulrich Hetmaniuk. SNL 1416.

\date Last modified on 04-Apr-06.
*/
template<class Scalar,class LocalOrdinal = int,class GlobalOrdinal = LocalOrdinal,class Node = Kokkos::DefaultNode::DefaultNodeType,class LocalMatVec = Kokkos::DefaultSparseMultiply<Scalar,LocalOrdinal,Node>,class LocalMatSolve = Kokkos::DefaultSparseSolve<Scalar,LocalOrdinal,Node> >
class Chebyshev : virtual public Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  // \name Constructors and Destructors
  //@{

  //! Chebyshev constructor with given Tpetra::RowMatrix input.
  explicit Chebyshev(const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> >& Matrix);

  //! Chebyshev destructor.
  virtual ~Chebyshev();

  //@}

  //@{ \name Preconditioner computation methods

  //! Sets all the parameters for the preconditioner
  void setParameters(Teuchos::ParameterList& params);

  //! Initialize
  void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(IsInitialized_);
  }

  //! Computes the preconditioner.
  void compute();

  //! If preconditioner is completed, this query returns true, otherwise it returns false.
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

  //! Uses AztecOO's CG to estimate lambda_min and lambda_max.
  static void CG(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator, 
                const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal, 
                const int MaximumIterations, 
                Scalar& lambda_min, Scalar& lambda_max);

  //@}

private:
  
  //! Copy constructor (should never be used)
  Chebyshev(const Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>& rhs);

  //! operator= (should never be used)
  Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>& operator=(const Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>& rhs);

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


//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::Chebyshev(const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> >& A)
: A_(A),
  Comm_(A->getRowMap()->getComm()),
  Time_( Teuchos::rcp( new Teuchos::Time("Tifpack::Chebyshev") ) ),
  PolyDegree_(1),
  EigRatio_(30.0),
  LambdaMin_(0.0),
  LambdaMax_(100.0),
  MinDiagonalValue_(0.0),
  ZeroStartingSolution_(true),
  Condest_(-1.0),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0),
  ComputeFlops_(0.0),
  ApplyFlops_(0.0),
  NumMyRows_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0)
{
  TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, 
      Teuchos::typeName(*this) << "::Chebyshev(): input matrix reference was null.");
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~Chebyshev() {
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::setParameters(Teuchos::ParameterList& List) {

  Tifpack::GetParameter(List, "chebyshev: ratio eigenvalue", EigRatio_);
  Tifpack::GetParameter(List, "chebyshev: min eigenvalue", LambdaMin_);
  Tifpack::GetParameter(List, "chebyshev: max eigenvalue", LambdaMax_);
  Tifpack::GetParameter(List, "chebyshev: degree",PolyDegree_);
  Tifpack::GetParameter(List, "chebyshev: min diagonal value", MinDiagonalValue_);
  Tifpack::GetParameter(List, "chebyshev: zero starting solution", ZeroStartingSolution_);

  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>* ID = 0;
  Tifpack::GetParameter(List, "chebyshev: operator inv diagonal", ID);

  if (ID != 0) {
    InvDiagonal_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(*ID) );
  }
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
const Teuchos::RCP<const Teuchos::Comm<int> > & 
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComm() const{
  return(Comm_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getMatrix() const {
  return(A_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getDomainMap() const {
  return A_->getDomainMap();
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRangeMap() const {
  return A_->getRangeMap();
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
bool Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApply() const {
  return true;
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
int Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
int Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
int Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getApplyTime() const {
  return(ApplyTime_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComputeFlops() const {
  return(ComputeFlops_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getApplyFlops() const {
  return(ApplyFlops_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getCondEst() const {
  return(Condest_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::computeCondEst(
                     CondestType CT,
                     LocalOrdinal MaxIters, 
                     magnitudeType Tol,
                     const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix) {
  if (!isComputed()) // cannot compute right now
    return(-1.0);

  // always compute it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Tifpack::Condest(*this, CT, MaxIters, Tol, matrix);

  return(Condest_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(
          const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                Teuchos::ETransp mode,
                 Scalar alpha,
                 Scalar beta) const {
  TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error, 
      "Tifpack::Chebyshev::apply() ERROR, not yet computed.");

  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Tifpack::Chebyshev::apply() ERROR: X.getNumVectors() != Y.getNumVectors().");

  Time_->start();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  Teuchos::RCP< const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xcopy;
  if (&(X.get2dView()[0][0]) == &(Y.get2dView()[0][0]))
    Xcopy = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > xView = Xcopy->get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > yView = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<const Scalar> invdiag = InvDiagonal_->get1dView();

  size_t nVecs = Y.getNumVectors();

  //--- Do a quick solve when the matrix is identity
  if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_)) {
    if (nVecs == 1) {
      Teuchos::ArrayRCP<Scalar> y = yView[0];
      Teuchos::ArrayRCP<const Scalar> x = xView[0];
      for (size_t i = 0; i < NumMyRows_; ++i)
        y[i] = x[i]*invdiag[i];
    }
    else {
      for (size_t i = 0; i < NumMyRows_; ++i) {
        const Scalar& coeff = invdiag[i];
        for (size_t k = 0; k < nVecs; ++k)
          yView[k][i] = xView[k][i] * coeff;
      }
    } // if (nVec == 1)
    return;
  } // if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_))

  //--- initialize coefficients
  // Note that delta stores the inverse of ML_Cheby::delta
  Scalar alpha_cheby = LambdaMax_ / EigRatio_;
  Scalar beta_cheby = 1.1 * LambdaMax_;
  Scalar delta = 2.0 / (beta_cheby - alpha_cheby);
  Scalar theta = 0.5 * (beta_cheby + alpha_cheby);
  Scalar s1 = theta * delta;

  //--- Define vectors
  // In ML_Cheby, V corresponds to pAux and W to dk
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> V(X);
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> W(X);

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > vView = V.get2dView();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > wView = W.get2dViewNonConst();

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();

  Scalar oneOverTheta = one/theta;

  // Do the smoothing when block scaling is turned OFF
  // --- Treat the initial guess
  if (ZeroStartingSolution_ == false) {
    A_->apply(Y, V);
    // compute W = invDiag * ( X - V )/ Theta
    if (nVecs == 1) {
      Teuchos::ArrayRCP<const Scalar> x = xView[0];
      Teuchos::ArrayRCP<Scalar> w = wView[0];
      Teuchos::ArrayRCP<const Scalar> v = vView[0];
      for (size_t i = 0; i < NumMyRows_; ++i)
        w[i] = invdiag[i] * (x[i] - v[i]) * oneOverTheta;
    }
    else {
      for (size_t k = 0; k < nVecs; ++k) {
        Teuchos::ArrayRCP<Scalar> wk = wView[k];
        Teuchos::ArrayRCP<const Scalar> vk = vView[k];
        for (size_t i = 0; i < NumMyRows_; ++i) {
          Scalar coeff = invdiag[i]*oneOverTheta;
          wk[i] = (xView[k][i] - (vk[i])) * coeff;
        }
      }
    } // if (nVec == 1)
    // Update the vector Y
    Y.update(one, W, one);
  }
  else {
    // compute W = invDiag * X / Theta
    if (nVecs == 1) {
      Teuchos::ArrayRCP<const Scalar> x= xView[0];
      Teuchos::ArrayRCP<Scalar> w = wView[0];
      Teuchos::ArrayRCP<Scalar> y = yView[0];
      for (size_t i = 0; i < NumMyRows_; ++i) {
        w[i] = invdiag[i] * x[i] * oneOverTheta;
        y[i] = w[i];
      }
    }
    else {
      for (size_t k = 0; k < nVecs; ++k) {
        for (size_t i = 0; i < NumMyRows_; ++i) {
          Scalar coeff = invdiag[i]*oneOverTheta;
          wView[k][i] = xView[k][i] * coeff;
          yView[k][i] = wView[k][i];
        }
      }
    } // if (nVec == 1)
  } // if (ZeroStartingSolution_ == false)

  //--- Apply the polynomial
  Scalar rhok = 1.0/s1, rhokp1;
  Scalar dtemp1, dtemp2;
  int degreeMinusOne = PolyDegree_ - 1;
  for (int deg = 0; deg < degreeMinusOne; ++deg) {
    A_->apply(Y, V);
    rhokp1 = 1.0 / (2.0*s1 - rhok);
    dtemp1 = rhokp1 * rhok;
    dtemp2 = 2.0 * rhokp1 * delta;
    rhok = rhokp1;
    // compute W = dtemp1 * W
    W.scale(dtemp1);
    // compute W = W + dtemp2 * invDiag * ( X - V )
    for (size_t k = 0; k < nVecs; ++k) {
      for (size_t i = 0; i < NumMyRows_; ++i) {
        Scalar coeff = invdiag[i]*dtemp2;
        wView[k][i] += (xView[k][i] - (vView[k][i])) * coeff;
      }
    }
    // Update the vector Y
    Y.update(one, W, one);
  } // for (deg = 0; deg < degreeMinusOne; ++deg)

  // Flops are updated in each of the following. 

  ++NumApply_;
  Time_->stop();
  ApplyTime_ += Time_->totalElapsedTime();
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::applyMat(
       const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode) const
{
  TEST_FOR_EXCEPTION(isComputed() == false, std::runtime_error,
     "Tifpack::Chebyshev::applyMat() ERROR: isComputed() must be true prior to calling applyMat().");
  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Tifpack::Chebyshev::applyMat() ERROR: X.getNumVectors() != Y.getNumVectors().");
  A_->apply(X, Y, mode);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::initialize() {
  IsInitialized_ = false;

  TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error, 
      "Tifpack::Chebyshev::initialize ERROR: A_ == Teuchos::null");

  TEST_FOR_EXCEPTION(A_->getGlobalNumRows() != A_->getGlobalNumCols(), std::runtime_error,
     "Tifpack::Chebyshev::initialize ERROR: only square matrices are supported");

  NumMyRows_ = A_->getNodeNumRows();
  NumGlobalRows_ = A_->getGlobalNumRows();
  NumGlobalNonzeros_ = A_->getGlobalNumEntries();

  ++NumInitialize_;
  Time_->stop();
  InitializeTime_ += Time_->totalElapsedTime();
  IsInitialized_ = true;
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::compute()
{
  if (!isInitialized()) {
    initialize();
  }

  Time_->start(true);

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  TEST_FOR_EXCEPTION(PolyDegree_ <= 0, std::runtime_error,
    "Tifpack::Chebyshev::compute() ERROR: PolyDegree_ must be at least one");
  
  if (InvDiagonal_ == Teuchos::null)
  {
    InvDiagonal_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A_->getRowMap()) );

    A_->getLocalDiagCopy(*InvDiagonal_);

    // Inverse diagonal elements
    // Replace zeros with 1.0
    Teuchos::ArrayRCP<Scalar> diagvals = InvDiagonal_->get1dViewNonConst();
    for (size_t i = 0 ; i < NumMyRows_ ; ++i) {
      Scalar& diag = diagvals[i];
      if (Teuchos::ScalarTraits<Scalar>::magnitude(diag) < MinDiagonalValue_)
        diag = MinDiagonalValue_;
      else
        diag = 1.0 / diag;
    }
  }
  // otherwise the inverse of the diagonal has been given by the user

  ComputeFlops_ += NumMyRows_;

  ++NumCompute_;
  Time_->stop();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::
PowerMethod(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator, 
            const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal, 
            const int MaximumIterations, 
            Scalar& lambda_max)
{
  // this is a simple power method
  lambda_max = 0.0;
  Teuchos::Array<Scalar> RQ_top(1), RQ_bottom(1);
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(Operator.getDomainMap());
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> y(Operator.getRangeMap());
  x.Random();
  Teuchos::Array<Scalar> norms(x.getNumVectors());
  x.norm2(norms());

  x.scale(1.0 / norms[0]);

  Scalar one = Teuchos::ScalarTraits<Scalar>::one();
  Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  for (int iter = 0; iter < MaximumIterations; ++iter)
  {
    Operator.apply(x, y);
    y.elementWiseMultiply(1.0, InvPointDiagonal, y, 0.0);
    y.dot(x, RQ_top());
    x.dot(x, RQ_bottom());
    lambda_max = RQ_top[0] / RQ_bottom[0];
    y.norm2(norms());
    TEST_FOR_EXCEPTION(norms[0] == zero, std::runtime_error, "Tifpack::Chebyshev::PowerMethod ERROR, norm == 0");
    x.update( one / norms[0], y, zero);
  }
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::
CG(const Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Operator, 
            const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& InvPointDiagonal, 
   const int MaximumIterations, 
   Scalar& lambda_min, Scalar& lambda_max)
{
#ifdef HAVE_TIFPACK_AZTECOO
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> x(Operator.getDomainMap());
  Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> y(Operator.getRangeMap());
  x.Random();
  y.putScalar(0.0);

  Tpetra::LinearProblem LP(const_cast<Tpetra::Operator*>(&Operator), &x, &y);
  AztecOO solver(LP);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, AZ_none);

  Tifpack_DiagPreconditioner diag(Operator.OperatorDomainMap(),
                                 Operator.OperatorRangeMap(),
                                 InvPointDiagonal);
  solver.SetPrecOperator(&diag);
  solver.Iterate(MaximumIterations, 1e-10);

  const double* status = solver.GetAztecStatus();

  lambda_min = status[AZ_lambda_min];
  lambda_max = status[AZ_lambda_max];

  return(0);
#else
  throw std::runtime_error("Tifpack::Chebyshev::CG: support for AztecOO not currently implemented.");
#endif
}

//==========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
std::string Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::description() const {
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }
  //
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols()
      << "}";
  return oss.str();
}

//==========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  const int myImageID = Comm_->getRank();
  Teuchos::OSTab tab(out);

  Scalar MinVal, MaxVal;
  if (IsComputed_) {
    Teuchos::ArrayRCP<const Scalar> DiagView = InvDiagonal_->get1dView();
    Scalar myMinVal = DiagView[0];
    Scalar myMaxVal = DiagView[0];
    for(typename Teuchos::ArrayRCP<Scalar>::size_type i=1; i<DiagView.size(); ++i) {
      if (myMinVal > DiagView[i]) myMinVal = DiagView[i];
      if (myMaxVal < DiagView[i]) myMaxVal = DiagView[i];
    }
    Teuchos::reduceAll(*Comm_, Teuchos::REDUCE_MIN, 1, &myMinVal, &MinVal);
    Teuchos::reduceAll(*Comm_, Teuchos::REDUCE_MAX, 1, &myMaxVal, &MaxVal);
  }

  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium: 
  //    high: 
  // extreme: 
  if (vl != VERB_NONE && myImageID == 0) {
    out << this->description() << endl;
    out << endl;
    out << "===============================================================================" << std::endl;
    out << "Degree of polynomial      = " << PolyDegree_ << std::endl;
    if   (ZeroStartingSolution_) { out << "Using zero starting solution" << endl; }
    else                         { out << "Using input starting solution" << endl; }
    if   (Condest_ == -1.0) { out << "Condition number estimate       = N/A" << endl; }
    else                    { out << "Condition number estimate       = " << Condest_ << endl; }
    if (IsComputed_) {
      out << "Minimum value on stored inverse diagonal = " << MinVal << std::endl;
      out << "Maximum value on stored inverse diagonal = " << MaxVal << std::endl;
    }
    out << std::endl;
    out << "Phase           # calls    Total Time (s)     Total MFlops      MFlops/s       " << endl;
    out << "------------    -------    ---------------    ---------------   ---------------" << endl;
    out << setw(12) << "initialize()" << setw(5) << getNumInitialize() << "    " << setw(15) << getInitializeTime() << endl;
    out << setw(12) << "compute()" << setw(5) << getNumCompute()    << "    " << setw(15) << getComputeTime() << "    " 
        << setw(15) << getComputeFlops() << "    " 
        << setw(15) << (getComputeTime() != 0.0 ? getComputeFlops() / getComputeTime() * 1.0e-6 : 0.0) << endl;
    out << setw(12) << "apply()" << setw(5) << getNumApply()    << "    " << setw(15) << getApplyTime() << "    " 
        << setw(15) << getApplyFlops() << "    " 
        << setw(15) << (getApplyTime() != 0.0 ? getApplyFlops() / getApplyTime() * 1.0e-6 : 0.0) << endl;
    out << "===============================================================================" << std::endl;
    out << endl;
  }
}

}//namespace Tifpack

#endif // TIFPACK_CHEBYSHEV_HPP

