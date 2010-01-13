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

#include <string>

#include "Tifpack_ConfigDefs.hpp"
#include "Tifpack_Condest.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_Parameters.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_RowMatrix.hpp"

namespace Tifpack {

//! Tifpack::Chebyshev: class for preconditioning with Chebyshev polynomials

/*!
  The Tifpack::Chebyshev class enables the construction of preconditioners
  based on Chebyshev polynomials for a Tpetra::CrsMatrix.
  Tifpack::Chebyshev is derived from the Tifpack::Preconditioner class, 
  which is itself derived from Tpetra::Operator.
  Therefore this object can be used as preconditioner everywhere an
  ApplyInverse() method is required in the preconditioning step.

  The class is an adaptation of the routine ML_Cheby in Smoother/ml_smoother.hpp

<P> (04/04/06) Flops are not counted in the routine ApplyInverse()

<P> (04/04/06) The switch to use the transpose matrix is not used in ApplyInverse()

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
template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class Chebyshev : virtual public Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  //@{ \name Constructors/Destructors
  //! Chebyshev constructor with given Tpetra::CrsMatrix.
  /*! Creates an instance of Chebyshev class.
   *
   * \param
   * Matrix - (In) Pointer to the CrsMatrix to precondition.
   */
  Chebyshev(const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix);

  //! Destructor.
  virtual ~Chebyshev() {}

  //@}

  /*! This flag can be used to apply the preconditioner to the transpose of
   * the input operator. 
   * 
   * \return Integer error code, set to 0 if successful.  
   * Set to -1 if this implementation does not support transpose.
    */
  virtual int SetUseTranspose(bool UseTranspose_in)
  {
    UseTranspose_ = UseTranspose_in;
    return(0);
  }

  //@}

  //@{ \name Mathematical functions.

  //! Applies the matrix to an Tpetra::MultiVector.
  /*! 
    \param 
    X - (In) A Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing the result.
    */
  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Tpetra::_MultiVector of dimension NumVectors containing result.
    */
  void applyInverse(
          const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}

  //@{ \name Atribute access functions

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

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getRangeMap() const;

  void Initialize();
  
  bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Returns \c true if the preconditioner has been successfully computed.
  bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Computes the preconditioner.
  void Compute();

  //@}
 
  //@{ \name Miscellaneous

  const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Matrix() const 
  {
    return(*Matrix_);
  }

  //! Computes the condition number estimates and returns the value.
  magnitudeType Condest(
       const Tifpack::CondestType CT = Tifpack::Cheap,
       const LocalOrdinal MaxIters = 1550,
       const magnitudeType Tol = 1e-9,
       Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in = 0);

  //! Returns the condition number estimate, or -1.0 if not computed.
  magnitudeType Condest() const
  {
    return(Condest_);
  }

  //! Sets all the parameters for the preconditioner
  void SetParameters(Teuchos::ParameterList& List);

  //! Prints object to an output stream
  std::ostream& Print(std::ostream & os) const;

  //@}

  //@{ \name Timing and flop count

  //! Returns the number of calls to Initialize().
  int NumInitialize() const
  {
    return(NumInitialize_);
  }

  //! Returns the number of calls to Compute().
  int NumCompute() const
  {
    return(NumCompute_);
  }

  //! Returns the number of calls to ApplyInverse().
  int NumApplyInverse() const
  {
    return(NumApplyInverse_);
  }

  //! Returns the time spent in Initialize().
  double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the time spent in Compute().
  double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the time spent in ApplyInverse().
  double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  //! Returns the number of flops in the initialization phase.
  double InitializeFlops() const
  {
    return(0.0);
  }

  //! Returns the number of flops in the computation phase.
  double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  //! Returns the number of flops for the application of the preconditioner.
  double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

  // @}
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

private:
  
  // @}
  // @{ \name Private methods
  
  //! Copy constructor (PRIVATE, should not be used)
  Chebyshev(const Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>& rhs);
  Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>& operator=(const Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>& rhs);

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
  //! Contains the degree of Chebyshev polynomial.
  int PolyDegree_;
  //! If true, use the tranpose of \c Matrix_.
  bool UseTranspose_;
  //! Contains the estimated condition number
  magnitudeType Condest_;
  //! If true, Compute() also computes the condition number estimate.
  bool ComputeCondest_;
  //! Contains the ratio such that [LambdaMax_ / EigRatio_, LambdaMax_]
  //! is the interval of interest for the Chebyshev polynomial.
  Scalar EigRatio_;
  //! Contains the label of this object.
  std::string Label_;
  //! Contains an approximation to the smallest eigenvalue.
  Scalar LambdaMin_;
  //! Contains an approximation to the largest eigenvalue.
  Scalar LambdaMax_;
  //! Contains the minimum value on the diagonal.
  Scalar MinDiagonalValue_;
  // @}

  // @{ Other data
  //! Number of local rows.
  size_t NumMyRows_;
  //! Number of local nonzeros.
  size_t NumMyNonzeros_;
  //! Number of global rows.
  size_t NumGlobalRows_;
  //! Number of global nonzeros.
  size_t NumGlobalNonzeros_;
  //! Pointers to the matrix to be preconditioned
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix_;
  //! Contains the inverse of diagonal elements of \c Matrix.
  mutable Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > InvDiagonal_;
  //! If \c true, the Operator_ is a Tpetra::RowMatrix.
  bool IsRowMatrix_;
  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  // @}

};


//==============================================================================
// NOTE: This constructor has been introduced because SWIG does not appear
//       to appreciate dynamic_cast. An instruction of type
//       Matrix_ = dynamic_cast<const Tpetra::RowMatrix*> in the
//       other construction does not work in PyTrilinos -- of course
//       it does in any C++ code (for a Tpetra::Operator that is also
//       an Tpetra::RowMatrix).
//
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Chebyshev(const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix) :
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  PolyDegree_(1),
  UseTranspose_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  EigRatio_(30.0),
  Label_(),
  LambdaMin_(0.0),
  LambdaMax_(100.0),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Matrix_(Matrix),
  IsRowMatrix_(true), 
  ZeroStartingSolution_(true)
{
}

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameters(Teuchos::ParameterList& List)
{

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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const
{
  return Matrix_->getDomainMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const
{
  return Matrix_->getRangeMap();
}

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(
       const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode) const
{
  TEST_FOR_EXCEPTION(IsComputed() == false, std::runtime_error,
     "Tifpack::Chebyshev::apply ERROR: IsComputed() must be true prior to calling apply.");

  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Tifpack::Chebyshev::apply ERROR: X.getNumVectors() != Y.getNumVectors().");

  Matrix_->apply(X, Y, mode);
}

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Initialize()
{
  IsInitialized_ = false;

  TEST_FOR_EXCEPTION(Matrix_ == Teuchos::null, std::runtime_error, "Tifpack::Chebyshev::Initialize ERROR: Matrix_ == Teuchos::null");

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Teuchos::Time("Chebyshev") );

  TEST_FOR_EXCEPTION(Matrix().getGlobalNumRows() != Matrix().getGlobalNumCols(), std::runtime_error,
     "Tifpack::Chebyshev::Initialize ERROR: only square matrices are supported");

  NumMyRows_ = Matrix_->getNodeNumRows();
  NumMyNonzeros_ = Matrix_->getNodeNumEntries();
  NumGlobalRows_ = Matrix_->getGlobalNumRows();
  NumGlobalNonzeros_ = Matrix_->getGlobalNumEntries();

  ++NumInitialize_;
  Time_->stop();
  InitializeTime_ += Time_->totalElapsedTime();
  IsInitialized_ = true;
}

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Compute()
{
  if (!IsInitialized()) {
    Initialize();
  }

  Time_->start();

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  TEST_FOR_EXCEPTION(PolyDegree_ <= 0, std::runtime_error,
    "Tifpack::Chebyshev::apply ERROR: PolyDegree_ must be at least one");
  
  if (IsRowMatrix_ && InvDiagonal_ == Teuchos::null)
  {
    InvDiagonal_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Matrix().getRowMap()) );

    Matrix().getLocalDiagCopy(*InvDiagonal_);

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

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::ostream& Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Print(std::ostream & os) const
{

  Scalar MinVal, MaxVal;

  const Teuchos::RCP<const Teuchos::Comm<int> >& tcomm = Matrix_->getRowMap()->getComm();
  if (IsComputed_) {
    Teuchos::ArrayRCP<const Scalar> DiagView = InvDiagonal_->get1dView();
    Scalar myMinVal = DiagView[0];
    Scalar myMaxVal = DiagView[0];
    for(typename Teuchos::ArrayRCP<Scalar>::size_type i=1; i<DiagView.size(); ++i) {
      if (myMinVal > DiagView[i]) myMinVal = DiagView[i];
      if (myMaxVal < DiagView[i]) myMaxVal = DiagView[i];
    }
    Teuchos::reduceAll(*tcomm, Teuchos::REDUCE_MIN, 1, &myMinVal, &MinVal);
    Teuchos::reduceAll(*tcomm, Teuchos::REDUCE_MAX, 1, &myMaxVal, &MaxVal);
  }

  if (tcomm->getRank()==0) {
    os << std::endl;
    os << "================================================================================" << std::endl;
    os << "Chebyshev" << std::endl;
    os << "Degree of polynomial      = " << PolyDegree_ << std::endl;
    os << "Condition number estimate = " << Condest() << std::endl;
    os << "Global number of rows     = " << Matrix_->getRangeMap()->getGlobalNumElements() << std::endl;
    if (IsComputed_) {
      os << "Minimum value on stored inverse diagonal = " << MinVal << std::endl;
      os << "Maximum value on stored inverse diagonal = " << MaxVal << std::endl;
    }
    if (ZeroStartingSolution_) 
      os << "Using zero starting solution" << std::endl;
    else
      os << "Using input starting solution" << std::endl;
    os << std::endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << std::endl;
    os << "-----           -------   --------------       ------------     --------" << std::endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize_ 
       << "  " << std::setw(15) << InitializeTime_ 
       << "              0.0              0.0" << std::endl;
    os << "Compute()       "   << std::setw(5) << NumCompute_ 
       << "  " << std::setw(15) << ComputeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
    if (ComputeTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << std::endl;
    else
      os << "  " << std::setw(15) << 0.0 << std::endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse_ 
       << "  " << std::setw(15) << ApplyInverseTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
    if (ApplyInverseTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << std::endl;
    else
      os << "  " << std::setw(15) << 0.0 << std::endl;
    os << "================================================================================" << std::endl;
    os << std::endl;
  }

  return(os);
}

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::magnitudeType
Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Condest(
       const Tifpack::CondestType CT,
       const LocalOrdinal MaxIters,
       const magnitudeType Tol,
       Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // always compute it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Tifpack::Condest(*this, CT, MaxIters, Tol, Matrix_in);

  return(Condest_);
}

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y, Teuchos::ETransp mode) const
{
  TEST_FOR_EXCEPTION(!IsComputed(), std::runtime_error, "Tifpack::Chebyshev::applyInverse ERROR, not yet computed.");

  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Tifpack::Chebyshev::applyInverse ERROR: X.getNumVectors() != Y.getNumVectors().");

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

  //--- Initialize coefficients
  // Note that delta stores the inverse of ML_Cheby::delta
  Scalar alpha = LambdaMax_ / EigRatio_;
  Scalar beta = 1.1 * LambdaMax_;
  Scalar delta = 2.0 / (beta - alpha);
  Scalar theta = 0.5 * (beta + alpha);
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
    Matrix_->apply(Y, V);
    // Compute W = invDiag * ( X - V )/ Theta
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
    // Compute W = invDiag * X / Theta
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
    Matrix_->apply(Y, V);
    rhokp1 = 1.0 / (2.0*s1 - rhok);
    dtemp1 = rhokp1 * rhok;
    dtemp2 = 2.0 * rhokp1 * delta;
    rhok = rhokp1;
    // Compute W = dtemp1 * W
    W.scale(dtemp1);
    // Compute W = W + dtemp2 * invDiag * ( X - V )
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

  ++NumApplyInverse_;
  Time_->stop();
  ApplyInverseTime_ += Time_->totalElapsedTime();
}

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
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

//==============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Chebyshev<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
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

}//namespace Tifpack

#endif // TIFPACK_CHEBYSHEV_HPP

