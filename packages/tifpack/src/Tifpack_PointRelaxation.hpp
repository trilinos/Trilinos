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
#include "Teuchos_RCP.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Operator.hpp"
#include "Tifpack_Condest.hpp"
#include "Tifpack_Preconditioner.hpp"
#include "Tifpack_Parameters.hpp"

#include "Tpetra_Vector.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_Import.hpp"

#include "Teuchos_ParameterList.hpp"

#include <string>
#include <sstream>
#include <iostream>

namespace Tifpack {
enum PointRelaxationType {
  JACOBI,
  GS,
  SGS
};

//! Tifpack::PointRelaxation: a class to define point relaxation preconditioners for Tpetra::RowMatrix objects.

/*! 
  The Tifpack::PointRelaxation class enables the construction of point
  relaxation preconditioners of a Tpetra::RowMatrix. Tifpack::PointRelaxation 
  is derived from the Tifpack::Preconditioner class, which is itself derived
  from Tpetra::Operator.
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

Using Tifpack::Jacobi, the user can apply the specified number of sweeps
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
Tifpack::GaussSeidel does not consider backward Gauss-Seidel methods.

<P>For a list of supported parameters, please refer to page \ref ifp_params.

<P>The complete list of supported parameters is reported in page \ref ifp_params. For a presentation of basic relaxation schemes, please refer to page
\ref Tifpack::PointRelaxation.

\author Michael Heroux, SNL 9214.

\date Last modified on 22-Jan-05.
  
*/
template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class PointRelaxation : virtual public Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  //@{ \name Constructors/Destructors
  //! Tifpack::PointRelaxation constructor with given Tpetra::RowMatrix.
  /*! Creates an instance of Tifpack::PointRelaxation class.
   *
   * \param
   * Matrix - (In) Pointer to matrix to precondition.
   */
  PointRelaxation(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix);

  //! Destructor.
  virtual ~PointRelaxation() {}

  //@}

  /*! This flag can be used to apply the preconditioner to the transpose of
   * the input operator. 
   * 
   * \return Integer error code, set to 0 if successful.  
   * Set to -1 if this implementation does not support transpose.
    */
  int SetUseTranspose(bool UseTranspose_in)
  {
    UseTranspose_ = UseTranspose_in;
    return(0);
  }

  //@}

  //@{ \name Mathematical functions.

  //! Applies the matrix to a Tpetra::MultiVector.
  /*! 
    \param 
    X - (In) A Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
  void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Tpetra::MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning This routine is NOT AztecOO complaint.
    */
  void applyInverse(
          const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //! Returns the infinity norm of the global matrix (not implemented)
  magnitudeType NormInf() const
  {
    return(-1);
  }
  //@}

  //@{ \name Atribute access functions

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const
  {
    return(UseTranspose_);
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const
  {
    return(false);
  }

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getRangeMap() const;

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

  //! Computes the preconditioners.
  void Compute();

  //@}
 
  //@{ \name Miscellaneous

  const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Matrix() const 
  {
    return(*Matrix_);
  }

  //! Computes the condition number estimate and returns the value.
  magnitudeType Condest(
       const Tifpack::CondestType CT = Tifpack::Cheap,
       const LocalOrdinal MaxIters = 1550,
       const magnitudeType Tol = 1e-9,
			 Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix = 0);

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

private:
 
  // @{ Application of the preconditioner
  
  //! Applies the Jacobi preconditioner to X, returns the result in Y.
  void ApplyInverseJacobi(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  //! Applies the Gauss-Seidel preconditioner to X, returns the result in Y.
  void ApplyInverseGS(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  void ApplyInverseGS_RowMatrix(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  void ApplyInverseGS_CrsMatrix(
        const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  //! Applies the symmetric Gauss-Seidel preconditioner to X, returns the result in Y.
  void ApplyInverseSGS(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  void ApplyInverseSGS_RowMatrix(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  void ApplyInverseSGS_CrsMatrix(
        const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;

  //@}

private:
  
  //! Copy constructor (PRIVATE, should not be used)
  PointRelaxation(const Tifpack::PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>& rhs);
 
  //! operator = (PRIVATE, should not be used)
  PointRelaxation& operator=(const Tifpack::PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>& rhs);

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
  magnitudeType Condest_;
  //! If true, Compute() also computes the condition number estimate.
  bool ComputeCondest_;
  //! Contains the label of this object.
  std::string Label_;
  int PrecType_;
  Scalar MinDiagonalValue_;
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
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix_;
  //! Importer for parallel GS and SGS
  Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > Importer_;
  //! Contains the diagonal elements of \c Matrix.
  mutable Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Diagonal_;
  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  //! Backward-Mode Gauss Seidel 
  bool DoBackwardGS_;
  // @}
};//class PointRelaxation

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::PointRelaxation(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix_in)
: IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  NumSweeps_(1),
  DampingFactor_(1.0),
  UseTranspose_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  PrecType_(Tifpack::JACOBI),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Matrix_(Matrix_in),
  IsParallel_(false),
  ZeroStartingSolution_(true),
  DoBackwardGS_(false)
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::SetParameters(Teuchos::ParameterList& List)
{
  Teuchos::ParameterList validparams;
  Tifpack::GetValidParameters(validparams);
  List.validateParameters(validparams);

  std::string PT;
  if (PrecType_ == Tifpack::JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == Tifpack::GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == Tifpack::SGS)
    PT = "symmetric Gauss-Seidel";

  Tifpack::GetParameter(List, "relaxation: type", PT);

  if (PT == "Jacobi")
    PrecType_ = Tifpack::JACOBI;
  else if (PT == "Gauss-Seidel")
    PrecType_ = Tifpack::GS;
  else if (PT == "symmetric Gauss-Seidel")
    PrecType_ = Tifpack::SGS;
  else {
    std::ostringstream osstr;
    osstr << "Tifpack::PointRelaxation::SetParameters: unsupported parameter-value for 'relaxation: type' (" << PT << ")";
    std::string str = osstr.str();
    throw std::runtime_error(str);
  }

  Tifpack::GetParameter(List, "relaxation: sweeps",NumSweeps_);
  Tifpack::GetParameter(List, "relaxation: damping factor", DampingFactor_);
  Tifpack::GetParameter(List, "relaxation: min diagonal value", MinDiagonalValue_);
  Tifpack::GetParameter(List, "relaxation: zero starting solution", ZeroStartingSolution_);
  Tifpack::GetParameter(List, "relaxation: backward mode",DoBackwardGS_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getDomainMap() const
{
  return Matrix_->getDomainMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::getRangeMap() const
{
  return Matrix_->getRangeMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::apply(
       const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode) const
{
  TEST_FOR_EXCEPTION(IsComputed() == false, std::runtime_error,
     "Tifpack::PointRelaxation::apply ERROR: IsComputed() must be true prior to calling apply.");

  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Tifpack::PointRelaxation::apply ERROR: X.getNumVectors() != Y.getNumVectors().");

  Matrix_->apply(X, Y, mode);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::applyInverse(
          const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                Teuchos::ETransp mode) const
{
  TEST_FOR_EXCEPTION(IsComputed() == false, std::runtime_error,
     "Tifpack::PointRelaxation::applyInverse ERROR: IsComputed() must be true prior to calling applyInverse.");

  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Tifpack::PointRelaxation::applyInverse ERROR: X.getNumVectors() != Y.getNumVectors().");

  Time_->start(true);

  // AztecOO gives X and Y pointing to the same memory location.
  // In that case we need to create an auxiliary vector, Xcopy
  Teuchos::RCP< const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Xcopy;
  if (&(X.get2dView()[0][0]) == &(Y.get2dView()[0][0]))
    Xcopy = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(X) );
  else
    Xcopy = Teuchos::rcp( &X, false );

  if (ZeroStartingSolution_)
    Y.putScalar(0.0);

  // Flops are updated in each of the following.
  switch (PrecType_) {
  case Tifpack::JACOBI:
    ApplyInverseJacobi(*Xcopy,Y);
    break;
  case Tifpack::GS:
    ApplyInverseGS(*Xcopy,Y);
    break;
  case Tifpack::SGS:
    ApplyInverseSGS(*Xcopy,Y);
    break;
  default:
    throw std::runtime_error("Tifpack::PointRelaxation::applyInverse internal logic error.");
  }

  ++NumApplyInverse_;
  Time_->stop();
  ApplyInverseTime_ += Time_->totalElapsedTime();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyInverseJacobi(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  int NumVectors = X.getNumVectors();
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> A_times_Y( Y.getMap(),NumVectors );

  for (int j = 0; j < NumSweeps_ ; j++) {
    apply(Y,A_times_Y);
    A_times_Y.update(1.0,X,-1.0);
    Y.elementWiseMultiply(DampingFactor_, *Diagonal_, A_times_Y, 1.0);
  }

  // Flops:
  // - matrix vector              (2 * NumGlobalNonzeros_)
  // - update                     (2 * NumGlobalRows_)
  // - Multiply:
  //   - DampingFactor            (NumGlobalRows_)
  //   - Diagonal                 (NumGlobalRows_)
  //   - A + B                    (NumGlobalRows_)
  //   - 1.0                      (NumGlobalRows_)
  ApplyInverseFlops_ += NumVectors * (6 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyInverseGS(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* CrsMatrix = dynamic_cast<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>*>(&*Matrix_);
  // try to pick the best option; performances may be improved
  // if several sweeps are used.
  if (CrsMatrix != 0)
  {
    ApplyInverseGS_CrsMatrix(*CrsMatrix, X, Y);
  }
  else {
    ApplyInverseGS_RowMatrix(X, Y);
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyInverseGS_RowMatrix(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  size_t NumVectors = X.getNumVectors();

  size_t maxLength = Matrix().getNodeMaxNumRowEntries();
  Teuchos::Array<LocalOrdinal> Indices(maxLength);
  Teuchos::Array<Scalar> Values(maxLength);

  Teuchos::RCP< Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Y2;
  if (IsParallel_)
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Importer_->getTargetMap(), NumVectors) );
  else
    Y2 = Teuchos::rcp( &Y, false );

  // extract views (for nicer and faster code)
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > y_ptr = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > y2_ptr = Y2->get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr =  X.get2dView();
  Teuchos::ArrayRCP<const Scalar> d_ptr = Diagonal_->get1dView();
 
  for (int j = 0; j < NumSweeps_ ; j++) {

    // data exchange is here, once per sweep
    if (IsParallel_)
      Y2->doImport(Y,*Importer_,Tpetra::INSERT);

    if(!DoBackwardGS_){      
      /* Forward Mode */
      for (int i = 0 ; i < NumMyRows_ ; ++i) {
        
        size_t NumEntries;
        Matrix_->getLocalRowCopy(i, Indices(), Values(), NumEntries);
        
        for (size_t m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;
          for (size_t k = 0 ; k < NumEntries ; ++k) {
            LocalOrdinal col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }
          
          y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
        }
      }
    }
    else {
      /* Backward Mode */
      for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {
        size_t NumEntries;
        Matrix_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

        for (size_t m = 0 ; m < NumVectors ; ++m) {
          
          Scalar dtemp = 0.0;
          for (size_t k = 0 ; k < NumEntries ; ++k) {
            LocalOrdinal col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
        }
      }
    }

    // using Export() sounded quite expensive   
    if (IsParallel_) {
      for (size_t m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
    }
  }

  ApplyInverseFlops_ += NumVectors * (4 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyInverseGS_CrsMatrix(
        const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  size_t NumVectors = X.getNumVectors();

  Teuchos::ArrayRCP<const LocalOrdinal> Indices;
  Teuchos::ArrayRCP<const Scalar> Values;

  Teuchos::RCP< Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Importer_->getTargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > y_ptr = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > y2_ptr = Y2->get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr =  X.get2dView();
  Teuchos::ArrayRCP<const Scalar> d_ptr = Diagonal_->get1dView();
  
  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      Y2->doImport(Y,*Importer_,Tpetra::INSERT);

    if(!DoBackwardGS_){  
      /* Forward Mode */
      for (int i = 0 ; i < NumMyRows_ ; ++i) {

        LocalOrdinal col;
        Scalar diag = d_ptr[i];
        
        A.getLocalRowView(i, Indices, Values);
        size_t NumEntries = Indices.size();
        
        for (size_t m = 0 ; m < NumVectors ; ++m) {
          
          Scalar dtemp = 0.0;
          
          for (size_t k = 0; k < NumEntries; ++k) {
            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }
          
          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
        }
      }
    }
    else {
      /* Backward Mode */
      for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

        LocalOrdinal col;
        Scalar diag = d_ptr[i];
        
        A.getLocalRowView(i, Indices, Values);
        size_t NumEntries = Indices.size();
        
        for (size_t m = 0 ; m < NumVectors ; ++m) {
          
          Scalar dtemp = 0.0;
          for (size_t k = 0; k < NumEntries; ++k) {
            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }
          
          y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
          
        }
      }
    }
    
    if (IsParallel_) {
      for (size_t m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
    }
  }

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyInverseSGS(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* CrsMatrix = dynamic_cast<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>*>(&*Matrix_);
  // try to pick the best option; performance may be improved
  // if several sweeps are used.
  if (CrsMatrix != 0)
  {
    ApplyInverseSGS_CrsMatrix(*CrsMatrix, X, Y);
  }
  else
    ApplyInverseSGS_RowMatrix(X, Y);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyInverseSGS_RowMatrix(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  size_t NumVectors = X.getNumVectors();
  size_t maxLength = Matrix().getNodeMaxNumRowEntries();
  Teuchos::Array<LocalOrdinal> Indices(maxLength);
  Teuchos::Array<Scalar> Values(maxLength);

  Teuchos::RCP< Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Importer_->getTargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > y_ptr = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > y2_ptr = Y2->get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr =  X.get2dView();
  Teuchos::ArrayRCP<const Scalar> d_ptr = Diagonal_->get1dView();
  
  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      Y2->doImport(Y,*Importer_,Tpetra::INSERT);

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      size_t NumEntries;
      Scalar diag = d_ptr[i];

      Matrix_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

      for (size_t m = 0 ; m < NumVectors ; ++m) {

        Scalar dtemp = 0.0;

        for (size_t k = 0 ; k < NumEntries ; ++k) {
          LocalOrdinal col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      size_t NumEntries;
      Scalar diag = d_ptr[i];

      Matrix_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

      for (size_t m = 0 ; m < NumVectors ; ++m) {

        Scalar dtemp = 0.0;
        for (size_t k = 0 ; k < NumEntries ; ++k) {
          LocalOrdinal col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    if (IsParallel_) {
      for (size_t m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
    }
  }

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::ApplyInverseSGS_CrsMatrix(
        const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& A,
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  size_t NumVectors = X.getNumVectors();

  Teuchos::ArrayRCP<const LocalOrdinal> Indices;
  Teuchos::ArrayRCP<const Scalar> Values;

  Teuchos::RCP< Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Y2;
  if (IsParallel_) {
    Y2 = Teuchos::rcp( new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Importer_->getTargetMap(), NumVectors) );
  }
  else
    Y2 = Teuchos::rcp( &Y, false );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > y_ptr = Y.get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > y2_ptr = Y2->get2dViewNonConst();
  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > x_ptr =  X.get2dView();
  Teuchos::ArrayRCP<const Scalar> d_ptr = Diagonal_->get1dView();
  
  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      Y2->doImport(Y,*Importer_,Tpetra::INSERT);

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      Scalar diag = d_ptr[i];

      A.getLocalRowView(i, Indices, Values);
      size_t NumEntries = Indices.size();

      for (size_t m = 0 ; m < NumVectors ; ++m) {

        Scalar dtemp = 0.0;

        for (size_t k = 0; k < NumEntries; ++k) {

          LocalOrdinal col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      Scalar diag = d_ptr[i];

      A.getLocalRowView(i, Indices, Values);
      size_t NumEntries = Indices.size();

      for (size_t m = 0 ; m < NumVectors ; ++m) {

        Scalar dtemp = 0.0;
        for (size_t k = 0; k < NumEntries; ++k) {

          LocalOrdinal col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) * diag;

      }
    }

    if (IsParallel_) {
      for (size_t m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
    }
  }

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
typename PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::magnitudeType
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Condest(
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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
std::ostream& PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Print(std::ostream & os) const
{
  Scalar MinVal = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar MaxVal = Teuchos::ScalarTraits<Scalar>::zero();

  const Teuchos::RCP<const Teuchos::Comm<int> >& tcomm = Matrix_->getRowMap()->getComm();
  if (IsComputed_) {
    Teuchos::ArrayRCP<Scalar> DiagView = Diagonal_->get1dViewNonConst();
    Scalar myMinVal = DiagView[0];
    Scalar myMaxVal = DiagView[0];
    for(typename Teuchos::ArrayRCP<Scalar>::size_type i=0; i<DiagView.size(); ++i) {
      if (myMinVal > DiagView[i]) myMinVal = DiagView[i];
      if (myMaxVal < DiagView[i]) myMaxVal = DiagView[i];
    }

    Teuchos::reduceAll(*tcomm, Teuchos::REDUCE_MIN, 1, &myMinVal, &MinVal);
    Teuchos::reduceAll(*tcomm, Teuchos::REDUCE_MAX, 1, &myMaxVal, &MaxVal);
  }

  if (!tcomm->getRank()) {
    os << std::endl;
    os << "================================================================================" << std::endl;
    os << "Tifpack_PointRelaxation" << std::endl;
    os << "Sweeps         = " << NumSweeps_ << std::endl;
    os << "damping factor = " << DampingFactor_ << std::endl;
    if (PrecType_ == Tifpack::JACOBI)
      os << "Type           = Jacobi" << std::endl;
    else if (PrecType_ == Tifpack::GS)
      os << "Type           = Gauss-Seidel" << std::endl;
    else if (PrecType_ == Tifpack::SGS)
      os << "Type           = symmetric Gauss-Seidel" << std::endl;
    if (DoBackwardGS_) 
      os << "Using backward mode (GS only)" << std::endl;
    if (ZeroStartingSolution_) 
      os << "Using zero starting solution" << std::endl;
    else
      os << "Using input starting solution" << std::endl;
    os << "Condition number estimate = " << Condest_ << std::endl;
    os << "Global number of rows            = " << Matrix_->getGlobalNumRows() << std::endl;
    if (IsComputed_) {
      os << "Minimum value on stored diagonal = " << MinVal << std::endl;
      os << "Maximum value on stored diagonal = " << MaxVal << std::endl;
    }
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

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Initialize()
{
  IsInitialized_ = false;

  TEST_FOR_EXCEPTION(Matrix_ == Teuchos::null, std::runtime_error,
    "Tifpack::PointRelaxation::Initialize ERROR, Matrix is NULL");

  if (Time_ == Teuchos::null)
    Time_ = Teuchos::rcp( new Teuchos::Time("Tifpack::PointRelaxation") );

//  size_t globalrows = Matrix().getGlobalNumRows();
//  size_t globalcols = Matrix().getGlobalNumCols();
//  TEST_FOR_EXCEPTION(globalrows != globalcols, std::runtime_error,
//   "Tifpack::PointRelaxation::Initialize ERROR, only square matrices are supported");

  NumMyRows_ = Matrix_->getNodeNumRows();
  NumMyNonzeros_ = Matrix_->getNodeNumEntries();
  NumGlobalRows_ = Matrix_->getGlobalNumRows();
  NumGlobalNonzeros_ = Matrix_->getGlobalNumEntries();

  const Teuchos::RCP<const Teuchos::Comm<int> >& tcomm = Matrix_->getRowMap()->getComm();
  if (tcomm->getSize() != 1)
    IsParallel_ = true;
  else
    IsParallel_ = false;

  ++NumInitialize_;
  InitializeTime_ += Time_->totalElapsedTime();
  IsInitialized_ = true;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Compute()
{
  if (!IsInitialized())
    Initialize();

  Time_->start(true);

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  TEST_FOR_EXCEPTION(NumSweeps_ < 0, std::runtime_error,
    "Tifpack::PointRelaxation::Compute, NumSweeps_ must be >= 0");
  
  Diagonal_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Matrix().getRowMap()) );

  TEST_FOR_EXCEPTION(Diagonal_ == Teuchos::null, std::runtime_error,
    "Tifpack::PointRelaxation::Compute, failed to create Diagonal_");

  Matrix().getLocalDiagCopy(*Diagonal_);

  // check diagonal elements, store the inverses, and verify that
  // no zeros are around. If an element is zero, then by default
  // its inverse is zero as well (that is, the row is ignored).
  Teuchos::ArrayRCP<Scalar> DiagView = Diagonal_->get1dViewNonConst();
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    Scalar& diag = DiagView[i];
    if (std::abs(diag) < MinDiagonalValue_)
      diag = MinDiagonalValue_;
    if (diag != 0.0)
      diag = 1.0 / diag;
  }
  ComputeFlops_ += NumMyRows_;

  //Marzio's comment:
  // We need to import data from external processors. Here I create a
  // Tpetra::Import object because I cannot assume that Matrix_ has one.
  // This is a bit of waste of resources (but the code is more robust).
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && ((PrecType_ == Tifpack::GS) || (PrecType_ == Tifpack::SGS))) {
    Importer_ = Teuchos::rcp( new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(Matrix().getColMap(),
                                  Matrix().getRowMap()) );
    TEST_FOR_EXCEPTION(Importer_ == Teuchos::null, std::runtime_error,
      "Tifpack::PointRelaxation::Compute ERROR failed to create Importer_");
  }

  ++NumCompute_;
  Time_->stop();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

}//namespace Tifpack

#endif // TIFPACK_POINTRELAXATION_HPP

