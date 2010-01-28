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
#include "Tifpack_Condest.hpp"
#include "Tifpack_Parameters.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <Teuchos_TestForException.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <string>
#include <iostream>
#include <sstream>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

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
template<class Scalar,class LocalOrdinal = int,class GlobalOrdinal = LocalOrdinal,class Node = Kokkos::DefaultNode::DefaultNodeType,class LocalMatVec = Kokkos::DefaultSparseMultiply<Scalar,LocalOrdinal,Node>,class LocalMatSolve = Kokkos::DefaultSparseSolve<Scalar,LocalOrdinal,Node> >
class PointRelaxation : virtual public Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

public:
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  // \name Constructors and Destructors
  //@{

  //! PointRelaxation constructor with given Tpetra::RowMatrix input.
  explicit PointRelaxation(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix);

  //! PointRelaxation Destructor.
  virtual ~PointRelaxation();

  //@}
  //@{ Construction methods
  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes five parameter names: level_fill, drop_tolerance,
     absolute_threshold, relative_threshold and overlap_mode. These names are
     case insensitive. For level_fill the ParameterEntry must have type int, the 
     threshold entries must have type double and overlap_mode must have type
     Tpetra::CombineMode.
  */
  void setParameters(Teuchos::ParameterList& params);

  //! Initialize L and U with values from user matrix A.
  /*! Copies values from the user's matrix into the nonzero pattern of L and U.
    \param In 
           A - User matrix to be factored.
    \warning The graph of A must be identical to the graph passed in to IlukGraph constructor.
             
   */
  void initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(IsInitialized_);
  }

  //! Compute IC factor U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the RILU(k) factors L and U using the current:
    <ol>
    <li> IlukGraph specifying the structure of L and U.
    <li> Value for the RILU(k) relaxation parameter.
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
    InitValues() must be called before the factorization can proceed.
   */
  void compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  inline bool isComputed() const {
    return(IsComputed_);
  }

  //@}

  //! @name Methods implementing Tpetra::Operator.
  //@{ 

  //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y - (InOut) A Tpetra::MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning This routine is NOT AztecOO compliant.
  */
  void apply(
      const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
            Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getDomainMap() const;

  //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& getRangeMap() const;

  bool hasTransposeApply() const;

  //@}

  //@{
  //! \name Mathematical functions.

  //! Applies the matrix to a Tpetra::MultiVector.
  /*! 
    \param 
    X - (In) A Tpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
  void applyMat(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //! Computes the estimated condition number and returns the value.
  magnitudeType computeCondEst(CondestType CT = Cheap, 
                               LocalOrdinal MaxIters = 1550,
                               magnitudeType Tol = 1e-9,
                               const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &matrix = Teuchos::null);

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

  //! Returns the number of calls to Compute().
  int getNumCompute() const;

  //! Returns the number of calls to apply().
  int getNumApply() const;

  //! Returns the time spent in initialize().
  double getInitializeTime() const;

  //! Returns the time spent in Compute().
  double getComputeTime() const;

  //! Returns the time spent in apply().
  double getApplyTime() const;


  // @}

  //! @name Overridden from Teuchos::Describable 
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:

  // @{ Internal methods

  //! Copy constructor (should never be used)
  PointRelaxation(const PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>& RHS);

  //! operator= (should never be used)
  PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>& operator=(const PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>& RHS);

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

  // @{ Internal data and parameters

  //! reference to the matrix to be preconditioned.
  const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! Reference to the communicator object.
  const Teuchos::RCP<const Teuchos::Comm<int> > Comm_;
  //! Importer for parallel GS and SGS
  Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > Importer_;
  //! Contains the diagonal elements of \c Matrix.
  mutable Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Diagonal_;
  //! Number of application of the preconditioner (should be greater than 0).
  int NumSweeps_;
  //! Which type of point relaxation approach to use
  int PrecType_;
  //! Minimum diagonal value
  Scalar MinDiagonalValue_;
  //! Damping factor.
  double DampingFactor_;
  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! If \c true, more than 1 processor is currently used.
  bool IsParallel_;
  //! If \c true, the starting solution is always the zero vector.
  bool ZeroStartingSolution_;
  //! Backward-Mode Gauss Seidel 
  bool DoBackwardGS_;
  //! Condition number estimate
  magnitudeType Condest_;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been computed successfully.
  bool IsComputed_;
  //! Contains the number of successful calls to initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to apply().
  mutable int NumApply_;
  //! Contains the time for all successful calls to initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to apply().
  mutable double ApplyTime_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyFlops_;
  //! Number of local rows.
  LocalOrdinal NumMyRows_;
  //! Number of global rows.
  global_size_t NumGlobalRows_;
  //! Number of global nonzeros.
  global_size_t NumGlobalNonzeros_;

  //@}

}; //class PointRelaxation


//==============================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::PointRelaxation(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A)
: A_(A),
  Comm_(A->getRowMap()->getComm()),
  NumSweeps_(1),
  PrecType_(Tifpack::JACOBI),
  MinDiagonalValue_(0.0),
  DampingFactor_(1.0),
  Time_( Teuchos::rcp( new Teuchos::Time("Tifpack::PointRelaxation") ) ),
  IsParallel_(false),
  ZeroStartingSolution_(true),
  DoBackwardGS_(false),
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
      Teuchos::typeName(*this) << "::PointRelaxation(): input matrix reference was null.");
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~PointRelaxation() {
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::setParameters(Teuchos::ParameterList& List)
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
    osstr << "Tifpack::PointRelaxation::setParameters: unsupported parameter-value for 'relaxation: type' (" << PT << ")";
    std::string str = osstr.str();
    throw std::runtime_error(str);
  }

  Tifpack::GetParameter(List, "relaxation: sweeps",NumSweeps_);
  Tifpack::GetParameter(List, "relaxation: damping factor", DampingFactor_);
  Tifpack::GetParameter(List, "relaxation: min diagonal value", MinDiagonalValue_);
  Tifpack::GetParameter(List, "relaxation: zero starting solution", ZeroStartingSolution_);
  Tifpack::GetParameter(List, "relaxation: backward mode",DoBackwardGS_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
const Teuchos::RCP<const Teuchos::Comm<int> > & 
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComm() const{
  return(Comm_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getMatrix() const {
  return(A_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getDomainMap() const
{
  return A_->getDomainMap();
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >&
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRangeMap() const
{
  return A_->getRangeMap();
}

//==============================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
bool PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApply() const {
  return true;
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
int PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
int PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
int PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getApplyTime() const {
  return(ApplyTime_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getComputeFlops() const {
  return(ComputeFlops_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
double PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getApplyFlops() const {
  return(ApplyFlops_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getCondEst() const {
  return(Condest_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
typename PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::magnitudeType
PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::computeCondEst(
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
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(
          const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                Teuchos::ETransp mode) const {
  TEST_FOR_EXCEPTION(isComputed() == false, std::runtime_error,
     "Tifpack::PointRelaxation::applyInverse ERROR: isComputed() must be true prior to calling applyInverse.");

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

  ++NumApply_;
  Time_->stop();
  ApplyTime_ += Time_->totalElapsedTime();
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::applyMat(
       const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
             Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
             Teuchos::ETransp mode) const
{
  TEST_FOR_EXCEPTION(isComputed() == false, std::runtime_error,
     "Tifpack::PointRelaxation::apply ERROR: isComputed() must be true prior to calling apply.");

  TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
     "Tifpack::PointRelaxation::apply ERROR: X.getNumVectors() != Y.getNumVectors().");

  A_->apply(X, Y, mode);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::initialize() {
  IsInitialized_ = false;

  TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error,
    "Tifpack::PointRelaxation::Initialize ERROR, Matrix is NULL");

//  size_t globalrows = A_->getGlobalNumRows();
//  size_t globalcols = A_->getGlobalNumCols();
//  TEST_FOR_EXCEPTION(globalrows != globalcols, std::runtime_error,
//   "Tifpack::PointRelaxation::Initialize ERROR, only square matrices are supported");

  NumMyRows_ = A_->getNodeNumRows();
  NumGlobalRows_ = A_->getGlobalNumRows();
  NumGlobalNonzeros_ = A_->getGlobalNumEntries();

  if (Comm_->getSize() != 1)
    IsParallel_ = true;
  else
    IsParallel_ = false;

  ++NumInitialize_;
  InitializeTime_ += Time_->totalElapsedTime();
  IsInitialized_ = true;
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::compute()
{
  if (!isInitialized())
    initialize();

  Time_->start(true);

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  TEST_FOR_EXCEPTION(NumSweeps_ < 0, std::runtime_error,
    "Tifpack::PointRelaxation::compute, NumSweeps_ must be >= 0");
  
  Diagonal_ = Teuchos::rcp( new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(A_->getRowMap()) );

  TEST_FOR_EXCEPTION(Diagonal_ == Teuchos::null, std::runtime_error,
    "Tifpack::PointRelaxation::compute, failed to create Diagonal_");

  A_->getLocalDiagCopy(*Diagonal_);

  // check diagonal elements, store the inverses, and verify that
  // no zeros are around. If an element is zero, then by default
  // its inverse is zero as well (that is, the row is ignored).
  Teuchos::ArrayRCP<Scalar> DiagView = Diagonal_->get1dViewNonConst();
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    Scalar& diag = DiagView[i];
    if (Teuchos::ScalarTraits<Scalar>::magnitude(diag) < MinDiagonalValue_)
      diag = MinDiagonalValue_;
    if (diag != Teuchos::ScalarTraits<Scalar>::zero())
      diag = Teuchos::ScalarTraits<Scalar>::one() / diag;
  }
  ComputeFlops_ += NumMyRows_;

  //Marzio's comment:
  // We need to import data from external processors. Here I create a
  // Tpetra::Import object because I cannot assume that A_ has one.
  // This is a bit of waste of resources (but the code is more robust).
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && ((PrecType_ == Tifpack::GS) || (PrecType_ == Tifpack::SGS))) {
    Importer_ = Teuchos::rcp( new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(A_->getColMap(),
                                  A_->getRowMap()) );
    TEST_FOR_EXCEPTION(Importer_ == Teuchos::null, std::runtime_error,
      "Tifpack::PointRelaxation::compute ERROR failed to create Importer_");
  }

  ++NumCompute_;
  Time_->stop();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::ApplyInverseJacobi(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  int NumVectors = X.getNumVectors();
  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> A_times_Y( Y.getMap(),NumVectors );

  for (int j = 0; j < NumSweeps_ ; j++) {
    applyMat(Y,A_times_Y);
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
  ApplyFlops_ += NumVectors * (6 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::ApplyInverseGS(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* CrsMatrix = dynamic_cast<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>*>(&*A_);
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

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::ApplyInverseGS_RowMatrix(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  size_t NumVectors = X.getNumVectors();

  size_t maxLength = A_->getNodeMaxNumRowEntries();
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
        A_->getLocalRowCopy(i, Indices(), Values(), NumEntries);
        
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
    else {
      /* Backward Mode */
      for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {
        size_t NumEntries;
        A_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

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

  ApplyFlops_ += NumVectors * (4 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::ApplyInverseGS_CrsMatrix(
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

  ApplyFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::ApplyInverseSGS(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* CrsMatrix = dynamic_cast<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>*>(&*A_);
  // try to pick the best option; performance may be improved
  // if several sweeps are used.
  if (CrsMatrix != 0)
  {
    ApplyInverseSGS_CrsMatrix(*CrsMatrix, X, Y);
  }
  else
    ApplyInverseSGS_RowMatrix(X, Y);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::ApplyInverseSGS_RowMatrix(
        const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
              Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  size_t NumVectors = X.getNumVectors();
  size_t maxLength = A_->getNodeMaxNumRowEntries();
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

      A_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

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

      A_->getLocalRowCopy(i, Indices(), Values(), NumEntries);

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

  ApplyFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
}

//==========================================================================
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::ApplyInverseSGS_CrsMatrix(
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

  ApplyFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
}

//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
std::string PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::description() const {
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
  if (PrecType_ == Tifpack::JACOBI)   oss << "Type = Jacobi, " << std::endl;
  else if (PrecType_ == Tifpack::GS)  oss << "Type = Gauss-Seidel, " << std::endl;
  else if (PrecType_ == Tifpack::SGS) oss << "Type = Sym. Gauss-Seidel, " << std::endl;
  //
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols()
      << "}";
  return oss.str();
}

//==========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void PointRelaxation<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
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

  Scalar MinVal = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar MaxVal = Teuchos::ScalarTraits<Scalar>::zero();

  if (IsComputed_) {
    Teuchos::ArrayRCP<Scalar> DiagView = Diagonal_->get1dViewNonConst();
    Scalar myMinVal = DiagView[0];
    Scalar myMaxVal = DiagView[0];
    for(typename Teuchos::ArrayRCP<Scalar>::size_type i=0; i<DiagView.size(); ++i) {
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
    out << "===============================================================================" << endl;
    out << "Sweeps         = " << NumSweeps_ << endl;
    out << "damping factor = " << DampingFactor_ << endl;
    if (PrecType_ == Tifpack::GS && DoBackwardGS_) {
      out << "Using backward mode (GS only)" << endl;
    }
    if (ZeroStartingSolution_) {
      out << "Using zero starting solution" << endl;
    }
    else {
      out << "Using input starting solution" << endl;
    }
    if (Condest_ == -1.0) {
      out << "Condition number estimate       = N/A" << endl;
    }
    else {
      out << "Condition number estimate       = " << Condest_ << endl;
    }
    if (IsComputed_) {
      out << "Minimum value on stored diagonal = " << MinVal << endl;
      out << "Maximum value on stored diagonal = " << MaxVal << endl;
    }
    out << endl;
    out << "Phase           # calls    Total Time (s)     Total MFlops      MFlops/s       " << endl;
    out << "------------    -------    ---------------    ---------------   ---------------" << endl;
    out << setw(12) << "initialize()" << setw(5) << getNumInitialize() << "    " << setw(15) << getInitializeTime() << endl;
    out << setw(12) << "compute()" << setw(5) << getNumCompute()    << "    " << setw(15) << getComputeTime() << "    " 
        << setw(15) << getComputeFlops() << "    " 
        << setw(15) << (getComputeTime() != 0.0 ? getComputeFlops() / getComputeTime() * 1.0e-6 : 0.0) << endl;
    out << setw(12) << "apply()" << setw(5) << getNumApply()    << "    " << setw(15) << getApplyTime() << "    " 
        << setw(15) << getApplyFlops() << "    " 
        << setw(15) << (getApplyTime() != 0.0 ? getApplyFlops() / getApplyTime() * 1.0e-6 : 0.0) << endl;
    out << "===============================================================================" << endl;
    out << endl;
  }
}

}//namespace Tifpack

#endif // TIFPACK_POINTRELAXATION_HPP

