/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef IFPACK2_SPARSECONTAINER_DEF_HPP
#define IFPACK2_SPARSECONTAINER_DEF_HPP

#include "Ifpack2_SparseContainer_decl.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#include "Teuchos_DefaultMpiComm.hpp"
#else
#include "Teuchos_DefaultSerialComm.hpp"
#endif

namespace Ifpack2 {

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::
SparseContainer (const size_t NumRows, const size_t NumVectors) :
  numRows_ (NumRows),
  NumVectors_ (NumVectors),
  IsInitialized_ (false),
  IsComputed_ (false),
#ifdef HAVE_MPI
  LocalComm_ (Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF))),
#else
  LocalComm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ())),
#endif // HAVE_MPI
  needPermutation_ (true)
{}

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::~SparseContainer()
{
  destroy();
}

//==============================================================================
template<class MatrixType, class InverseType>
size_t SparseContainer<MatrixType,InverseType>::getNumRows() const
{
  if (isInitialized()) return numRows_;
  else return 0;
}

//==============================================================================
// Sets the number of vectors for X/Y
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::setNumVectors(const size_t NumVectors_in)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    NumVectors_in <= 0, std::runtime_error, "Ifpack2::SparseContainer::"
    "setNumVectors: The input argument must be positive, but you specified "
    "NumVectors_in = " << NumVectors_in << " <= 0.");

  if (! IsInitialized_  || NumVectors_ != NumVectors_in) {
    NumVectors_=NumVectors_in;
  }
}

//==============================================================================
// Returns the ID associated to local row i.
template<class MatrixType, class InverseType>
typename MatrixType::local_ordinal_type & SparseContainer<MatrixType,InverseType>::ID(const size_t i)
{
  return GID_[i];
}

//==============================================================================
// Returns  true is the container has been successfully initialized.
template<class MatrixType, class InverseType>
bool SparseContainer<MatrixType,InverseType>::isInitialized() const
{
  return IsInitialized_;
}

//==============================================================================
// Returns  true is the container has been successfully computed.
template<class MatrixType, class InverseType>
bool SparseContainer<MatrixType,InverseType>::isComputed() const
{
  return IsComputed_;
}

//==============================================================================
// Sets all necessary parameters.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::setParameters(const Teuchos::ParameterList& List)
{
  List_ = List;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::initialize()
{
  if(IsInitialized_) destroy();
  IsInitialized_=false;

  // Allocations
  Map_ = Teuchos::rcp( new Tpetra::Map<InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>(numRows_,0,LocalComm_) );
  Matrix_ = Teuchos::rcp( new Tpetra::CrsMatrix<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>(Map_,0) );
  GID_.resize(numRows_);
  setNumVectors(NumVectors_);

  // create the inverse
  Inverse_ = Teuchos::rcp( new InverseType(Matrix_) );
  TEUCHOS_TEST_FOR_EXCEPTION( Inverse_ == Teuchos::null, std::runtime_error, "Ifpack2::SparseContainer::initialize error in inverse constructor.");
  Inverse_->setParameters(List_);

  // Note from IFPACK: Call Inverse_->Initialize() in Compute(). This saves
  // some time, because the diagonal blocks can be extracted faster,
  // and only once.
  IsInitialized_ = true;
}

//==============================================================================
// Finalizes the linear system matrix and prepares for the application of the inverse.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
compute (const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix)
{
  IsComputed_=false;
  this->initialize ();

  // extract the submatrices
  this->extract (Matrix);

  // initialize & compute the inverse operator
  Inverse_->initialize ();
  Inverse_->compute ();

  // Determine whether apply() and weightedApply() need to permute X and Y.
  // If ID(i) == i for all i in [0, numRows-1], then we don't need to permute.
  // bool needPermutation = false;
  // const size_t numRows = Inverse_->getRangeMap ()->getNodeNumElements ();
  // if (numRows == this->GID_.size ()) {
  //   for (size_t i = 0; i < numRows; ++i) {
  //     if (this->ID(i) != i) {
  //       needPermutation = true;
  //       break;
  //     }
  //   }
  // }
  // needPermutation_ = needPermutation;
  needPermutation_ = true;

  IsComputed_ = true;
}

//==============================================================================
// Computes Y = alpha * M^{-1} X + beta*Y
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
apply (const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
       Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
       Teuchos::ETransp mode,
       MatrixScalar alpha,
       MatrixScalar beta)
{
  using Teuchos::ArrayRCP;
  using Teuchos::Range1D;
  using Teuchos::RCP;
  typedef Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> MV;
  typedef Tpetra::Map<MatrixLocalOrdinal, MatrixGlobalOrdinal, MatrixNode> map_type;
  const size_t numVecs = X.getNumVectors ();

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, "Ifpack2::SparseContainer::apply: "
    "You must have called the compute() method before you may call apply().  "
    "You may call the apply() method as many times as you want after calling "
    "compute() once, but you must have called compute() at least once.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVecs != Y.getNumVectors (), std::runtime_error,
    "Ifpack2::SparseContainer::apply: X and Y have different numbers of "
    "vectors.  X has " << X.getNumVectors ()
    << ", but Y has " << X.getNumVectors () << ".");

  if (numVecs == 0) {
    return; // done! nothing to do
  }

  // The operator Inverse_ works on a permuted subset of the local
  // parts of X and Y.  The subset and permutation are defined by the
  // GID_ array, which the ID(j) member function accesses.  If the
  // permutation is trivial and the subset is exactly equal to the
  // local indices, then we can use the local parts of X and Y
  // exactly, without needing to permute.  Otherwise, we have to use
  // temporary storage to permute X and Y.

  // Get Maps that describe the local parts of X and Y.  "Local" means
  // "the part owned by the calling MPI process only."  These Maps are
  // trivially "globally distributed" over MPI_COMM_SELF.  While this
  // does technically mean the same thing as "locally replicated" for
  // MPI_COMM_SELF, we say "globally distributed" because this more
  // generally includes communicators other than MPI_COMM_SELF, in
  // case we wish to expand that definition in the future.
  RCP<const map_type> map_X_local =
    rcp (new map_type (X.getLocalLength (), X.getMap ()->getIndexBase (),
                       LocalComm_, Tpetra::GloballyDistributed,
                       X.getMap ()->getNode ()));
  RCP<const map_type> map_Y_local =
    rcp (new map_type (Y.getLocalLength (), Y.getMap ()->getIndexBase (),
                       LocalComm_, Tpetra::GloballyDistributed,
                       Y.getMap ()->getNode ()));

  // Get views of the local parts of X and Y.
  RCP<const MV> X_local = X.offsetView (map_X_local, 0);
  RCP<MV> Y_local = Y.offsetViewNonConst (map_Y_local, 0);

  if (! needPermutation_) {
    Inverse_->apply (*X_local, *Y_local, mode, alpha, beta);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      this->GID_.size () != numRows_, std::logic_error, "this->GID_.size() = "
      << this->GID_.size () << " != numRows_ = " << numRows_ << ".  "
      "Please report this bug to the Ifpack2 developers.");

    // Create temporary permuted versions of the input and output.
    // (Re)allocate X_ and/or Y_ only if necessary.  We'll use them to
    // store the permuted versions of X resp. Y.  The current Maps of
    // X_ and Y_ don't matter, since we're treating them as local MVs.
    // Note that X_local_perm has the domain Map of the operator,
    // which may be a permuted subset of X_local's Map.  Similarly,
    // Y_local_perm has the range Map of the operator, which may be a
    // permuted subset of Y_local's Map.
    //
    // numRows_ here gives the number of rows in the row Map of the
    // local Inverse_ operator.
    //
    // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
    // here that the row Map and the range Map of that operator are
    // the same.
    //
    // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
    // really belongs in Tpetra.

    if (X_.is_null () || X_->getLocalLength () != numRows_ ||
        X_->getNumVectors () < numVecs) {
      X_ = rcp (new MV (Inverse_->getDomainMap (), numVecs));
    }
    RCP<MV> X_local_perm = X_->subViewNonConst (Range1D (0, numVecs - 1));
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const MatrixScalar> X_local_j =
        X_local->getVector (j)->get1dView ();
      ArrayRCP<MatrixScalar> X_local_perm_j =
        X_local_perm->getVectorNonConst (j)->get1dViewNonConst ();
      for (size_t i = 0; i < numRows_; ++i) {
        const size_t i_perm = this->ID (i);
        X_local_perm_j[i] = X_local_j[i_perm];
      }
    }

    // We must permute the output multivector even on input to
    // Inverse_->apply(), since the Inverse_ operator might use it as
    // an initial guess for a linear solve.  We have no way of knowing
    // whether it does or does not.
    if (Y_.is_null () || Y_->getLocalLength () != numRows_ ||
        Y_->getNumVectors () < numVecs) {
      Y_ = rcp (new MV (Inverse_->getRangeMap (), numVecs));
    }
    RCP<MV> Y_local_perm = Y_->subViewNonConst (Range1D (0, numVecs - 1));
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const MatrixScalar> Y_local_j =
        Y_local->getVector (j)->get1dView ();
      ArrayRCP<MatrixScalar> Y_local_perm_j =
        Y_local_perm->getVectorNonConst (j)->get1dViewNonConst ();
      for (size_t i = 0; i < numRows_; ++i) {
        const size_t i_perm = this->ID (i);
        Y_local_perm_j[i] = Y_local_j[i_perm];
      }
    }

    // Apply the local Inverse_ operator to the permuted input and output.
    Inverse_->apply (*X_local_perm, *Y_local_perm, mode, alpha, beta);

    // Copy the permuted subset output vector Y_local_perm into the
    // original output multivector (view) Y_local.
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<MatrixScalar> Y_local_j =
        Y_local->getVectorNonConst (j)->get1dViewNonConst ();
      ArrayRCP<const MatrixScalar> Y_local_perm_j =
        Y_local_perm->getVector (j)->get1dView ();
      for (size_t i = 0; i < numRows_; ++i) {
        const size_t i_perm = this->ID (i);
        Y_local_j[i_perm] = Y_local_perm_j[i];
      }
    }
  }
}


//==============================================================================
// Computes Y = alpha * diag(D) * M^{-1} (diag(D) X) + beta*Y
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
weightedApply (const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
               Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
               const Tpetra::Vector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& D,
               Teuchos::ETransp mode,
               MatrixScalar alpha,
               MatrixScalar beta)
{
  using Teuchos::ArrayRCP;
  using Teuchos::Range1D;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using std::cerr;
  using std::endl;
  typedef Teuchos::ScalarTraits<MatrixScalar> STS;
  typedef Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> MV;
  typedef Tpetra::Map<MatrixLocalOrdinal, MatrixGlobalOrdinal, MatrixNode> map_type;
  const size_t numVecs = X.getNumVectors ();

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, "Ifpack2::SparseContainer::"
    "weightedApply: You must have called the compute() method before you may "
    "call apply().  You may call the apply() method as many times as you want "
    "after calling compute() once, but you must have called compute() at least "
    "once.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    Inverse_.is_null (), std::logic_error, "Ifpack2::SparseContainer::"
    "weightedApply: Inverse_ is null.  Please report this bug to the Ifpack2 "
    "developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    numVecs != Y.getNumVectors (), std::runtime_error,
    "Ifpack2::SparseContainer::weightedApply: X and Y have different numbers "
    "of vectors.  X has " << X.getNumVectors () << ", but Y has "
    << X.getNumVectors () << ".");
  if (numVecs == 0) {
    return; // done! nothing to do
  }

  // The operator Inverse_ works on a permuted subset of the local
  // parts of X and Y.  The subset and permutation are defined by
  // the GID_ array, which the ID(j) member function accesses.  If
  // the permutation is trivial and the subset is exactly equal to
  // the local indices, then we can use the local parts of X and Y
  // exactly, without needing to permute.  Otherwise, we have to use
  // temporary storage to permute X and Y.

  // Get Maps that describe the local parts of X and Y.  "Local" means
  // "the part owned by the calling MPI process only."  These Maps are
  // trivially "globally distributed" over MPI_COMM_SELF.  While this
  // does technically mean the same thing as "locally replicated" for
  // MPI_COMM_SELF, we say "globally distributed" because this more
  // generally includes communicators other than MPI_COMM_SELF, in
  // case we wish to expand that definition in the future.
  RCP<const map_type> map_X_local =
    rcp (new map_type (X.getLocalLength (), X.getMap ()->getIndexBase (),
                       LocalComm_, Tpetra::GloballyDistributed,
                       X.getMap ()->getNode ()));
  RCP<const map_type> map_Y_local =
    rcp (new map_type (Y.getLocalLength (), Y.getMap ()->getIndexBase (),
                       LocalComm_, Tpetra::GloballyDistributed,
                       Y.getMap ()->getNode ()));

  // Get views of the local parts of X and Y.
  RCP<const MV> X_local = X.offsetView (map_X_local, 0);
  RCP<MV> Y_local = Y.offsetViewNonConst (map_Y_local, 0);

  RCP<MV> X_local_perm;
  RCP<MV> Y_local_perm;
  if (needPermutation_) {
    // Create temporary permuted versions of the input and output.
    // (Re)allocate X_ and/or Y_ only if necessary.  We'll use them to
    // store the permuted versions of X resp. Y.  The current Maps of
    // X_ and Y_ don't matter, since we're treating them as local MVs.
    // Note that X_local_perm has the domain Map of the operator,
    // which may be a permuted subset of X_local's Map.  Similarly,
    // Y_local_perm has the range Map of the operator, which may be a
    // permuted subset of Y_local's Map.
    //
    // numRows_ here gives the number of rows in the row Map of the
    // local Inverse_ operator.
    //
    // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
    // here that the row Map and the range Map of that operator are
    // the same.
    //
    // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
    // really belongs in Tpetra.

    if (X_.is_null () || X_->getLocalLength () != numRows_ ||
        X_->getNumVectors () < numVecs) {
      X_ = rcp (new MV (Inverse_->getDomainMap (), numVecs));
    }
    X_local_perm = X_->subViewNonConst (Range1D (0, numVecs - 1));
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const MatrixScalar> X_local_j =
        X_local->getVector (j)->get1dView ();
      ArrayRCP<MatrixScalar> X_local_perm_j =
        X_local_perm->getVectorNonConst (j)->get1dViewNonConst ();
      for (size_t i = 0; i < numRows_; ++i) {
        const size_t i_perm = this->ID (i);
        X_local_perm_j[i] = X_local_j[i_perm];
      }
    }

    // We must permute the output multivector even on input to
    // Inverse_->apply(), since Inverse_ might use it as an initial
    // guess for a linear solve.  We have no way of knowing whether it
    // does or does not.
    if (Y_.is_null () || Y_->getLocalLength () != numRows_ ||
        Y_->getNumVectors () < numVecs) {
      Y_ = rcp (new MV (Inverse_->getRangeMap (), numVecs));
    }
    Y_local_perm = Y_->subViewNonConst (Range1D (0, numVecs - 1));
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<const MatrixScalar> Y_local_j =
        Y_local->getVector (j)->get1dView ();
      ArrayRCP<MatrixScalar> Y_local_perm_j =
        Y_local_perm->getVectorNonConst (j)->get1dViewNonConst ();
      for (size_t i = 0; i < numRows_; ++i) {
        const size_t i_perm = this->ID (i);
        Y_local_perm_j[i] = Y_local_j[i_perm];
      }
    }
  } else {
    // Don't worry about the const cast, since we won't actually write
    // to X_local_perm in this case.  The only reason for it is to
    // avoid code duplication.
    X_local_perm = rcp_const_cast<MV> (X_local);
    Y_local_perm = Y_local;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    X_local_perm.is_null (), std::logic_error, "Ifpack2::SparseContainer::"
    "weightedApply: X_local_perm is null.  Please report this bug to the "
    "Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local_perm.is_null (), std::logic_error, "Ifpack2::SparseContainer::"
    "weightedApply: Y_local_perm is null.  Please report this bug to the "
    "Ifpack2 developers.");

  // Create a temporary vector (if necessary) X_scaled for diag(D) * X.
  RCP<MV> X_scaled;
  if (needPermutation_) {
    X_scaled = rcp (new MV (Inverse_->getDomainMap (), numVecs));
  } else {
    // We're not using X_ for the permuted X, so we can use it as a
    // temporary here.  If necessary, reallocate X_ so that it has at
    // least as many columns as X.
    if (X_.is_null () || X_->getNumVectors () < numVecs) {
      X_ = rcp (new MV (Inverse_->getDomainMap (), numVecs));
    }
    X_scaled = X_->subViewNonConst (Range1D (0, numVecs - 1));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    X_scaled.is_null (), std::logic_error, "Ifpack2::SparseContainer::"
    "weightedApply: X_scaled is null.  Please report this bug to the Ifpack2 "
    "developers.");

  // X_scaled := diag(D) * X
  X_scaled->elementWiseMultiply (STS::one (), D, *X_local_perm, STS::zero ());

  // Y_temp will hold the result of M^{-1}*X_scaled.  If beta == 0, we
  // can write the result of Inverse_->apply() directly to
  // Y_local_perm, so Y_temp may alias Y_local_perm.  Otherwise, if
  // beta != 0, we need temporary storage for M^{-1}*X_scaled, so
  // Y_temp must be different than Y_local_perm.
  RCP<MV> Y_temp;
  if (beta == STS::zero ()) {
    Y_temp = Y_local_perm;
  } else if (needPermutation_) {
    Y_temp = rcp (new MV (Inverse_->getRangeMap (), numVecs));
  } else {
    // We're not using Y_ for the permuted Y, so we can use it as a
    // temporary here.  If necessary, reallocate Y_ so that it has at
    // least as many columns as Y.
    if (Y_.is_null () || Y_->getNumVectors () < numVecs) {
      Y_ = rcp (new MV (Inverse_->getRangeMap (), numVecs));
    }
    Y_temp = Y_->subViewNonConst (Range1D (0, numVecs - 1));
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_temp.is_null (), std::logic_error, "Ifpack2::SparseContainer::"
    "weightedApply: Y_temp is null.  Please report this bug to the Ifpack2 "
    "developers.");

  // Y_temp := M^{-1} * X_scaled
  Inverse_->apply (*X_scaled, *Y_temp, mode);
  // Y_local_perm := beta * Y_local_perm + alpha * diag(D) * Y_temp
  Y_local_perm->elementWiseMultiply (alpha, D, *Y_temp, beta);

  if (needPermutation_) {
    // Copy the permuted subset output vector Y_local_perm into the
    // original output multivector (view) Y_local.
    for (size_t j = 0; j < numVecs; ++j) {
      ArrayRCP<MatrixScalar> Y_local_j =
        Y_local->getVectorNonConst (j)->get1dViewNonConst ();
      ArrayRCP<const MatrixScalar> Y_local_perm_j =
        Y_local_perm->getVector (j)->get1dView ();
      for (size_t i = 0; i < numRows_; ++i) {
        const size_t i_perm = this->ID (i);
        Y_local_j[i_perm] = Y_local_perm_j[i];
      }
    }
  }
}

//==============================================================================
// Destroys all data.
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::destroy()
{
  IsInitialized_ = false;
  IsComputed_ = false;
}

//==============================================================================
template<class MatrixType, class InverseType>
std::ostream& SparseContainer<MatrixType,InverseType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}

//==============================================================================
template<class MatrixType, class InverseType>
std::string SparseContainer<MatrixType,InverseType>::description() const
{
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

  oss << "}";
  return oss.str();
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << endl;
  os << "Ifpack2_SparseContainer" << endl;
  os << "Number of rows          = " << numRows_ << endl;
  os << "Number of vectors       = " << NumVectors_ << endl;
  os << "isInitialized()         = " << IsInitialized_ << endl;
  os << "isComputed()            = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

//==============================================================================
// Extract the submatrices identified by the ID set int ID().
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::extract(const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix_in)
{
  size_t MatrixInNumRows= Matrix_in->getNodeNumRows();

  // Sanity checking
  for (size_t j = 0 ; j < numRows_ ; ++j) {
    TEUCHOS_TEST_FOR_EXCEPTION( GID_[j] < 0 || (size_t) GID_[j] >= MatrixInNumRows, std::runtime_error, "Ifpack2::SparseContainer::applyInverse compute has not been called.");
  }

  int Length = Matrix_in->getNodeMaxNumRowEntries();
  Teuchos::Array<MatrixScalar>         Values;
  Teuchos::Array<MatrixLocalOrdinal>   Indices;
  Teuchos::Array<InverseScalar>        Values_insert;
  Teuchos::Array<InverseGlobalOrdinal> Indices_insert;

  Values.resize(Length);
  Indices.resize(Length);
  Values_insert.resize(Length);
  Indices_insert.resize(Length);

  for (size_t j = 0 ; j < numRows_ ; ++j) {
    MatrixLocalOrdinal LRID = ID(j);
    size_t NumEntries;

    Matrix_in->getLocalRowCopy(LRID,Indices(),Values(),NumEntries);

    size_t num_entries_found=0;
    for (size_t k = 0 ; k < NumEntries ; ++k) {
      MatrixLocalOrdinal LCID = Indices[k];

      // skip off-processor elements
      if ((size_t)LCID >= MatrixInNumRows)
        continue;

      // for local column IDs, look for each ID in the list
      // of columns hosted by this object
      // FIXME: use STL
      InverseLocalOrdinal jj = -1;
      for (size_t kk = 0 ; kk < numRows_ ; ++kk)
        if (ID(kk) == LCID)
          jj = kk;

      if (jj != -1) {
        Indices_insert[num_entries_found] = Map_->getGlobalElement(jj);
        Values_insert[num_entries_found]  = Values[k];
        num_entries_found++;
      }

    }
    Matrix_->insertGlobalValues(j,Indices_insert(0,num_entries_found),Values_insert(0,num_entries_found));
  }

  Matrix_->fillComplete();
}


} // namespace Ifpack2
#endif // IFPACK2_SPARSECONTAINER_HPP
