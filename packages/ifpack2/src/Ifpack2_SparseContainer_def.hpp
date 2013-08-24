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
  IsInitialized_ (false),
  IsComputed_ (false),
#ifdef HAVE_MPI
  LocalComm_ (Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF)))
#else
  LocalComm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ()))
#endif // HAVE_MPI
{
  (void) NumVectors;
}

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
  using Teuchos::rcp;
  typedef Tpetra::Map<InverseLocalOrdinal, InverseGlobalOrdinal,
                      InverseNode> map_type;
  typedef Tpetra::CrsMatrix<InverseScalar, InverseLocalOrdinal,
                            InverseGlobalOrdinal, InverseNode> crs_matrix_type;
  if (IsInitialized_) {
    destroy ();
  }
  IsInitialized_ = false;

  Map_ = rcp (new map_type (numRows_, 0, LocalComm_));
  Matrix_ = rcp (new crs_matrix_type (Map_, 0));
  GID_.resize (numRows_);

  // Create the inverse operator on the local block.  We give it the
  // matrix here, but don't call its initialize() or compute() methods
  // yet, since we haven't initialized the matrix yet.
  Inverse_ = rcp (new InverseType (Matrix_));
  Inverse_->setParameters (List_);

  // Note from IFPACK: Call Inverse_->initialize() in compute(). This
  // saves some time, because the diagonal blocks can be extracted
  // faster, and only once.
  IsInitialized_ = true;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
compute (const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix)
{
  IsComputed_ = false;
  if (! this->isInitialized ()) {
    this->initialize ();
  }

  // Extract the submatrix.
  this->extract (Matrix);

  // initialize & compute the inverse operator
  Inverse_->initialize ();
  Inverse_->compute ();

  IsComputed_ = true;
}


//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
applyImpl (const Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>& X,
           Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode>& Y,
           Teuchos::ETransp mode,
           InverseScalar alpha,
           InverseScalar beta) const
{
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
  TEUCHOS_TEST_FOR_EXCEPTION(
    Inverse_->getDomainMap ()->getNodeNumElements () != X.getLocalLength (),
    std::logic_error, "Ifpack2::SparseContainer::apply: Inverse_ "
    "operator and X have incompatible dimensions (" <<
    Inverse_->getDomainMap ()->getNodeNumElements () << " resp. "
    << X.getLocalLength () << ").  Please report this bug to "
    "the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    Inverse_->getRangeMap ()->getNodeNumElements () != Y.getLocalLength (),
    std::logic_error, "Ifpack2::SparseContainer::apply: Inverse_ "
    "operator and Y have incompatible dimensions (" <<
    Inverse_->getRangeMap ()->getNodeNumElements () << " resp. "
    << Y.getLocalLength () << ").  Please report this bug to "
    "the Ifpack2 developers.");

  Inverse_->apply (X, Y, mode, alpha, beta);
}


//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
apply (const Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& X,
       Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode>& Y,
       Teuchos::ETransp mode,
       MatrixScalar alpha,
       MatrixScalar beta) const
{
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // The InverseType template parameter might have different template
  // parameters (Scalar, LO, GO, and/or Node) than MatrixType.  For
  // example, MatrixType (a global object) might use a bigger GO
  // (global ordinal) type than InverseType (which deals with the
  // diagonal block, a local object).  This means that we might have
  // to convert X and Y to the Tpetra::MultiVector specialization that
  // InverseType wants.  This class' X_ and Y_ internal fields are of
  // the right type for InverseType, so we can use those as targets.

  // Tpetra::MultiVector specialization corresponding to MatrixType.
  typedef Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> MV_mat;
  // Tpetra::MultiVector specialization corresponding to InverseType.
  typedef Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> MV_inv;
  MultiVectorLocalGatherScatter<MV_mat, MV_inv> mvgs;
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
  // local indices, then we could use the local parts of X and Y
  // exactly, without needing to permute.  Otherwise, we have to use
  // temporary storage to permute X and Y.  For now, we always use
  // temporary storage.
  //
  // Create temporary permuted versions of the input and output.
  // (Re)allocate X_ and/or Y_ only if necessary.  We'll use them to
  // store the permuted versions of X resp. Y.  Note that X_local has
  // the domain Map of the operator, which may be a permuted subset of
  // the local Map corresponding to X.getMap().  Similarly, Y_local
  // has the range Map of the operator, which may be a permuted subset
  // of the local Map corresponding to Y.getMap().  numRows_ here
  // gives the number of rows in the row Map of the local Inverse_
  // operator.
  //
  // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
  // here that the row Map and the range Map of that operator are
  // the same.
  //
  // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
  // really belongs in Tpetra.

  if (X_.is_null ()) {
    X_ = rcp (new MV_inv (Inverse_->getDomainMap (), numVecs));
  }
  RCP<MV_inv> X_local = X_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    X_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::SparseContainer::apply: "
    "X_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gather (*X_local, X, this->GID_);

  // We must gather the output multivector Y even on input to
  // Inverse_->apply(), since the Inverse_ operator might use it as an
  // initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  if (Y_.is_null ()) {
    Y_ = rcp (new MV_inv (Inverse_->getRangeMap (), numVecs));
  }
  RCP<MV_inv> Y_local = Y_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::SparseContainer::apply: "
    "Y_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gather (*Y_local, Y, this->GID_);

  // Apply the local operator:
  // Y_local := beta*Y_local + alpha*M^{-1}*X_local
  this->applyImpl (*X_local, *Y_local, mode, as<InverseScalar> (alpha),
                   as<InverseScalar> (beta));

  // Scatter the permuted subset output vector Y_local back into the
  // original output multivector Y.
  mvgs.scatter (Y, *Y_local, this->GID_);
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
               MatrixScalar beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::Range1D;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using std::cerr;
  using std::endl;
  typedef Teuchos::ScalarTraits<MatrixScalar> STS;

  // The InverseType template parameter might have different template
  // parameters (Scalar, LO, GO, and/or Node) than MatrixType.  For
  // example, MatrixType (a global object) might use a bigger GO
  // (global ordinal) type than InverseType (which deals with the
  // diagonal block, a local object).  This means that we might have
  // to convert X and Y to the Tpetra::MultiVector specialization that
  // InverseType wants.  This class' X_ and Y_ internal fields are of
  // the right type for InverseType, so we can use those as targets.

  // Tpetra::MultiVector specialization corresponding to MatrixType.
  typedef Tpetra::MultiVector<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> MV_mat;
  // Tpetra::MultiVector specialization corresponding to InverseType.
  typedef Tpetra::MultiVector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> MV_inv;
  // Tpetra::Vector specialization corresponding to InverseType.
  typedef Tpetra::Vector<InverseScalar,InverseLocalOrdinal,InverseGlobalOrdinal,InverseNode> V_inv;
  MultiVectorLocalGatherScatter<MV_mat, MV_inv> mvgs;
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
  // parts of X and Y.  The subset and permutation are defined by the
  // GID_ array, which the ID(j) member function accesses.  If the
  // permutation is trivial and the subset is exactly equal to the
  // local indices, then we could use the local parts of X and Y
  // exactly, without needing to permute.  Otherwise, we have to use
  // temporary storage to permute X and Y.  For now, we always use
  // temporary storage.
  //
  // Create temporary permuted versions of the input and output.
  // (Re)allocate X_ and/or Y_ only if necessary.  We'll use them to
  // store the permuted versions of X resp. Y.  Note that X_local has
  // the domain Map of the operator, which may be a permuted subset of
  // the local Map corresponding to X.getMap().  Similarly, Y_local
  // has the range Map of the operator, which may be a permuted subset
  // of the local Map corresponding to Y.getMap().  numRows_ here
  // gives the number of rows in the row Map of the local Inverse_
  // operator.
  //
  // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
  // here that the row Map and the range Map of that operator are
  // the same.
  //
  // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
  // really belongs in Tpetra.

  if (X_.is_null ()) {
    X_ = rcp (new MV_inv (Inverse_->getDomainMap (), numVecs));
  }
  RCP<MV_inv> X_local = X_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    X_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::SparseContainer::weightedApply: "
    "X_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gather (*X_local, X, this->GID_);

  // We must gather the output multivector Y even on input to
  // Inverse_->apply(), since the Inverse_ operator might use it as an
  // initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  if (Y_.is_null ()) {
    Y_ = rcp (new MV_inv (Inverse_->getRangeMap (), numVecs));
  }
  RCP<MV_inv> Y_local = Y_;
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::SparseContainer::weightedApply: "
    "Y_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gather (*Y_local, Y, this->GID_);

  // Apply the diagonal scaling D to the input X.  It's our choice
  // whether the result has the original input Map of X, or the
  // permuted subset Map of X_local.  If the latter, we also need to
  // gather D into the permuted subset Map.  We choose the latter, to
  // save memory and computation.  Thus, we do the following:
  //
  // 1. Gather D into a temporary vector D_local.
  // 2. Create a temporary X_scaled to hold diag(D_local) * X_local.
  // 3. Compute X_scaled := diag(D_loca) * X_local.

  RCP<V_inv> D_local = rcp (new V_inv (Inverse_->getDomainMap ()));
  TEUCHOS_TEST_FOR_EXCEPTION(
    D_local->getLocalLength () != numRows_, std::logic_error,
    "Ifpack2::SparseContainer::weightedApply: "
    "D_local has length " << X_local->getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gather (*D_local, D, this->GID_);
  RCP<MV_inv> X_scaled = rcp (new MV_inv (Inverse_->getDomainMap (), numVecs));
  X_scaled->elementWiseMultiply (STS::one (), *D_local, *X_local, STS::zero ());

  // Y_temp will hold the result of M^{-1}*X_scaled.  If beta == 0, we
  // can write the result of Inverse_->apply() directly to Y_local, so
  // Y_temp may alias Y_local.  Otherwise, if beta != 0, we need
  // temporary storage for M^{-1}*X_scaled, so Y_temp must be
  // different than Y_local.
  RCP<MV_inv> Y_temp;
  if (beta == STS::zero ()) {
    Y_temp = Y_local;
  } else {
    Y_temp = rcp (new MV_inv (Inverse_->getRangeMap (), numVecs));
  }

  // Apply the local operator: Y_temp := M^{-1} * X_scaled
  Inverse_->apply (*X_scaled, *Y_temp, mode);
  // Y_local := beta * Y_local + alpha * diag(D_local) * Y_temp.
  //
  // Note that we still use the permuted subset scaling D_local here,
  // because Y_temp has the same permuted subset Map.  That's good, in
  // fact, because it's a subset: less data to read and multiply.
  Y_local->elementWiseMultiply (alpha, *D_local, *Y_temp, beta);

  // Copy the permuted subset output vector Y_local into the original
  // output multivector Y.
  mvgs.scatter (Y, *Y_local, this->GID_);
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
  os << "isInitialized()         = " << IsInitialized_ << endl;
  os << "isComputed()            = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
extract (const Teuchos::RCP<const Tpetra::RowMatrix<MatrixScalar,MatrixLocalOrdinal,MatrixGlobalOrdinal,MatrixNode> >& Matrix_in)
{
  const size_t MatrixInNumRows = Matrix_in->getNodeNumRows ();

  // Sanity checking
  for (size_t j = 0; j < numRows_; ++j) {
    TEUCHOS_TEST_FOR_EXCEPTION(
      GID_[j] < 0 || (size_t) GID_[j] >= MatrixInNumRows,
      std::runtime_error, "Ifpack2::SparseContainer::extract: "
      "GID_[j=" << j << "] = " << GID_[j] << ", which is out of the valid "
      "range.  This probably means that compute() has not yet been called.");
  }

  int Length = Matrix_in->getNodeMaxNumRowEntries();
  Teuchos::Array<MatrixScalar>         Values;
  Teuchos::Array<MatrixLocalOrdinal>   Indices;
  Teuchos::Array<InverseScalar>        Values_insert;
  Teuchos::Array<InverseGlobalOrdinal> Indices_insert;

  Values.resize (Length);
  Indices.resize (Length);
  Values_insert.resize (Length);
  Indices_insert.resize (Length);

  const InverseLocalOrdinal INVALID =
    Teuchos::OrdinalTraits<InverseLocalOrdinal>::invalid ();
  for (size_t j = 0; j < numRows_; ++j) {
    const MatrixLocalOrdinal LRID = this->ID (j);
    size_t NumEntries;
    Matrix_in->getLocalRowCopy (LRID, Indices (), Values (), NumEntries);

    size_t num_entries_found = 0;
    for (size_t k = 0; k < NumEntries; ++k) {
      const MatrixLocalOrdinal LCID = Indices[k];

      // Skip off-process elements
      //
      // FIXME (mfh 24 Aug 2013) This assumes the following:
      //
      // 1. The column and row Maps begin with the same set of
      //    on-process entries.
      //
      // 2. All off-process indices in the column Map of the input
      //    matrix occur after that initial set.
      if (static_cast<size_t> (LCID) >= MatrixInNumRows) {
        continue;
      }
      // for local column IDs, look for each ID in the list
      // of columns hosted by this object
      InverseLocalOrdinal jj = INVALID;
      for (size_t kk = 0; kk < numRows_; ++kk) {
        if (ID(kk) == LCID) {
          jj = kk;
        }
      }

      if (jj != INVALID) {
        Indices_insert[num_entries_found] = Map_->getGlobalElement(jj);
        Values_insert[num_entries_found]  = Values[k];
        num_entries_found++;
      }
    }
    Matrix_->insertGlobalValues (j, Indices_insert (0, num_entries_found),
                                 Values_insert (0, num_entries_found));
  }

  Matrix_->fillComplete();
}


} // namespace Ifpack2
#endif // IFPACK2_SPARSECONTAINER_HPP
