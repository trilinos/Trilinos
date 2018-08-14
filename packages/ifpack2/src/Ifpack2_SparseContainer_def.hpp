/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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
#include "Teuchos_TestForException.hpp"

namespace Ifpack2 {

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType,InverseType>::
SparseContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<Teuchos::Array<local_ordinal_type> >& partitions,
                 const Teuchos::RCP<const import_type>& importer,
                 int OverlapLevel,
                 scalar_type DampingFactor) :
  Container<MatrixType> (matrix, partitions, importer, OverlapLevel,
                         DampingFactor),
  IsInitialized_ (false),
  IsComputed_ (false),
#ifdef HAVE_MPI
  localComm_ (Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF)))
#else
  localComm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ()))
#endif // HAVE_MPI
{
  global_ordinal_type indexBase = 0;
  for(int i = 0; i < this->numBlocks_; i++)
    localMaps_.emplace_back(this->blockRows_[i], indexBase, localComm_);
}

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType, InverseType>::
SparseContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<local_ordinal_type>& localRows) :
  Container<MatrixType> (matrix, localRows),
  IsInitialized_(false),
  IsComputed_(false),
#ifdef HAVE_MPI
  localComm_ (Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_SELF)))
#else
  localComm_ (Teuchos::rcp(new Teuchos::SerialComm<int>()))
#endif // HAVE_MPI
{
  global_ordinal_type indexBase = 0;
  localMaps_.emplace_back(this->blockRows_[0], indexBase, localComm_);
}

//==============================================================================
template<class MatrixType,class InverseType>
SparseContainer<MatrixType,InverseType>::~SparseContainer()
{
  for(auto inv : Inverses_)
    delete inv.get();
}

//==============================================================================
template<class MatrixType, class InverseType>
bool SparseContainer<MatrixType,InverseType>::isInitialized() const
{
  return IsInitialized_;
}

//==============================================================================
template<class MatrixType, class InverseType>
bool SparseContainer<MatrixType,InverseType>::isComputed() const
{
  return IsComputed_;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::setParameters(const Teuchos::ParameterList& List)
{
  List_ = List;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::initialize ()
{
  using Teuchos::RCP;
  // We assume that if you called this method, you intend to recompute
  // everything.  Thus, we release references to all the internal
  // objects.  We do this first to save memory.  (In an RCP
  // assignment, the right-hand side and left-hand side coexist before
  // the left-hand side's reference count gets updated.)
  IsInitialized_ = false;
  IsComputed_ = false;

  // (Re)create the CrsMatrices that will contain the
  // local matrices to use for solves.
  diagBlocks_.reserve(this->numBlocks_);
  Inverses_.reserve(this->numBlocks_);
  for(int i = 0; i < this->numBlocks_; i++)
  {
    // Create a local map for the block, with same size as block has rows.
    // Note: this map isn't needed elsewhere in SparseContainer, but the
    // diagBlocks_[...] will keep it alive
    RCP<InverseMap> tempMap(new InverseMap(this->blockRows_[i], 0, localComm_));
    diagBlocks_.emplace_back(new InverseCrs(tempMap, 0));
    // Create the inverse operator using the local matrix.  We give it
    // the matrix here, but don't call its initialize() or compute()
    // methods yet, since we haven't initialized the matrix yet.
    Inverses_.push_back(ptr(new InverseType(diagBlocks_[i])));
  }
  IsInitialized_ = true;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::compute ()
{
  IsComputed_ = false;
  if (! this->isInitialized ()) {
    this->initialize ();
  }

  // Extract the submatrix.
  this->extract ();

  // The inverse operator already has a pointer to the submatrix.  Now
  // that the submatrix is filled in, we can initialize and compute
  // the inverse operator.
  for(int i = 0; i < this->numBlocks_; i++)
    Inverses_[i]->initialize ();
  for(int i = 0; i < this->numBlocks_; i++)
    Inverses_[i]->compute ();
  IsComputed_ = true;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType, InverseType>::clearBlocks ()
{
  for(auto inv : Inverses_)
    delete inv.get();
  Inverses_.clear();
  diagBlocks_.clear();
  Container<MatrixType>::clearBlocks();
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
applyImpl (inverse_mv_type& X,
           inverse_mv_type& Y,
           int blockIndex,
           int /* stride */,
           Teuchos::ETransp mode,
           InverseScalar alpha,
           InverseScalar beta) const
{
    TEUCHOS_TEST_FOR_EXCEPTION(
      Inverses_[blockIndex]->getDomainMap()->getNodeNumElements() != X.getLocalLength(),
      std::logic_error, "Ifpack2::SparseContainer::apply: Inverse_ "
      "operator and X have incompatible dimensions (" <<
      Inverses_[blockIndex]->getDomainMap()->getNodeNumElements() << " resp. "
      << X.getLocalLength() << ").  Please report this bug to "
      "the Ifpack2 developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      Inverses_[blockIndex]->getRangeMap()->getNodeNumElements() != Y.getLocalLength(),
      std::logic_error, "Ifpack2::SparseContainer::apply: Inverse_ "
      "operator and Y have incompatible dimensions (" <<
      Inverses_[blockIndex]->getRangeMap()->getNodeNumElements() << " resp. "
      << Y.getLocalLength() << ").  Please report this bug to "
      "the Ifpack2 developers.");
    Inverses_[blockIndex]->apply(X, Y, mode, alpha, beta);
}

template<class MatrixType, class InverseType>
void SparseContainer<MatrixType, InverseType>::
apply (HostView& X,
       HostView& Y,
       int blockIndex,
       int stride,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::ArrayView;
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

  // Tpetra::MultiVector specialization corresponding to InverseType.
  Details::MultiVectorLocalGatherScatter<mv_type, inverse_mv_type> mvgs;
  size_t numVecs = X.extent(1);

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, "Ifpack2::SparseContainer::apply: "
    "You must have called the compute() method before you may call apply().  "
    "You may call the apply() method as many times as you want after calling "
    "compute() once, but you must have called compute() at least once.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent(1) != Y.extent(1), std::runtime_error,
    "Ifpack2::SparseContainer::apply: X and Y have different numbers of "
    "vectors.  X has " << X.extent(1)
    << ", but Y has " << Y.extent(1) << ".");

  if (numVecs == 0) {
    return; // done! nothing to do
  }

  const local_ordinal_type numRows_ = this->blockRows_[blockIndex];

  // The operator Inverse_ works on a permuted subset of the local
  // parts of X and Y.  The subset and permutation are defined by the
  // index array returned by getLocalRows().  If the permutation is
  // trivial and the subset is exactly equal to the local indices,
  // then we could use the local parts of X and Y exactly, without
  // needing to permute.  Otherwise, we have to use temporary storage
  // to permute X and Y.  For now, we always use temporary storage.
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

  if(invX.size() == 0)
  {
    for(local_ordinal_type i = 0; i < this->numBlocks_; i++)
      invX.emplace_back(Inverses_[i]->getDomainMap(), numVecs);
    for(local_ordinal_type i = 0; i < this->numBlocks_; i++)
      invY.emplace_back(Inverses_[i]->getDomainMap(), numVecs);
  }
  inverse_mv_type& X_local = invX[blockIndex];
  TEUCHOS_TEST_FOR_EXCEPTION(
    X_local.getLocalLength() != (size_t) numRows_, std::logic_error,
    "Ifpack2::SparseContainer::apply: "
    "X_local has length " << X_local.getLocalLength() << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  const ArrayView<const local_ordinal_type> localRows = this->getLocalRows(blockIndex);
  mvgs.gatherMVtoView(X_local, X, localRows);

  // We must gather the output multivector Y even on input to
  // Inverse_->apply(), since the Inverse_ operator might use it as an
  // initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  inverse_mv_type& Y_local = invY[blockIndex];
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local.getLocalLength () != (size_t) numRows_, std::logic_error,
    "Ifpack2::SparseContainer::apply: "
    "Y_local has length " << Y_local.getLocalLength () << ", which does "
    "not match numRows_ = " << numRows_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gatherMVtoView(Y_local, Y, localRows);

  // Apply the local operator:
  // Y_local := beta*Y_local + alpha*M^{-1}*X_local
  this->applyImpl(X_local, Y_local, blockIndex, stride, mode, as<InverseScalar>(alpha),
                  as<InverseScalar>(beta));

  // Scatter the permuted subset output vector Y_local back into the
  // original output multivector Y.
  mvgs.scatterMVtoView(Y, Y_local, localRows);
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType, InverseType>::
weightedApply (HostView& X,
               HostView& Y,
               HostView& D,
               int blockIndex,
               int /* stride */,
               Teuchos::ETransp mode,
               scalar_type alpha,
               scalar_type beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Range1D;
  using Teuchos::Ptr;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using std::cerr;
  using std::endl;
  typedef Teuchos::ScalarTraits<scalar_type> STS;

  // The InverseType template parameter might have different template
  // parameters (Scalar, LO, GO, and/or Node) than MatrixType.  For
  // example, MatrixType (a global object) might use a bigger GO
  // (global ordinal) type than InverseType (which deals with the
  // diagonal block, a local object).  This means that we might have
  // to convert X and Y to the Tpetra::MultiVector specialization that
  // InverseType wants.  This class' X_ and Y_ internal fields are of
  // the right type for InverseType, so we can use those as targets.

  // Tpetra::Vector specialization corresponding to InverseType.
  typedef Tpetra::Vector<InverseScalar, InverseLocalOrdinal, InverseGlobalOrdinal, InverseNode> inverse_vector_type;

  Details::MultiVectorLocalGatherScatter<mv_type, inverse_mv_type> mvgs;
  const size_t numVecs = X.extent(1);

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error, "Ifpack2::SparseContainer::"
    "weightedApply: You must have called the compute() method before you may "
    "call apply().  You may call the apply() method as many times as you want "
    "after calling compute() once, but you must have called compute() at least "
    "once.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent(1) != Y.extent(1), std::runtime_error,
    "Ifpack2::SparseContainer::weightedApply: X and Y have different numbers "
    "of vectors.  X has " << X.extent(1) << ", but Y has "
    << Y.extent(1) << ".");

  if (numVecs == 0) {
    return; // done! nothing to do
  }

  // The operator Inverse_ works on a permuted subset of the local
  // parts of X and Y.  The subset and permutation are defined by the
  // index array returned by getLocalRows().  If the permutation is
  // trivial and the subset is exactly equal to the local indices,
  // then we could use the local parts of X and Y exactly, without
  // needing to permute.  Otherwise, we have to use temporary storage
  // to permute X and Y.  For now, we always use temporary storage.
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
  const local_ordinal_type numRows = this->blockRows_[blockIndex];

  if(invX.size() == 0)
  {
    for(local_ordinal_type i = 0; i < this->numBlocks_; i++)
      invX.emplace_back(Inverses_[i]->getDomainMap(), numVecs);
    for(local_ordinal_type i = 0; i < this->numBlocks_; i++)
      invY.emplace_back(Inverses_[i]->getDomainMap(), numVecs);
  }
  inverse_mv_type& X_local = invX[blockIndex];
  const ArrayView<const local_ordinal_type> localRows = this->getLocalRows(blockIndex);
  mvgs.gatherMVtoView(X_local, X, localRows);

  // We must gather the output multivector Y even on input to
  // Inverse_->apply(), since the Inverse_ operator might use it as an
  // initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  inverse_mv_type Y_local = invY[blockIndex];
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local.getLocalLength() != (size_t) numRows, std::logic_error,
    "Ifpack2::SparseContainer::weightedApply: "
    "Y_local has length " << X_local.getLocalLength() << ", which does "
    "not match numRows_ = " << numRows << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gatherMVtoView(Y_local, Y, localRows);

  // Apply the diagonal scaling D to the input X.  It's our choice
  // whether the result has the original input Map of X, or the
  // permuted subset Map of X_local.  If the latter, we also need to
  // gather D into the permuted subset Map.  We choose the latter, to
  // save memory and computation.  Thus, we do the following:
  //
  // 1. Gather D into a temporary vector D_local.
  // 2. Create a temporary X_scaled to hold diag(D_local) * X_local.
  // 3. Compute X_scaled := diag(D_loca) * X_local.

  inverse_vector_type D_local(Inverses_[blockIndex]->getDomainMap());
  TEUCHOS_TEST_FOR_EXCEPTION(
    D_local.getLocalLength() != (size_t) this->blockRows_[blockIndex], std::logic_error,
    "Ifpack2::SparseContainer::weightedApply: "
    "D_local has length " << X_local.getLocalLength () << ", which does "
    "not match numRows_ = " << this->blockRows_[blockIndex] << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gatherMVtoView(D_local, D, localRows);
  inverse_mv_type X_scaled(Inverses_[blockIndex]->getDomainMap(), numVecs);
  X_scaled.elementWiseMultiply(STS::one(), D_local, X_local, STS::zero());

  // Y_temp will hold the result of M^{-1}*X_scaled.  If beta == 0, we
  // can write the result of Inverse_->apply() directly to Y_local, so
  // Y_temp may alias Y_local.  Otherwise, if beta != 0, we need
  // temporary storage for M^{-1}*X_scaled, so Y_temp must be
  // different than Y_local.
  Ptr<inverse_mv_type> Y_temp;
  bool deleteYT = false;
  if (beta == STS::zero ()) {
    Y_temp = ptr(&Y_local);
  } else {
    Y_temp = ptr(new inverse_mv_type(Inverses_[blockIndex]->getRangeMap(), numVecs));
    deleteYT = true;
  }
  // Apply the local operator: Y_temp := M^{-1} * X_scaled
  Inverses_[blockIndex]->apply(X_scaled, *Y_temp, mode);
  // Y_local := beta * Y_local + alpha * diag(D_local) * Y_tmp.
  //
  // Note that we still use the permuted subset scaling D_local here,
  // because Y_temp has the same permuted subset Map.  That's good, in
  // fact, because it's a subset: less data to read and multiply.
  Y_local.elementWiseMultiply(alpha, D_local, *Y_temp, beta);
  if(deleteYT)
    delete Y_temp.get();
  // Copy the permuted subset output vector Y_local into the original
  // output multivector Y.
  mvgs.scatterMVtoView(Y, Y_local, localRows);
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
  oss << "\"Ifpack2::SparseContainer\": {";
  if (isInitialized()) {
    if (isComputed()) {
      oss << "status = initialized, computed";
    }
    else {
      oss << "status = initialized, not computed";
    }
  }
  else {
    oss << "status = not initialized, not computed";
  }
  for(int i = 0; i < this->numBlocks_; i++)
  {
    oss << ", Block Inverse " << i << ": {";
    oss << Inverses_[i]->description();
    oss << "}";
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
  os << "Ifpack2::SparseContainer" << endl;
  for(int i = 0; i < this->numBlocks_; i++)
  {
    os << "Block " << i << " rows: = " << this->blockRows_[i] << endl;
  }
  os << "isInitialized()         = " << IsInitialized_ << endl;
  os << "isComputed()            = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
extract ()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;

  auto& A = *this->inputMatrix_;
  const size_t MatrixInNumRows = A.getNodeNumRows ();

  // Sanity checking
  for(local_ordinal_type i = 0; i < this->numBlocks_; i++)
  {
    const ArrayView<const local_ordinal_type> localRows = this->getLocalRows(i);
    for (local_ordinal_type j = 0; j < localRows.size(); j++)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        localRows[j] < 0 ||
        static_cast<size_t> (localRows[j]) >= MatrixInNumRows,
        std::runtime_error, "Ifpack2::SparseContainer::extract: localRows[j="
        << j << "] = " << localRows[j] << ", which is out of the valid range.  "
        "This probably means that compute() has not yet been called.");
    }
  }

  const size_t maxNumEntriesInRow = A.getNodeMaxNumRowEntries();
  Array<scalar_type> Values;
  Array<local_ordinal_type> Indices;
  Array<InverseScalar> Values_insert;
  Array<InverseGlobalOrdinal> Indices_insert;

  Values.resize (maxNumEntriesInRow);
  Indices.resize (maxNumEntriesInRow);
  Values_insert.resize (maxNumEntriesInRow);
  Indices_insert.resize (maxNumEntriesInRow);

  const InverseLocalOrdinal INVALID = Teuchos::OrdinalTraits<InverseLocalOrdinal>::invalid ();
  for(local_ordinal_type i = 0; i < this->numBlocks_; i++)
  {
    const local_ordinal_type numRows_ = this->blockRows_[i];
    const ArrayView<const local_ordinal_type> localRows = this->getLocalRows(i);
    for (local_ordinal_type j = 0; j < numRows_; j++)
    {
      const local_ordinal_type localRow = this->partitions_[this->partitionIndices_[i] + j];
      size_t numEntries;
      A.getLocalRowCopy(localRow, Indices(), Values(), numEntries);
      size_t num_entries_found = 0;
      for(size_t k = 0; k < numEntries; ++k)
      {
        const local_ordinal_type localCol = Indices[k];
        // Skip off-process elements
        //
        // FIXME (mfh 24 Aug 2013) This assumes the following:
        //
        // 1. The column and row Maps begin with the same set of
        //    on-process entries, in the same order.  That is,
        //    on-process row and column indices are the same.
        // 2. All off-process indices in the column Map of the input
        //    matrix occur after that initial set.
        if(static_cast<size_t> (localCol) >= MatrixInNumRows)
          continue;
        // for local column IDs, look for each ID in the list
        // of columns hosted by this object
        InverseLocalOrdinal jj = INVALID;
        for(local_ordinal_type kk = 0; kk < numRows_; kk++)
        {
          if(localRows[kk] == localCol)
            jj = kk;
        }
        if (jj != INVALID)
        {
          Indices_insert[num_entries_found] = localMaps_[i].getGlobalElement(jj);
          Values_insert[num_entries_found] = Values[k];
          num_entries_found++;
        }
      }
      diagBlocks_[i]->insertGlobalValues(j, Indices_insert (0, num_entries_found),
                                        Values_insert (0, num_entries_found));
    }
    // FIXME (mfh 24 Aug 2013) If we generalize the current set of
    // assumptions on the column and row Maps (see note above), we may
    // need to supply custom domain and range Maps to fillComplete().
    diagBlocks_[i]->fillComplete ();
  }
}

template<typename MatrixType, typename InverseType>
std::string SparseContainer<MatrixType, InverseType>::getName()
{
  typedef ILUT<Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> > ILUTInverse;
#ifdef HAVE_IFPACK2_AMESOS2
  typedef Details::Amesos2Wrapper<Tpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type>> AmesosInverse;
  if(std::is_same<InverseType, ILUTInverse>::value)
  {
    return "SparseILUT";
  }
  else if(std::is_same<InverseType, AmesosInverse>::value)
  {
    return "SparseAmesos";
  }
  else
  {
    throw std::logic_error("InverseType for SparseContainer must be Ifpack2::ILUT or Details::Amesos2Wrapper");
  }
#else
  // Macros can't have commas in their arguments, so we have to
  // compute the bool first argument separately.
  constexpr bool inverseTypeIsILUT = std::is_same<InverseType, ILUTInverse>::value;
  TEUCHOS_TEST_FOR_EXCEPTION(! inverseTypeIsILUT, std::logic_error,
    "InverseType for SparseContainer must be Ifpack2::ILUT<ROW>");
  return "SparseILUT";    //the only supported sparse container specialization if no Amesos2
#endif
}

} // namespace Ifpack2

// For ETI
#include "Ifpack2_ILUT.hpp"
#ifdef HAVE_IFPACK2_AMESOS2
#include "Ifpack2_Details_Amesos2Wrapper.hpp"
#endif

// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.

#ifdef HAVE_IFPACK2_AMESOS2
#  define IFPACK2_SPARSECONTAINER_INSTANT(S,LO,GO,N) \
    template class Ifpack2::SparseContainer< Tpetra::RowMatrix<S, LO, GO, N>, \
                                             Ifpack2::ILUT<Tpetra::RowMatrix<S,LO,GO,N> > >; \
    template class Ifpack2::SparseContainer< Tpetra::RowMatrix<S, LO, GO, N>, \
                                             Ifpack2::Details::Amesos2Wrapper<Tpetra::RowMatrix<S,LO,GO,N> > >;
#else
#  define IFPACK2_SPARSECONTAINER_INSTANT(S,LO,GO,N) \
    template class Ifpack2::SparseContainer< Tpetra::RowMatrix<S,LO,GO,N>, \
                                             Ifpack2::ILUT<Tpetra::RowMatrix<S, LO, GO, N> > >;
#endif
#endif // IFPACK2_SPARSECONTAINER_DEF_HPP
