// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
                 const Teuchos::Array<Teuchos::Array<LO> >& partitions,
                 const Teuchos::RCP<const import_type>&,
                 bool pointIndexed) :
  ContainerImpl<MatrixType, InverseScalar> (matrix, partitions, pointIndexed),
#ifdef HAVE_MPI
  localComm_ (Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF)))
#else
  localComm_ (Teuchos::rcp (new Teuchos::SerialComm<int> ()))
#endif // HAVE_MPI
{}

//==============================================================================
template<class MatrixType, class InverseType>
SparseContainer<MatrixType, InverseType>::
~SparseContainer() {}

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
  // We assume that if you called this method, you intend to recompute
  // everything.  Thus, we release references to all the internal
  // objects.  We do this first to save memory.  (In an RCP
  // assignment, the right-hand side and left-hand side coexist before
  // the left-hand side's reference count gets updated.)
  Teuchos::RCP<Teuchos::Time> timer =
    Teuchos::TimeMonitor::getNewCounter ("Ifpack2::SparseContainer::initialize");
  Teuchos::TimeMonitor timeMon (*timer);

  //Will create the diagonal blocks and their inverses
  //in extract()
  diagBlocks_.assign(this->numBlocks_, Teuchos::null);
  Inverses_.assign(this->numBlocks_, Teuchos::null);

  // Extract the submatrices.
  this->extractGraph();

  // Initialize the inverse operator.
  for(int i = 0; i < this->numBlocks_; i++)
  {
    Inverses_[i]->setParameters(List_);
    Inverses_[i]->initialize ();
  }

  this->IsInitialized_ = true;
  this->IsComputed_ = false;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::compute ()
{
  Teuchos::RCP<Teuchos::Time> timer =
    Teuchos::TimeMonitor::getNewCounter ("Ifpack2::SparseContainer::compute");
  Teuchos::TimeMonitor timeMon (*timer);

  this->IsComputed_ = false;
  if (!this->isInitialized ()) {
    this->initialize ();
  }

  // Extract the submatrices values
  this->extractValues();

  // Compute the inverse operator.
  for(int i = 0; i < this->numBlocks_; i++) {
    Inverses_[i]->compute ();
  }

  this->IsComputed_ = true;
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType, InverseType>::clearBlocks ()
{
  for(auto inv : Inverses_)
    delete inv.get();
  Inverses_.clear();
  diagBlocks_.clear();
  ContainerImpl<MatrixType, InverseScalar>::clearBlocks();
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
solveBlockMV(const inverse_mv_type& X,
             inverse_mv_type& Y,
             int blockIndex,
             Teuchos::ETransp mode,
             InverseScalar alpha,
             InverseScalar beta) const
{
    TEUCHOS_TEST_FOR_EXCEPTION(
      Inverses_[blockIndex]->getDomainMap()->getLocalNumElements() != X.getLocalLength(),
      std::logic_error, "Ifpack2::SparseContainer::apply: Inverse_ "
      "operator and X have incompatible dimensions (" <<
      Inverses_[blockIndex]->getDomainMap()->getLocalNumElements() << " resp. "
      << X.getLocalLength() << ").  Please report this bug to "
      "the Ifpack2 developers.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      Inverses_[blockIndex]->getRangeMap()->getLocalNumElements() != Y.getLocalLength(),
      std::logic_error, "Ifpack2::SparseContainer::apply: Inverse_ "
      "operator and Y have incompatible dimensions (" <<
      Inverses_[blockIndex]->getRangeMap()->getLocalNumElements() << " resp. "
      << Y.getLocalLength() << ").  Please report this bug to "
      "the Ifpack2 developers.");
    Inverses_[blockIndex]->apply(X, Y, mode, alpha, beta);
}

template<class MatrixType, class InverseType>
void SparseContainer<MatrixType, InverseType>::
apply (ConstHostView X,
       HostView Y,
       int blockIndex,
       Teuchos::ETransp mode,
       SC alpha,
       SC beta) const
{
  using Teuchos::ArrayView;
  using Teuchos::as;

  // The InverseType template parameter might have different template
  // parameters (Scalar, LO, GO, and/or Node) than MatrixType.  For
  // example, MatrixType (a global object) might use a bigger GO
  // (global ordinal) type than InverseType (which deals with the
  // diagonal block, a local object).  This means that we might have
  // to convert X and Y to the Tpetra::MultiVector specialization that
  // InverseType wants.  This class' X_ and Y_ internal fields are of
  // the right type for InverseType, so we can use those as targets.
  Teuchos::RCP<Teuchos::Time> timer =
    Teuchos::TimeMonitor::getNewCounter ("Ifpack2::SparseContainer::apply");
  Teuchos::TimeMonitor timeMon (*timer);


  // Tpetra::MultiVector specialization corresponding to InverseType.
  Details::MultiVectorLocalGatherScatter<mv_type, inverse_mv_type> mvgs;
  size_t numVecs = X.extent(1);

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! this->IsComputed_, std::runtime_error, "Ifpack2::SparseContainer::apply: "
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

  const LO numRows = this->blockSizes_[blockIndex];

  // The operator Inverse_ works on a permuted subset of the local
  // parts of X and Y.  The subset and permutation are defined by the
  // index array returned by getBlockRows().  If the permutation is
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
  // of the local Map corresponding to Y.getMap().  numRows here
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
    for(LO i = 0; i < this->numBlocks_; i++)
      invX.emplace_back(Inverses_[i]->getDomainMap(), numVecs);
    for(LO i = 0; i < this->numBlocks_; i++)
      invY.emplace_back(Inverses_[i]->getDomainMap(), numVecs);
  }
  inverse_mv_type& X_local = invX[blockIndex];
  TEUCHOS_TEST_FOR_EXCEPTION(
    X_local.getLocalLength() != size_t(numRows * this->scalarsPerRow_), std::logic_error,
    "Ifpack2::SparseContainer::apply: "
    "X_local has length " << X_local.getLocalLength() << ", which does "
    "not match numRows = " << numRows * this->scalarsPerRow_ << ".  Please report this bug to "
    "the Ifpack2 developers.");
  const ArrayView<const LO> blockRows = this->getBlockRows(blockIndex);
  if(this->scalarsPerRow_ == 1)
    mvgs.gatherMVtoView(X_local, X, blockRows);
  else
    mvgs.gatherMVtoViewBlock(X_local, X, blockRows, this->scalarsPerRow_);

  // We must gather the output multivector Y even on input to
  // Inverse_->apply(), since the Inverse_ operator might use it as an
  // initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  inverse_mv_type& Y_local = invY[blockIndex];
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local.getLocalLength () != size_t(numRows * this->scalarsPerRow_), std::logic_error,
    "Ifpack2::SparseContainer::apply: "
    "Y_local has length " << Y_local.getLocalLength () << ", which does "
    "not match numRows = " << numRows * this->scalarsPerRow_ << ".  Please report this bug to "
    "the Ifpack2 developers.");

  if(this->scalarsPerRow_ == 1)
    mvgs.gatherMVtoView(Y_local, Y, blockRows);
  else
    mvgs.gatherMVtoViewBlock(Y_local, Y, blockRows, this->scalarsPerRow_);

  // Apply the local operator:
  // Y_local := beta*Y_local + alpha*M^{-1}*X_local
  this->solveBlockMV(X_local, Y_local, blockIndex, mode,
      InverseScalar(alpha), InverseScalar(beta));


  // Scatter the permuted subset output vector Y_local back into the
  // original output multivector Y.
  if(this->scalarsPerRow_ == 1)
    mvgs.scatterMVtoView(Y, Y_local, blockRows);
  else
    mvgs.scatterMVtoViewBlock(Y, Y_local, blockRows, this->scalarsPerRow_);
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType, InverseType>::
weightedApply (ConstHostView X,
               HostView Y,
               ConstHostView D,
               int blockIndex,
               Teuchos::ETransp mode,
               SC alpha,
               SC beta) const
{
  using Teuchos::ArrayView;
  using Teuchos::Range1D;
  using std::cerr;
  using std::endl;
  typedef Teuchos::ScalarTraits<InverseScalar> STS;

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
    ! this->IsComputed_, std::runtime_error, "Ifpack2::SparseContainer::"
    "weightedApply: You must have called the compute() method before you may "
    "call apply().  You may call the apply() method as many times as you want "
    "after calling compute() once, but you must have called compute() at least "
    "once.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent(1) != Y.extent(1), std::runtime_error,
    "Ifpack2::SparseContainer::weightedApply: X and Y have different numbers "
    "of vectors.  X has " << X.extent(1) << ", but Y has "
    << Y.extent(1) << ".");

  //bmk 7-2019: BlockRelaxation already checked this, but if that changes...
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->scalarsPerRow_ > 1, std::logic_error, "Ifpack2::SparseContainer::weightedApply: "
    "Use of block rows isn't allowed in overlapping Jacobi (the only method that uses weightedApply");

  if (numVecs == 0) {
    return; // done! nothing to do
  }

  // The operator Inverse_ works on a permuted subset of the local
  // parts of X and Y.  The subset and permutation are defined by the
  // index array returned by getBlockRows().  If the permutation is
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
  // of the local Map corresponding to Y.getMap().  numRows here
  // gives the number of rows in the row Map of the local Inverse_
  // operator.
  //
  // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
  // here that the row Map and the range Map of that operator are
  // the same.
  //
  // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
  // really belongs in Tpetra.
  const LO numRows = this->blockSizes_[blockIndex];

  if(invX.size() == 0)
  {
    for(LO i = 0; i < this->numBlocks_; i++)
      invX.emplace_back(Inverses_[i]->getDomainMap(), numVecs);
    for(LO i = 0; i < this->numBlocks_; i++)
      invY.emplace_back(Inverses_[i]->getDomainMap(), numVecs);
  }
  inverse_mv_type& X_local = invX[blockIndex];
  const ArrayView<const LO> blockRows = this->getBlockRows(blockIndex);
  mvgs.gatherMVtoView(X_local, X, blockRows);

  // We must gather the output multivector Y even on input to
  // Inverse_->apply(), since the Inverse_ operator might use it as an
  // initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  inverse_mv_type Y_local = invY[blockIndex];
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y_local.getLocalLength() != size_t(numRows), std::logic_error,
    "Ifpack2::SparseContainer::weightedApply: "
    "Y_local has length " << X_local.getLocalLength() << ", which does "
    "not match numRows = " << numRows << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gatherMVtoView(Y_local, Y, blockRows);

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
    D_local.getLocalLength() != size_t(this->blockSizes_[blockIndex]), std::logic_error,
    "Ifpack2::SparseContainer::weightedApply: "
    "D_local has length " << X_local.getLocalLength () << ", which does "
    "not match numRows = " << this->blockSizes_[blockIndex] << ".  Please report this bug to "
    "the Ifpack2 developers.");
  mvgs.gatherMVtoView(D_local, D, blockRows);
  inverse_mv_type X_scaled(Inverses_[blockIndex]->getDomainMap(), numVecs);
  X_scaled.elementWiseMultiply(STS::one(), D_local, X_local, STS::zero());

  // Y_temp will hold the result of M^{-1}*X_scaled.  If beta == 0, we
  // can write the result of Inverse_->apply() directly to Y_local, so
  // Y_temp may alias Y_local.  Otherwise, if beta != 0, we need
  // temporary storage for M^{-1}*X_scaled, so Y_temp must be
  // different than Y_local.
  inverse_mv_type* Y_temp;
  if (InverseScalar(beta) == STS::zero ()) {
    Y_temp = &Y_local;
  } else {
    Y_temp = new inverse_mv_type(Inverses_[blockIndex]->getRangeMap(), numVecs);
  }
  // Apply the local operator: Y_temp := M^{-1} * X_scaled
  Inverses_[blockIndex]->apply(X_scaled, *Y_temp, mode);
  // Y_local := beta * Y_local + alpha * diag(D_local) * Y_tmp.
  //
  // Note that we still use the permuted subset scaling D_local here,
  // because Y_temp has the same permuted subset Map.  That's good, in
  // fact, because it's a subset: less data to read and multiply.
  Y_local.elementWiseMultiply(alpha, D_local, *Y_temp, beta);
  if(Y_temp != &Y_local)
    delete Y_temp;
  // Copy the permuted subset output vector Y_local into the original
  // output multivector Y.
  mvgs.scatterMVtoView(Y, Y_local, blockRows);
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
  if (this->isInitialized()) {
    if (this->isComputed()) {
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
    os << "Block " << i << " rows: = " << this->blockSizes_[i] << endl;
  }
  os << "isInitialized()         = " << this->IsInitialized_ << endl;
  os << "isComputed()            = " << this->IsComputed_ << endl;
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
  using Teuchos::RCP;
  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  //To extract diagonal blocks, need to translate local rows to local columns.
  //Strategy: make a lookup table that translates local cols in the matrix to offsets in blockRows_:
  //blockOffsets_[b] <= offset < blockOffsets_[b+1]: tests whether the column is in block b.
  //offset - blockOffsets_[b]: gives the column within block b.
  //
  //This provides the block and col within a block in O(1).
  Teuchos::Array<InverseGlobalOrdinal> indicesToInsert;
  Teuchos::Array<InverseScalar> valuesToInsert;
  if(this->scalarsPerRow_ > 1)
  {
    Array<LO> colToBlockOffset(this->inputBlockMatrix_->getLocalNumCols(), INVALID);
    for(int i = 0; i < this->numBlocks_; i++)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockSize = this->blockSizes_[i];
      LO blockPointSize = this->blockSizes_[i] * this->scalarsPerRow_;
      LO blockEnd = blockStart + blockSize;
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(LO j = 0; j < blockSize; j++)
      {
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      //First, count the number of entries in each row of block i
      //(to allocate it with StaticProfile)
      Array<size_t> rowEntryCounts(blockPointSize, 0);
      //blockRow counts the BlockCrs LIDs that are going into this block
      //Rows are inserted into the CrsMatrix in sequential order
      using inds_type = typename block_crs_matrix_type::local_inds_host_view_type;
      using vals_type = typename block_crs_matrix_type::values_host_view_type;
      for(LO blockRow = 0; blockRow < blockSize; blockRow++)
      {
        //get a raw view of the whole block row
        inds_type indices;
        vals_type values;
        LO inputRow = this->blockRows_[blockStart + blockRow];
        this->inputBlockMatrix_->getLocalRowView(inputRow, indices, values);
        LO numEntries = (LO) indices.size();
        for(LO br = 0; br < this->bcrsBlockSize_; br++)
        {
          for(LO k = 0; k < numEntries; k++)
          {
            LO colOffset = colToBlockOffset[indices[k]];
            if(blockStart <= colOffset && colOffset < blockEnd)
            {
              rowEntryCounts[blockRow * this->bcrsBlockSize_ + br] += this->bcrsBlockSize_;
            }
          }
        }
      }
      //now that row sizes are known, can allocate the diagonal matrix
      RCP<InverseMap> tempMap(new InverseMap(blockPointSize, 0, this->localComm_));
      diagBlocks_[i] = rcp(new InverseCrs(tempMap, rowEntryCounts));
      Inverses_[i] = rcp(new InverseType(diagBlocks_[i]));
      //insert the actual entries, one row at a time
      for(LO blockRow = 0; blockRow < blockSize; blockRow++)
      {
        //get a raw view of the whole block row
        inds_type indices;
        vals_type values;
        LO inputRow = this->blockRows_[blockStart + blockRow];
        this->inputBlockMatrix_->getLocalRowView(inputRow, indices, values);
        LO numEntries = (LO) indices.size();
        for(LO br = 0; br < this->bcrsBlockSize_; br++)
        {
          indicesToInsert.clear();
          valuesToInsert.clear();
          for(LO k = 0; k < numEntries; k++)
          {
            LO colOffset = colToBlockOffset[indices[k]];
            if(blockStart <= colOffset && colOffset < blockEnd)
            {
              LO blockCol = colOffset - blockStart;
              //bc iterates over the columns in (block) entry k
              for(LO bc = 0; bc < this->bcrsBlockSize_; bc++)
              {
                indicesToInsert.push_back(blockCol * this->bcrsBlockSize_ + bc);
                valuesToInsert.push_back(values[k * this->bcrsBlockSize_ * this->bcrsBlockSize_ + bc * this->bcrsBlockSize_ + br]);
              }
            }
          }
          InverseGlobalOrdinal rowToInsert = blockRow * this->bcrsBlockSize_ + br;
          if(indicesToInsert.size())
            diagBlocks_[i]->insertGlobalValues(rowToInsert, indicesToInsert(), valuesToInsert());
        }
      }
      diagBlocks_[i]->fillComplete();
    }
  }
  else
  {
    //get the mapping from point-indexed matrix columns to offsets in blockRows_
    //(this includes regular CrsMatrix columns, in which case scalarsPerRow_ == 1)
    Array<LO> colToBlockOffset(this->inputMatrix_->getLocalNumCols() * this->bcrsBlockSize_, INVALID);
    for(int i = 0; i < this->numBlocks_; i++)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockSize = this->blockSizes_[i];
      LO blockEnd = blockStart + blockSize;
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(LO j = 0; j < blockSize; j++)
      {
        //translateRowToCol will return the corresponding split column
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      Teuchos::Array<size_t> rowEntryCounts(blockSize, 0);
      for(LO j = 0; j < blockSize; j++)
      {
        rowEntryCounts[j] = this->getInputRowView(this->blockRows_[blockStart + j]).size();
      }
      RCP<InverseMap> tempMap(new InverseMap(blockSize, 0, this->localComm_));
      diagBlocks_[i] = rcp(new InverseCrs(tempMap, rowEntryCounts));
      Inverses_[i] = rcp(new InverseType(diagBlocks_[i]));
      for(LO blockRow = 0; blockRow < blockSize; blockRow++)
      {
        valuesToInsert.clear();
        indicesToInsert.clear();
        //get a view of the split row
        LO inputSplitRow = this->blockRows_[blockStart + blockRow];
        auto rowView = this->getInputRowView(inputSplitRow);
        for(size_t k = 0; k < rowView.size(); k++)
        {
          LO colOffset = colToBlockOffset[rowView.ind(k)];
          if(blockStart <= colOffset && colOffset < blockEnd)
          {
            LO blockCol = colOffset - blockStart;
            indicesToInsert.push_back(blockCol);
            valuesToInsert.push_back(rowView.val(k));
          }
        }
        if(indicesToInsert.size())
          diagBlocks_[i]->insertGlobalValues(blockRow, indicesToInsert(), valuesToInsert());
      }
      diagBlocks_[i]->fillComplete ();
    }
  }
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
extractGraph ()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  //To extract diagonal blocks, need to translate local rows to local columns.
  //Strategy: make a lookup table that translates local cols in the matrix to offsets in blockRows_:
  //blockOffsets_[b] <= offset < blockOffsets_[b+1]: tests whether the column is in block b.
  //offset - blockOffsets_[b]: gives the column within block b.
  //
  //This provides the block and col within a block in O(1).
  Teuchos::RCP<Teuchos::Time> timer =
    Teuchos::TimeMonitor::getNewCounter ("Ifpack2::SparseContainer::extractGraph");
  Teuchos::TimeMonitor timeMon (*timer);

  Teuchos::Array<InverseGlobalOrdinal> indicesToInsert;
  if(this->scalarsPerRow_ > 1)
  {
    Array<LO> colToBlockOffset(this->inputBlockMatrix_->getLocalNumCols(), INVALID);
    for(int i = 0; i < this->numBlocks_; i++)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockSize = this->blockSizes_[i];
      LO blockPointSize = this->blockSizes_[i] * this->scalarsPerRow_;
      LO blockEnd = blockStart + blockSize;
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(LO j = 0; j < blockSize; j++)
      {
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      //First, count the number of entries in each row of block i
      //(to allocate it with StaticProfile)
      Array<size_t> rowEntryCounts(blockPointSize, 0);
      //blockRow counts the BlockCrs LIDs that are going into this block
      //Rows are inserted into the CrsMatrix in sequential order
      using inds_type = typename block_crs_matrix_type::local_inds_host_view_type;
      for(LO blockRow = 0; blockRow < blockSize; blockRow++)
      {
        //get a raw view of the whole block row
        inds_type indices;
        LO inputRow = this->blockRows_[blockStart + blockRow];
        this->inputBlockMatrix_->getGraph()->getLocalRowView(inputRow, indices);
        LO numEntries = (LO) indices.size();
        for(LO br = 0; br < this->bcrsBlockSize_; br++)
        {
          for(LO k = 0; k < numEntries; k++)
          {
            LO colOffset = colToBlockOffset[indices[k]];
            if(blockStart <= colOffset && colOffset < blockEnd)
            {
              rowEntryCounts[blockRow * this->bcrsBlockSize_ + br] += this->bcrsBlockSize_;
            }
          }
        }
      }
      //now that row sizes are known, can allocate the diagonal matrix
      RCP<InverseMap> tempMap(new InverseMap(blockPointSize, 0, this->localComm_));
      auto diagGraph = rcp(new InverseGraph(tempMap, rowEntryCounts));
      //insert the actual entries, one row at a time
      for(LO blockRow = 0; blockRow < blockSize; blockRow++)
      {
        //get a raw view of the whole block row
        inds_type indices;
        LO inputRow = this->blockRows_[blockStart + blockRow];
        this->inputBlockMatrix_->getGraph()->getLocalRowView(inputRow, indices);
        LO numEntries = (LO) indices.size();
        for(LO br = 0; br < this->bcrsBlockSize_; br++)
        {
          indicesToInsert.clear();
          for(LO k = 0; k < numEntries; k++)
          {
            LO colOffset = colToBlockOffset[indices[k]];
            if(blockStart <= colOffset && colOffset < blockEnd)
            {
              LO blockCol = colOffset - blockStart;
              //bc iterates over the columns in (block) entry k
              for(LO bc = 0; bc < this->bcrsBlockSize_; bc++)
              {
                indicesToInsert.push_back(blockCol * this->bcrsBlockSize_ + bc);
              }
            }
          }
          InverseGlobalOrdinal rowToInsert = blockRow * this->bcrsBlockSize_ + br;
          if(indicesToInsert.size())
            diagGraph->insertGlobalIndices(rowToInsert, indicesToInsert());
        }
      }
      diagGraph->fillComplete();

      // create matrix block
      diagBlocks_[i] = rcp(new InverseCrs(diagGraph));
      Inverses_[i] = rcp(new InverseType(diagBlocks_[i]));
    }
  }
  else
  {
    //get the mapping from point-indexed matrix columns to offsets in blockRows_
    //(this includes regular CrsMatrix columns, in which case scalarsPerRow_ == 1)
    Array<LO> colToBlockOffset(this->inputMatrix_->getLocalNumCols() * this->bcrsBlockSize_, INVALID);
    for(int i = 0; i < this->numBlocks_; i++)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockSize = this->blockSizes_[i];
      LO blockEnd = blockStart + blockSize;
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(LO j = 0; j < blockSize; j++)
      {
        //translateRowToCol will return the corresponding split column
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      Teuchos::Array<size_t> rowEntryCounts(blockSize, 0);
      for(LO j = 0; j < blockSize; j++)
      {
        rowEntryCounts[j] = this->getInputRowView(this->blockRows_[blockStart + j]).size();
      }
      RCP<InverseMap> tempMap(new InverseMap(blockSize, 0, this->localComm_));
      auto diagGraph = rcp(new InverseGraph(tempMap, rowEntryCounts));
      for(LO blockRow = 0; blockRow < blockSize; blockRow++)
      {
        indicesToInsert.clear();
        //get a view of the split row
        LO inputSplitRow = this->blockRows_[blockStart + blockRow];
        auto rowView = this->getInputRowView(inputSplitRow);
        for(size_t k = 0; k < rowView.size(); k++)
        {
          LO colOffset = colToBlockOffset[rowView.ind(k)];
          if(blockStart <= colOffset && colOffset < blockEnd)
          {
            LO blockCol = colOffset - blockStart;
            indicesToInsert.push_back(blockCol);
          }
        }
        if(indicesToInsert.size())
          diagGraph->insertGlobalIndices(blockRow, indicesToInsert());
      }
      diagGraph->fillComplete();

      // create matrix block
      diagBlocks_[i] = rcp(new InverseCrs(diagGraph));
      Inverses_[i] = rcp(new InverseType(diagBlocks_[i]));
    }
  }
}

//==============================================================================
template<class MatrixType, class InverseType>
void SparseContainer<MatrixType,InverseType>::
extractValues ()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  //To extract diagonal blocks, need to translate local rows to local columns.
  //Strategy: make a lookup table that translates local cols in the matrix to offsets in blockRows_:
  //blockOffsets_[b] <= offset < blockOffsets_[b+1]: tests whether the column is in block b.
  //offset - blockOffsets_[b]: gives the column within block b.
  //
  //This provides the block and col within a block in O(1).
  Teuchos::RCP<Teuchos::Time> timer =
    Teuchos::TimeMonitor::getNewCounter ("Ifpack2::SparseContainer::extractValues");
  Teuchos::TimeMonitor timeMon (*timer);

  Teuchos::Array<InverseGlobalOrdinal> indicesToInsert;
  Teuchos::Array<InverseScalar> valuesToInsert;
  if(this->scalarsPerRow_ > 1)
  {
    Array<LO> colToBlockOffset(this->inputBlockMatrix_->getLocalNumCols(), INVALID);
    for(int i = 0; i < this->numBlocks_; i++)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockSize = this->blockSizes_[i];
      LO blockEnd = blockStart + blockSize;
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(LO j = 0; j < blockSize; j++)
      {
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      using inds_type = typename block_crs_matrix_type::local_inds_host_view_type;
      using vals_type = typename block_crs_matrix_type::values_host_view_type;
      //insert the actual entries, one row at a time
      diagBlocks_[i]->resumeFill();
      for(LO blockRow = 0; blockRow < blockSize; blockRow++)
      {
        //get a raw view of the whole block row
        inds_type indices;
        vals_type values;
        LO inputRow = this->blockRows_[blockStart + blockRow];
        this->inputBlockMatrix_->getLocalRowView(inputRow, indices, values);
        LO numEntries = (LO) indices.size();
        for(LO br = 0; br < this->bcrsBlockSize_; br++)
        {
          indicesToInsert.clear();
          valuesToInsert.clear();
          for(LO k = 0; k < numEntries; k++)
          {
            LO colOffset = colToBlockOffset[indices[k]];
            if(blockStart <= colOffset && colOffset < blockEnd)
            {
              LO blockCol = colOffset - blockStart;
              //bc iterates over the columns in (block) entry k
              for(LO bc = 0; bc < this->bcrsBlockSize_; bc++)
              {
                indicesToInsert.push_back(blockCol * this->bcrsBlockSize_ + bc);
                valuesToInsert.push_back(values[k * this->bcrsBlockSize_ * this->bcrsBlockSize_ + bc * this->bcrsBlockSize_ + br]);
              }
            }
          }
          InverseGlobalOrdinal rowToInsert = blockRow * this->bcrsBlockSize_ + br;
          if(indicesToInsert.size())
            diagBlocks_[i]->replaceGlobalValues(rowToInsert, indicesToInsert(), valuesToInsert());
        }
      }
      diagBlocks_[i]->fillComplete();
    }
  }
  else
  {
    //get the mapping from point-indexed matrix columns to offsets in blockRows_
    //(this includes regular CrsMatrix columns, in which case scalarsPerRow_ == 1)
    Array<LO> colToBlockOffset(this->inputMatrix_->getLocalNumCols() * this->bcrsBlockSize_, INVALID);
    for(int i = 0; i < this->numBlocks_; i++)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockSize = this->blockSizes_[i];
      LO blockEnd = blockStart + blockSize;
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(LO j = 0; j < blockSize; j++)
      {
        //translateRowToCol will return the corresponding split column
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      diagBlocks_[i]->resumeFill();
      for(LO blockRow = 0; blockRow < blockSize; blockRow++)
      {
        valuesToInsert.clear();
        indicesToInsert.clear();
        //get a view of the split row
        LO inputSplitRow = this->blockRows_[blockStart + blockRow];
        auto rowView = this->getInputRowView(inputSplitRow);
        for(size_t k = 0; k < rowView.size(); k++)
        {
          LO colOffset = colToBlockOffset[rowView.ind(k)];
          if(blockStart <= colOffset && colOffset < blockEnd)
          {
            LO blockCol = colOffset - blockStart;
            indicesToInsert.push_back(blockCol);
            valuesToInsert.push_back(rowView.val(k));
          }
        }
        if(indicesToInsert.size())
          diagBlocks_[i]->replaceGlobalValues(blockRow, indicesToInsert, valuesToInsert);
      }
      diagBlocks_[i]->fillComplete ();
    }
  }
}

template<typename MatrixType, typename InverseType>
std::string SparseContainer<MatrixType, InverseType>::getName()
{
  typedef ILUT<Tpetra::RowMatrix<SC, LO, GO, NO> > ILUTInverse;
#ifdef HAVE_IFPACK2_AMESOS2
  typedef Details::Amesos2Wrapper<Tpetra::RowMatrix<SC, LO, GO, NO>> AmesosInverse;
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
