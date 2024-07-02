// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_CONTAINER_DEF_HPP
#define IFPACK2_CONTAINER_DEF_HPP

#include <Ifpack2_Details_MultiVectorLocalGatherScatter.hpp>
#include <Teuchos_Time.hpp>

namespace Ifpack2 {

//Implementation of Ifpack2::Container

template<class MatrixType>
Container<MatrixType>::Container(
    const Teuchos::RCP<const row_matrix_type>& matrix,
    const Teuchos::Array<Teuchos::Array<LO> >& partitions,
    bool pointIndexed) :
  inputMatrix_ (matrix),
  inputCrsMatrix_ (Teuchos::rcp_dynamic_cast<const crs_matrix_type>(inputMatrix_)),
  inputBlockMatrix_ (Teuchos::rcp_dynamic_cast<const block_crs_matrix_type>(inputMatrix_)),
  pointIndexed_(pointIndexed),
  IsInitialized_(false),
  IsComputed_(false)
{
  using Teuchos::Ptr;
  using Teuchos::RCP;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  NumLocalRows_ = inputMatrix_->getLocalNumRows();
  NumGlobalRows_ = inputMatrix_->getGlobalNumRows();
  NumGlobalNonzeros_ = inputMatrix_->getGlobalNumEntries();
  IsParallel_ = inputMatrix_->getRowMap()->getComm()->getSize() != 1;
  hasBlockCrs_ = !inputBlockMatrix_.is_null();
  if(hasBlockCrs_)
    bcrsBlockSize_ = inputBlockMatrix_->getBlockSize();
  else
    bcrsBlockSize_ = 1;
  if(hasBlockCrs_ && !pointIndexed_)
    scalarsPerRow_ = bcrsBlockSize_;
  else
    scalarsPerRow_ = 1;
  setBlockSizes(partitions);
  //Sanity check the partitions
  #ifdef HAVE_IFPACK2_DEBUG
  // Check whether the input set of local row indices is correct.
  const map_type& rowMap = *inputMatrix_->getRowMap();
  for(int i = 0; i < numBlocks_; i++)
  {
    Teuchos::ArrayView<const LO> blockRows = getBlockRows(i);
    for(LO j = 0; j < blockSizes_[i]; j++)
    {
      LO row = blockRows[j];
      if(pointIndexed)
      {
        //convert the point row to the corresponding block row
        row /= bcrsBlockSize_;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(
          !rowMap.isNodeLocalElement(row),
          std::invalid_argument, "Ifpack2::Container: "
          "On process " << rowMap.getComm()->getRank() << " of "
          << rowMap.getComm()->getSize() << ", in the given set of local row "
          "indices blockRows = " << Teuchos::toString(blockRows) << ", the following "
          "entries is not valid local row index on the calling process: "
          << row << ".");
    }
  }
  #endif
}

template<class MatrixType>
Container<MatrixType>::
~Container() {}

template<class MatrixType>
Teuchos::ArrayView<const typename MatrixType::local_ordinal_type>
Container<MatrixType>::getBlockRows(int blockIndex) const
{
  return Teuchos::ArrayView<const LO>
    (&blockRows_[blockOffsets_[blockIndex]], blockSizes_[blockIndex]);
}

template<class MatrixType>
void Container<MatrixType>::setBlockSizes(const Teuchos::Array<Teuchos::Array<LO> >& partitions)
{
  //First, create a grand total of all the rows in all the blocks
  //Note: If partitioner allowed overlap, this could be greater than the # of local rows
  LO totalBlockRows = 0;
  numBlocks_ = partitions.size();
  blockSizes_.resize(numBlocks_);
  blockOffsets_.resize(numBlocks_);
  maxBlockSize_ = 0;
  for(int i = 0; i < numBlocks_; i++)
  {
    LO rowsInBlock = partitions[i].size();
    blockSizes_[i] = rowsInBlock;
    blockOffsets_[i] = totalBlockRows;
    totalBlockRows += rowsInBlock;
    maxBlockSize_ = std::max(maxBlockSize_, rowsInBlock * scalarsPerRow_);
  }
  blockRows_.resize(totalBlockRows);
  //set blockRows_: each entry is the partition/block of the row
  LO iter = 0;
  for(int i = 0; i < numBlocks_; i++)
  {
    for(int j = 0; j < blockSizes_[i]; j++)
    {
      blockRows_[iter++] = partitions[i][j];
    }
  }
}

template<class MatrixType>
void Container<MatrixType>::getMatDiag() const
{
  if(Diag_.is_null())
  {
    Diag_ = rcp(new vector_type(inputMatrix_->getDomainMap()));
    inputMatrix_->getLocalDiagCopy(*Diag_);
  }
}

template<class MatrixType>
bool Container<MatrixType>::isInitialized() const {
  return IsInitialized_;
}

template<class MatrixType>
bool Container<MatrixType>::isComputed () const {
  return IsComputed_;
}

template<class MatrixType>
void Container<MatrixType>::
applyMV(const mv_type& X, mv_type& Y) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
}

template<class MatrixType>
void Container<MatrixType>::
weightedApplyMV(const mv_type& X, mv_type& Y, vector_type& W) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
}

template<class MatrixType>
std::string Container<MatrixType>::
getName()
{
  return "Generic";
}

template<class MatrixType>
void Container<MatrixType>::DoGSBlock(ConstHostView X, HostView Y, HostView Y2, HostView Resid,
    SC dampingFactor, LO i) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
}

template <class MatrixType>
void Container<MatrixType>::DoJacobi(ConstHostView X, HostView Y, SC dampingFactor) const
{
  using STS = Teuchos::ScalarTraits<ISC>;
  const ISC one = STS::one();
  // use blockRows_ and blockSizes_
  size_t numVecs = X.extent(1);
  // Non-overlapping Jacobi
  for (LO i = 0; i < numBlocks_; i++)
  {
    // may happen that a partition is empty
    if(blockSizes_[i] != 1 || hasBlockCrs_)
    {
      if(blockSizes_[i] == 0 )
        continue;
      apply(X, Y, i, Teuchos::NO_TRANS, dampingFactor, one);
    }
    else    // singleton, can't access Containers_[i] as it was never filled and may be null.
    {
      LO LRID = blockRows_[blockOffsets_[i]];
      getMatDiag();
      auto diagView = Diag_->getLocalViewHost(Tpetra::Access::ReadOnly);
      ISC d = one * (static_cast<ISC> (dampingFactor)) / diagView(LRID, 0);
      for(size_t nv = 0; nv < numVecs; nv++)
      {
        ISC x = X(LRID, nv);
        Y(LRID, nv) += x * d;
      }
    }
  }
}

template <class MatrixType>
void Container<MatrixType>::DoOverlappingJacobi(ConstHostView X, HostView Y, ConstHostView W, SC dampingFactor, bool nonsymScaling) const
{
  using STS = Teuchos::ScalarTraits<SC>;
  // Overlapping Jacobi
  for(LO i = 0; i < numBlocks_; i++)
  {
    // may happen that a partition is empty
    if(blockSizes_[i] == 0)
      continue;
    if(blockSizes_[i] != 1) {
      if (!nonsymScaling) 
        weightedApply(X, Y, W, i, Teuchos::NO_TRANS, dampingFactor, STS::one());
      else {
        // A crummy way of doing nonsymmetric scaling. We effectively
        // first reverse scale x, which later gets scaled inside weightedApply
        // so the net effect is that x is not scaled.
        // This was done to keep using weightedApply() that is defined in
        // many spots in the code.
        HostView tempo("", X.extent(0), X.extent(1));
        size_t numVecs = X.extent(1);
        LO  bOffset = blockOffsets_[i];
        for (LO ii = 0; ii < blockSizes_[i]; ii++) {
          LO LRID = blockRows_[bOffset++];
          for (size_t jj = 0; jj < numVecs; jj++) tempo(LRID,jj)=X(LRID,jj)/ W(LRID,0);
        }
        weightedApply(tempo, Y, W, i, Teuchos::NO_TRANS, dampingFactor, STS::one());
      }
    }
    else    // singleton, can't access Containers_[i] as it was never filled and may be null.
    {
      const ISC one = STS::one();
      size_t numVecs = X.extent(1);
      LO LRID = blockRows_[blockOffsets_[i]];
      getMatDiag();
      auto diagView = Diag_->getLocalViewHost(Tpetra::Access::ReadOnly);
      ISC recip  = one / diagView(LRID, 0);
      ISC wval   = W(LRID,0);
      ISC combo  = wval*recip;
      ISC d = combo*(static_cast<ISC> (dampingFactor));
      for(size_t nv = 0; nv < numVecs; nv++)
      {
        ISC x = X(LRID, nv);
        Y(LRID, nv) = x * d + Y(LRID, nv);
      }
    }
  }
}

//Do Gauss-Seidel with just block i
//This is used 3 times: once in DoGaussSeidel and twice in DoSGS
template<class MatrixType, typename LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::DoGSBlock(
    ConstHostView X, HostView Y, HostView Y2, HostView Resid,
    SC dampingFactor, LO i) const
{
  using Teuchos::ArrayView;
  using STS = Teuchos::ScalarTraits<ISC>;
  size_t numVecs = X.extent(1);
  const ISC one = STS::one();
  if(this->blockSizes_[i] == 0)
    return; // Skip empty partitions
  if(this->hasBlockCrs_ && !this->pointIndexed_)
  {
    //Use efficient blocked version
    ArrayView<const LO> blockRows = this->getBlockRows(i);
    const size_t localNumRows = this->blockSizes_[i];
    using inds_type = typename block_crs_matrix_type::local_inds_host_view_type;
    using vals_type = typename block_crs_matrix_type::values_host_view_type;
    for(size_t j = 0; j < localNumRows; j++)
    {
      LO row = blockRows[j]; // Containers_[i]->ID (j);
      vals_type values;
      inds_type colinds;
      this->inputBlockMatrix_->getLocalRowView(row, colinds, values);
      LO numEntries = (LO) colinds.size();
      for(size_t m = 0; m < numVecs; m++)
      {
        for (int localR = 0; localR < this->bcrsBlockSize_; localR++)
          Resid(row * this->bcrsBlockSize_ + localR, m) = X(row * this->bcrsBlockSize_ + localR, m);
        for (LO k = 0; k < numEntries; ++k)
        {
          const LO col = colinds[k];
          for(int localR = 0; localR < this->bcrsBlockSize_; localR++)
          {
            for(int localC = 0; localC < this->bcrsBlockSize_; localC++)
            {
              Resid(row * this->bcrsBlockSize_ + localR, m) -=
                values[k * this->bcrsBlockSize_ * this->bcrsBlockSize_ + localR + localC * this->bcrsBlockSize_]
                * Y2(col * this->bcrsBlockSize_ + localC, m); }
          }
        }
      }
    }
    // solve with this block
    //
    // Note: I'm abusing the ordering information, knowing that X/Y
    // and Y2 have the same ordering for on-proc unknowns.
    //
    // Note: Add flop counts for inverse apply
    apply(Resid, Y2, i, Teuchos::NO_TRANS, dampingFactor, one);
  }
  else if(!this->hasBlockCrs_ && this->blockSizes_[i] == 1)
  {
    // singleton, can't access Containers_[i] as it was never filled and may be null.
    // a singleton calculation (just using matrix diagonal) is exact, all residuals should be zero.
    LO LRID = this->blockOffsets_[i];  // by definition, a singleton 1 row in block.
    //Use the KokkosSparse internal matrix for low-overhead values/indices access
    //But, can only do this if the matrix is accessible directly from host, since it's not a DualView
    using container_exec_space = typename ContainerImpl<MatrixType, LocalScalarType>::crs_matrix_type::execution_space;
    container_exec_space().fence();
    auto localA = this->inputCrsMatrix_->getLocalMatrixHost();
    using size_type = typename crs_matrix_type::local_matrix_host_type::size_type;
    const auto& rowmap = localA.graph.row_map;
    const auto& entries = localA.graph.entries;
    const auto& values = localA.values;
    this->getMatDiag();
    auto diagView = this->Diag_->getLocalViewHost(Tpetra::Access::ReadOnly);
    ISC d = (static_cast<ISC> (dampingFactor)) / diagView(LRID, 0);
    for(size_t m = 0; m < numVecs; m++)
    {
      // ISC x = X(LRID, m);
      ISC r = X(LRID, m);
      for(size_type k = rowmap(LRID); k < rowmap(LRID + 1); k++) {
        const LO col = entries(k);
        r -= values(k) * Y2(col, m);
      }

      ISC newy = r * d;
      Y2(LRID, m) += newy;
    }
  }
  else if(!this->inputCrsMatrix_.is_null() &&
      std::is_same<typename crs_matrix_type::device_type::memory_space, Kokkos::HostSpace>::value)
  {
    //Use the KokkosSparse internal matrix for low-overhead values/indices access
    //But, can only do this if the matrix is accessible directly from host, since it's not a DualView
    using container_exec_space = typename ContainerImpl<MatrixType, LocalScalarType>::crs_matrix_type::execution_space;
    container_exec_space().fence();
    auto localA = this->inputCrsMatrix_->getLocalMatrixHost();
    using size_type = typename crs_matrix_type::local_matrix_host_type::size_type;
    const auto& rowmap = localA.graph.row_map;
    const auto& entries = localA.graph.entries;
    const auto& values = localA.values;
    ArrayView<const LO> blockRows = this->getBlockRows(i);
    for(size_t j = 0; j < size_t(blockRows.size()); j++)
    {
      const LO row = blockRows[j];
      for(size_t m = 0; m < numVecs; m++)
      {
        ISC r = X(row, m);
        for(size_type k = rowmap(row); k < rowmap(row + 1); k++)
        {
          const LO col = entries(k);
          r -= values(k) * Y2(col, m);
        }
        Resid(row, m) = r;
      }
    }
    // solve with this block
    //
    // Note: I'm abusing the ordering information, knowing that X/Y
    // and Y2 have the same ordering for on-proc unknowns.
    //
    // Note: Add flop counts for inverse apply
    apply(Resid, Y2, i, Teuchos::NO_TRANS, dampingFactor, one);
  }
  else
  {
    //Either a point-indexed block matrix, or a normal row matrix
    //that doesn't support getLocalMatrix
    ArrayView<const LO> blockRows = this->getBlockRows(i);
    for(size_t j = 0; j < size_t(blockRows.size()); j++)
    {
      const LO row = blockRows[j];
      auto rowView = getInputRowView(row);
      for(size_t m = 0; m < numVecs; m++)
      {
        Resid(row, m) = X(row, m);
        for (size_t k = 0; k < rowView.size(); ++k)
        {
          const LO col = rowView.ind(k);
          Resid(row, m) -= rowView.val(k) * Y2(col, m);
        }
      }
    }
    // solve with this block
    //
    // Note: I'm abusing the ordering information, knowing that X/Y
    // and Y2 have the same ordering for on-proc unknowns.
    //
    // Note: Add flop counts for inverse apply
    apply(Resid, Y2, i, Teuchos::NO_TRANS, dampingFactor, one);
  }
}

template<class MatrixType>
void Container<MatrixType>::
DoGaussSeidel(ConstHostView X, HostView Y, HostView Y2, SC dampingFactor) const
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  //This function just extracts the diagonal if it hasn't already.
  getMatDiag();
  auto numVecs = X.extent(1);
  // X = RHS, Y = initial guess
  HostView Resid("", X.extent(0), X.extent(1));
  for(LO i = 0; i < numBlocks_; i++)
  {
    DoGSBlock(X, Y, Y2, Resid, dampingFactor, i);
  }
  if(IsParallel_)
  {
    auto numMyRows = inputMatrix_->getLocalNumRows();
    for (size_t m = 0; m < numVecs; ++m)
    {
      for (size_t i = 0; i < numMyRows * bcrsBlockSize_; ++i)
      {
        Y(i, m) = Y2(i, m);
      }
    }
  }
}

template<class MatrixType>
void Container<MatrixType>::
DoSGS(ConstHostView X, HostView Y, HostView Y2, SC dampingFactor) const
{
  // X = RHS, Y = initial guess
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Ptr;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  auto numVecs = X.extent(1);
  HostView Resid("", X.extent(0), X.extent(1));
  // Forward Sweep
  for(LO i = 0; i < numBlocks_; i++)
  {
    DoGSBlock(X, Y, Y2, Resid, dampingFactor, i);
  }
  static_assert(std::is_signed<LO>::value,
      "Local ordinal must be signed (unsigned breaks reverse iteration to 0)");
  // Reverse Sweep
  for(LO i = numBlocks_ - 1; i >= 0; --i)
  {
    DoGSBlock(X, Y, Y2, Resid, dampingFactor, i);
  }
  if(IsParallel_)
  {
    auto numMyRows = inputMatrix_->getLocalNumRows();
    for (size_t m = 0; m < numVecs; ++m)
    {
      for (size_t i = 0; i < numMyRows * bcrsBlockSize_; ++i)
      {
        Y(i, m) = Y2(i, m);
      }
    }
  }
}

template<class MatrixType>
void Container<MatrixType>::
clearBlocks()
{
  numBlocks_ = 0;
  blockRows_.clear();
  blockSizes_.clear();
  blockOffsets_.clear();
  Diag_ = Teuchos::null;      //Diag_ will be recreated if needed
}

//Implementation of Ifpack2::ContainerImpl

template<class MatrixType, class LocalScalarType>
ContainerImpl<MatrixType, LocalScalarType>::
ContainerImpl(
      const Teuchos::RCP<const row_matrix_type>& matrix,
      const Teuchos::Array<Teuchos::Array<LO> >& partitions,
      bool pointIndexed)
  : Container<MatrixType>(matrix, partitions, pointIndexed) {}

template<class MatrixType, class LocalScalarType>
ContainerImpl<MatrixType, LocalScalarType>::
~ContainerImpl() {}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
setParameters (const Teuchos::ParameterList& List) {}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
applyInverseJacobi (const mv_type& /* X */, mv_type& /* Y */,
                         SC dampingFactor,
                         bool /* zeroStartingSolution = false */,
                         int /* numSweeps = 1 */) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
applyMV (const mv_type& X, mv_type& Y) const
{
  ConstHostView XView = X.getLocalViewHost(Tpetra::Access::ReadOnly);
  HostView YView = Y.getLocalViewHost(Tpetra::Access::ReadWrite);
  this->apply (XView, YView, 0);
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
weightedApplyMV (const mv_type& X,
                 mv_type& Y,
                 vector_type& W) const
{
  ConstHostView XView = X.getLocalViewHost(Tpetra::Access::ReadOnly);
  HostView YView = Y.getLocalViewHost(Tpetra::Access::ReadWrite);
  ConstHostView WView = W.getLocalViewHost(Tpetra::Access::ReadOnly);
  weightedApply (XView, YView, WView, 0);
}

template<class MatrixType, class LocalScalarType>
std::string ContainerImpl<MatrixType, LocalScalarType>::
getName()
{
  return "Generic";
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
solveBlock(ConstHostSubviewLocal X,
           HostSubviewLocal Y,
           int blockIndex,
           Teuchos::ETransp mode,
           const LSC alpha,
           const LSC beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Not implemented.");
}

template<class MatrixType, class LocalScalarType>
typename ContainerImpl<MatrixType, LocalScalarType>::LO
ContainerImpl<MatrixType, LocalScalarType>::
translateRowToCol(LO row)
{
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  GO GO_INVALID = Teuchos::OrdinalTraits<GO>::invalid();
  const map_type& globalRowMap = *(this->inputMatrix_->getRowMap());
  const map_type& globalColMap = *(this->inputMatrix_->getColMap());
  LO rowLID = row;
  LO dofOffset = 0;
  if(this->pointIndexed_)
  {
    rowLID = row / this->bcrsBlockSize_;
    dofOffset = row % this->bcrsBlockSize_;
  }
  GO diagGID = globalRowMap.getGlobalElement(rowLID);
  TEUCHOS_TEST_FOR_EXCEPTION(
    diagGID == GO_INVALID,
    std::runtime_error, "Ifpack2::Container::translateRowToCol: "
    "On process " << this->inputMatrix_->getRowMap()->getComm()->getRank() <<
    ", at least one row index in the set of local "
    "row indices given to the constructor is not a valid local row index in "
    "the input matrix's row Map on this process.  This should be impossible "
    "because the constructor checks for this case.  Here is the complete set "
    "of invalid local row indices: " << rowLID << ".  "
    "Please report this bug to the Ifpack2 developers.");
  //now, can translate diagGID (both a global row AND global col ID) to local column
  LO colLID = globalColMap.getLocalElement(diagGID);
  TEUCHOS_TEST_FOR_EXCEPTION(
    colLID == LO_INVALID,
    std::runtime_error, "Ifpack2::Container::translateRowToCol: "
    "On process " << this->inputMatrix_->getRowMap()->getComm()->getRank() << ", "
    "at least one row index in the set of row indices given to the constructor "
    "does not have a corresponding column index in the input matrix's column "
    "Map.  This probably means that the column(s) in question is/are empty on "
    "this process, which would make the submatrix to extract structurally "
    "singular. The invalid global column index is " << diagGID << ".");
  //colLID could identify a block column - translate to split column if needed
  if(this->pointIndexed_)
    return colLID * this->bcrsBlockSize_ + dofOffset;
  return colLID;
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
apply (ConstHostView X,
       HostView Y,
       int blockIndex,
       Teuchos::ETransp mode,
       SC alpha,
       SC beta) const
{
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // The local operator might have a different Scalar type than
  // MatrixType.  This means that we might have to convert X and Y to
  // the Tpetra::MultiVector specialization that the local operator
  // wants.  This class' X_ and Y_ internal fields are of the right
  // type for the local operator, so we can use those as targets.

  Details::MultiVectorLocalGatherScatter<mv_type, local_mv_type> mvgs;

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! this->IsComputed_, std::runtime_error, "Ifpack2::Container::apply: "
    "You must have called the compute() method before you may call apply().  "
    "You may call the apply() method as many times as you want after calling "
    "compute() once, but you must have called compute() at least once.");

  const size_t numVecs = X.extent(1);

  if(numVecs == 0) {
    return; // done! nothing to do
  }

  // The local operator works on a permuted subset of the local parts
  // of X and Y.  The subset and permutation are defined by the index
  // array returned by getBlockRows().  If the permutation is trivial
  // and the subset is exactly equal to the local indices, then we
  // could use the local parts of X and Y exactly, without needing to
  // permute.  Otherwise, we have to use temporary storage to permute
  // X and Y.  For now, we always use temporary storage.
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

  if(X_localBlocks_.size() == 0 || X.extent(1) != X_local_.extent(1))
  {
    //need to resize (or create for the first time) the three scratch arrays
    X_localBlocks_.clear();
    Y_localBlocks_.clear();
    X_localBlocks_.reserve(this->numBlocks_);
    Y_localBlocks_.reserve(this->numBlocks_);

    X_local_ = HostViewLocal("X_local", this->blockRows_.size() * this->scalarsPerRow_, numVecs);
    Y_local_ = HostViewLocal("Y_local", this->blockRows_.size() * this->scalarsPerRow_, numVecs);

    //create all X_local and Y_local managed Views at once, are
    //reused in subsequent apply() calls
    for(int i = 0; i < this->numBlocks_; i++)
    {
      auto blockBounds = std::make_pair(this->blockOffsets_[i] * this->scalarsPerRow_,
          (this->blockOffsets_[i] + this->blockSizes_[i]) * this->scalarsPerRow_);
      X_localBlocks_.emplace_back(X_local_, blockBounds, Kokkos::ALL());
      Y_localBlocks_.emplace_back(Y_local_, blockBounds, Kokkos::ALL());
    }
  }

  const ArrayView<const LO> blockRows = this->getBlockRows(blockIndex);

  if(this->scalarsPerRow_ == 1)
    mvgs.gatherViewToView (X_localBlocks_[blockIndex], X, blockRows);
  else
    mvgs.gatherViewToViewBlock (X_localBlocks_[blockIndex], X, blockRows, this->scalarsPerRow_);

  // We must gather the contents of the output multivector Y even on
  // input to solveBlock(), since the inverse operator might use it as
  // an initial guess for a linear solve.  We have no way of knowing
  // whether it does or does not.

  if(this->scalarsPerRow_ == 1)
    mvgs.gatherViewToView (Y_localBlocks_[blockIndex], Y, blockRows);
  else
    mvgs.gatherViewToViewBlock (Y_localBlocks_[blockIndex], Y, blockRows, this->scalarsPerRow_);

  // Apply the local operator:
  // Y_local := beta*Y_local + alpha*M^{-1}*X_local
  this->solveBlock (X_localBlocks_[blockIndex], Y_localBlocks_[blockIndex], blockIndex, mode,
                   LSC(alpha), LSC(beta));

  // Scatter the permuted subset output vector Y_local back into the
  // original output multivector Y.
  if(this->scalarsPerRow_ == 1)
    mvgs.scatterViewToView (Y, Y_localBlocks_[blockIndex], blockRows);
  else
    mvgs.scatterViewToViewBlock (Y, Y_localBlocks_[blockIndex], blockRows, this->scalarsPerRow_);
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
weightedApply(ConstHostView X,
              HostView Y,
              ConstHostView D,
              int blockIndex,
              Teuchos::ETransp mode,
              SC alpha,
              SC beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::Range1D;
  using Teuchos::Ptr;
  using Teuchos::ptr;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using std::endl;
  using STS = Teuchos::ScalarTraits<SC>;

  // The local operator template parameter might have a different
  // Scalar type than MatrixType.  This means that we might have to
  // convert X and Y to the Tpetra::MultiVector specialization that
  // the local operator wants.  This class' X_ and Y_ internal fields
  // are of the right type for the local operator, so we can use those
  // as targets.

  const char prefix[] = "Ifpack2::Container::weightedApply: ";
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! this->IsComputed_, std::runtime_error, prefix << "You must have called the "
    "compute() method before you may call this method.  You may call "
    "weightedApply() as many times as you want after calling compute() once, "
    "but you must have called compute() at least once first.");
  
  //bmk 7-2019: BlockRelaxation already checked this, but if that changes...
  TEUCHOS_TEST_FOR_EXCEPTION(
    this->scalarsPerRow_ > 1, std::logic_error, prefix << "Use of block rows isn't allowed "
    "in overlapping Jacobi (the only method that uses weightedApply");

  const size_t numVecs = X.extent(1);

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent(1) != Y.extent(1), std::runtime_error,
    prefix << "X and Y have different numbers of vectors (columns).  X has "
    << X.extent(1) << ", but Y has " << Y.extent(1) << ".");

  if(numVecs == 0) {
    return; // done! nothing to do
  }

  const size_t numRows = this->blockSizes_[blockIndex];

  // The local operator works on a permuted subset of the local parts
  // of X and Y.  The subset and permutation are defined by the index
  // array returned by getBlockRows().  If the permutation is trivial
  // and the subset is exactly equal to the local indices, then we
  // could use the local parts of X and Y exactly, without needing to
  // permute.  Otherwise, we have to use temporary storage to permute
  // X and Y.  For now, we always use temporary storage.
  //
  // Create temporary permuted versions of the input and output.
  // (Re)allocate X_ and/or Y_ only if necessary.  We'll use them to
  // store the permuted versions of X resp. Y.  Note that X_local has
  // the domain Map of the operator, which may be a permuted subset of
  // the local Map corresponding to X.getMap().  Similarly, Y_local
  // has the range Map of the operator, which may be a permuted subset
  // of the local Map corresponding to Y.getMap().  numRows_ here
  // gives the number of rows in the row Map of the local operator.
  //
  // FIXME (mfh 20 Aug 2013) There might be an implicit assumption
  // here that the row Map and the range Map of that operator are
  // the same.
  //
  // FIXME (mfh 20 Aug 2013) This "local permutation" functionality
  // really belongs in Tpetra.
  if(X_localBlocks_.size() == 0 || X.extent(1) != X_local_.extent(1))
  {
    //need to resize (or create for the first time) the three scratch arrays
    X_localBlocks_.clear();
    Y_localBlocks_.clear();
    X_localBlocks_.reserve(this->numBlocks_);
    Y_localBlocks_.reserve(this->numBlocks_);

    X_local_ = HostViewLocal("X_local", this->blockRows_.size() * this->scalarsPerRow_, numVecs);
    Y_local_ = HostViewLocal("Y_local", this->blockRows_.size() * this->scalarsPerRow_, numVecs);

    //create all X_local and Y_local managed Views at once, are
    //reused in subsequent apply() calls
    for(int i = 0; i < this->numBlocks_; i++)
    {
      auto blockBounds = std::make_pair(this->blockOffsets_[i] * this->scalarsPerRow_,
          (this->blockOffsets_[i] + this->blockSizes_[i]) * this->scalarsPerRow_);
      X_localBlocks_.emplace_back(X_local_, blockBounds, Kokkos::ALL());
      Y_localBlocks_.emplace_back(Y_local_, blockBounds, Kokkos::ALL());
    }
  }
  if(int(weightedApplyScratch_.extent(0)) != 3 * this->maxBlockSize_ ||
      weightedApplyScratch_.extent(1) != numVecs)
  {
    weightedApplyScratch_ = HostViewLocal("weightedApply scratch", 3 * this->maxBlockSize_, numVecs);
  }

  ArrayView<const LO> blockRows = this->getBlockRows(blockIndex);

  Details::MultiVectorLocalGatherScatter<mv_type, local_mv_type> mvgs;

  //note: BlockCrs w/ weighted Jacobi isn't allowed, so no need to use block gather/scatter
  mvgs.gatherViewToView (X_localBlocks_[blockIndex], X, blockRows);
  // We must gather the output multivector Y even on input to
  // solveBlock(), since the local operator might use it as an initial
  // guess for a linear solve.  We have no way of knowing whether it
  // does or does not.

  mvgs.gatherViewToView (Y_localBlocks_[blockIndex], Y, blockRows);

  // Apply the diagonal scaling D to the input X.  It's our choice
  // whether the result has the original input Map of X, or the
  // permuted subset Map of X_local.  If the latter, we also need to
  // gather D into the permuted subset Map.  We choose the latter, to
  // save memory and computation.  Thus, we do the following:
  //
  // 1. Gather D into a temporary vector D_local.
  // 2. Create a temporary X_scaled to hold diag(D_local) * X_local.
  // 3. Compute X_scaled := diag(D_loca) * X_local.
  auto maxBS = this->maxBlockSize_;
  auto bs = this->blockSizes_[blockIndex] * this->scalarsPerRow_;

  HostSubviewLocal D_local(weightedApplyScratch_, std::make_pair(0, bs), std::make_pair(0, 1));
  mvgs.gatherViewToView (D_local, D, blockRows);
  HostSubviewLocal X_scaled(weightedApplyScratch_, std::make_pair(maxBS, maxBS + bs), Kokkos::ALL());
  for(size_t j = 0; j < numVecs; j++)
    for(size_t i = 0; i < numRows; i++)
      X_scaled(i, j) = X_localBlocks_[blockIndex](i, j) * D_local(i, 0);

  HostSubviewLocal Y_temp(weightedApplyScratch_, std::make_pair(maxBS * 2, maxBS * 2 + bs), Kokkos::ALL());
  // Apply the local operator: Y_temp := M^{-1} * X_scaled
  this->solveBlock (X_scaled, Y_temp, blockIndex, mode, STS::one(), STS::zero());
  // Y_local := beta * Y_local + alpha * diag(D_local) * Y_temp.
  //
  // Note that we still use the permuted subset scaling D_local here,
  // because Y_temp has the same permuted subset Map.  That's good, in
  // fact, because it's a subset: less data to read and multiply.
  LISC a = alpha;
  LISC b = beta;
  for(size_t j = 0; j < numVecs; j++)
    for(size_t i = 0; i < numRows; i++)
      Y_localBlocks_[blockIndex](i, j) = b * Y_localBlocks_[blockIndex](i, j) + a * Y_temp(i, j) * D_local(i, 0);

  // Copy the permuted subset output vector Y_local into the original
  // output multivector Y.
  mvgs.scatterViewToView (Y, Y_localBlocks_[blockIndex], blockRows);
}

template<class MatrixType, class LocalScalarType>
Details::StridedRowView<
  typename ContainerImpl<MatrixType, LocalScalarType>::SC,
  typename ContainerImpl<MatrixType, LocalScalarType>::LO,
  typename ContainerImpl<MatrixType, LocalScalarType>::GO,
  typename ContainerImpl<MatrixType, LocalScalarType>::NO>
ContainerImpl<MatrixType, LocalScalarType>::
getInputRowView(LO row) const
{  

  typedef typename MatrixType::nonconst_local_inds_host_view_type nonconst_local_inds_host_view_type;
  typedef typename MatrixType::nonconst_values_host_view_type nonconst_values_host_view_type;

  typedef typename MatrixType::local_inds_host_view_type local_inds_host_view_type;
  typedef typename MatrixType::values_host_view_type values_host_view_type;
  using IST = typename row_matrix_type::impl_scalar_type;

  if(this->hasBlockCrs_)
  {
    using h_inds_type = typename block_crs_matrix_type::local_inds_host_view_type;
    using h_vals_type = typename block_crs_matrix_type::values_host_view_type;
    h_inds_type colinds;
    h_vals_type values;
    this->inputBlockMatrix_->getLocalRowView(row / this->bcrsBlockSize_, colinds, values);
    size_t numEntries = colinds.size();
    // CMS: Can't say I understand what this really does
    //return StridedRowView(values + row % this->bcrsBlockSize_, colinds, this->bcrsBlockSize_, numEntries * this->bcrsBlockSize_);
    h_vals_type subvals = Kokkos::subview(values,std::pair<size_t,size_t>(row % this->bcrsBlockSize_,values.size()));
    return StridedRowView(subvals, colinds, this->bcrsBlockSize_, numEntries * this->bcrsBlockSize_);
  }
  else if(!this->inputMatrix_->supportsRowViews())
  {
    size_t maxEntries = this->inputMatrix_->getLocalMaxNumRowEntries();
    Teuchos::Array<LO> inds(maxEntries);
    Teuchos::Array<SC> vals(maxEntries);
    nonconst_local_inds_host_view_type inds_v(inds.data(),maxEntries);
    nonconst_values_host_view_type vals_v(reinterpret_cast<IST*>(vals.data()),maxEntries);
    size_t numEntries;
    this->inputMatrix_->getLocalRowCopy(row, inds_v, vals_v, numEntries);
    vals.resize(numEntries); inds.resize(numEntries);
    return StridedRowView(vals, inds);
  }
  else
  {
    // CMS - This is dangerous and might not work.
    local_inds_host_view_type colinds;
    values_host_view_type values;
    this->inputMatrix_->getLocalRowView(row, colinds, values);
    return StridedRowView(values, colinds, 1, colinds.size());
  }
}

template<class MatrixType, class LocalScalarType>
void ContainerImpl<MatrixType, LocalScalarType>::
clearBlocks()
{
  X_localBlocks_.clear();
  Y_localBlocks_.clear();
  X_local_ = HostViewLocal();
  Y_local_ = HostViewLocal();
  Container<MatrixType>::clearBlocks();
}

namespace Details {

//Implementation of Ifpack2::Details::StridedRowView
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StridedRowView<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
StridedRowView(h_vals_type vals_, h_inds_type inds_, int blockSize_, size_t nnz_)
  : vals(vals_), inds(inds_), blockSize(blockSize_), nnz(nnz_)
{}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
StridedRowView<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
StridedRowView(Teuchos::Array<SC>& vals_, Teuchos::Array<LO>& inds_)
  : vals(), inds(), blockSize(1), nnz(vals_.size())
{
  valsCopy.swap(vals_);
  indsCopy.swap(inds_);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar StridedRowView<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
val(size_t i) const
{
  #ifdef HAVE_IFPACK2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(i >= nnz, std::runtime_error,
      "Out-of-bounds access into Ifpack2::Container::StridedRowView");
  #endif
  if(vals.size() > 0)
  {
    if(blockSize == 1)
      return vals[i];
    else
      return vals[i * blockSize];
  }
  else
    return valsCopy[i];
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
LocalOrdinal StridedRowView<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
ind(size_t i) const
{
  #ifdef HAVE_IFPACK2_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(i >= nnz, std::runtime_error,
        "Out-of-bounds access into Ifpack2::Container::StridedRowView");
  #endif
  //inds is smaller than vals by a factor of the block size (dofs/node)
    if(inds.size() > 0)
  {
    if(blockSize == 1)
      return inds[i];
    else
      return inds[i / blockSize] * blockSize + i % blockSize;
  }
  else
    return indsCopy[i];
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
size_t StridedRowView<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
size() const
{
  return nnz;
}
}

}

template <class MatrixType>
std::ostream& operator<<(std::ostream& os, const Ifpack2::Container<MatrixType>& obj)
{
  return obj.print(os);
}

#define IFPACK2_CONTAINER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::Container<Tpetra::RowMatrix<S, LO, GO, N>>; \
  template class Ifpack2::ContainerImpl<Tpetra::RowMatrix<S, LO, GO, N>, S>; \
  template class Ifpack2::Details::StridedRowView<S, LO, GO, N>; \
  template std::ostream& operator<< <Tpetra::RowMatrix<S, LO, GO, N>>( \
      std::ostream& os, const Ifpack2::Container<Tpetra::RowMatrix<S, LO, GO, N>>& obj);

#endif

