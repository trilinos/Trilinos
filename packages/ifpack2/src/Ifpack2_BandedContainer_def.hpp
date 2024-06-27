// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_BANDEDCONTAINER_DEF_HPP
#define IFPACK2_BANDEDCONTAINER_DEF_HPP

#include "Teuchos_LAPACK.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include <iostream>
#include <sstream>

#ifdef HAVE_MPI
#  include <mpi.h>
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI

namespace Ifpack2 {

template<class MatrixType, class LocalScalarType>
BandedContainer<MatrixType, LocalScalarType>::
BandedContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                 const Teuchos::Array<Teuchos::Array<LO> >& partitions,
                 const Teuchos::RCP<const import_type>&,
                 bool pointIndexed) :
  ContainerImpl<MatrixType, LocalScalarType>(matrix, partitions, pointIndexed),
  ipiv_(this->blockRows_.size() * this->scalarsPerRow_),
  kl_(this->numBlocks_, -1),
  ku_(this->numBlocks_, -1),
  scalarOffsets_(this->numBlocks_)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! matrix->hasColMap (), std::invalid_argument, "Ifpack2::BandedContainer: "
    "The constructor's input matrix must have a column Map.");
}

template<class MatrixType, class LocalScalarType>
BandedContainer<MatrixType, LocalScalarType>::
~BandedContainer () {}

template<class MatrixType, class LocalScalarType>
void BandedContainer<MatrixType, LocalScalarType>::
setParameters (const Teuchos::ParameterList& List)
{
  if(List.isParameter("relaxation: banded container superdiagonals"))
  {
    int ku = List.get<int>("relaxation: banded container superdiagonals");
    for(LO b = 0; b < this->numBlocks_; b++)
      ku_[b] = ku;
  }
  if(List.isParameter("relaxation: banded container subdiagonals"))
  {
    int kl = List.get<int>("relaxation: banded container subdiagonals");
    for(LO b = 0; b < this->numBlocks_; b++)
      kl_[b] = kl;
  }
}

template<class MatrixType, class LocalScalarType>
void BandedContainer<MatrixType, LocalScalarType>::
computeBandwidth()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  //now, for any block where kl_ or ku_ has not already been set, compute the actual bandwidth
  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  size_t colToOffsetSize = this->inputMatrix_->getLocalNumCols();
  if(this->pointIndexed_)
    colToOffsetSize *= this->bcrsBlockSize_;
  Array<LO> colToBlockOffset(colToOffsetSize, INVALID);
  //Same logic as extract() to find entries in blocks efficiently
  //(but here, just to find bandwidth, no scalars used)
  for(int i = 0; i < this->numBlocks_; i++)
  {
    //maxSub, maxSuper are the maximum lower and upper bandwidth
    LO maxSub = 0;
    LO maxSuper = 0;
    if(this->scalarsPerRow_ > 1)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockEnd = (i == this->numBlocks_ - 1) ? this->blockRows_.size() : this->blockOffsets_[i + 1];
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(size_t j = 0; j < size_t(blockRows.size()); j++)
      {
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }

      using h_inds_type = typename block_crs_matrix_type::local_inds_host_view_type;
      using h_vals_type = typename block_crs_matrix_type::values_host_view_type;
      for(LO blockRow = 0; blockRow < LO(blockRows.size()); blockRow++)
      {
        //get a raw view of the whole block row
        h_inds_type indices;
        h_vals_type values;
        LO inputRow = this->blockRows_[blockStart + blockRow];
        this->inputBlockMatrix_->getLocalRowView(inputRow, indices, values);
        LO numEntries = (LO) indices.size();
        for(LO k = 0; k < numEntries; k++)
        {
          LO colOffset = colToBlockOffset[indices[k]];
          if(blockStart <= colOffset && colOffset < blockEnd)
          {
            //This entry does appear in the diagonal block.
            //(br, bc) identifies the scalar's position in the BlockCrs block.
            //Convert this to (r, c) which is its position in the container block.
            LO blockCol = colOffset - blockStart;
            for(LO bc = 0; bc < this->bcrsBlockSize_; bc++)
            {
              for(LO br = 0; br < this->bcrsBlockSize_; br++)
              {
                LO r = this->bcrsBlockSize_ * blockRow + br;
                LO c = this->bcrsBlockSize_ * blockCol + bc;
                if(r - c > maxSub)
                  maxSub = r - c;
                if(c - r > maxSuper)
                  maxSuper = c - r;
              }
            }
          }
        }
      }
    }
    else
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockEnd = (i == this->numBlocks_ - 1) ? this->blockRows_.size() : this->blockOffsets_[i + 1];
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(size_t j = 0; j < size_t(blockRows.size()); j++)
      {
        //translateRowToCol will return the corresponding split column
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      for(LO blockRow = 0; blockRow < LO(blockRows.size()); blockRow++)
      {
        //get a view of the general row
        LO inputSplitRow = this->blockRows_[blockStart + blockRow];
        auto rowView = this->getInputRowView(inputSplitRow);
        for(size_t k = 0; k < rowView.size(); k++)
        {
          LO colOffset = colToBlockOffset[rowView.ind(k)];
          if(blockStart <= colOffset && colOffset < blockEnd)
          {
            LO blockCol = colOffset - blockStart;
            maxSub = std::max(maxSub, blockRow - blockCol);
            maxSuper = std::max(maxSuper, blockCol - blockRow);
          }
        }
      }
    }
    kl_[i] = maxSub;
    ku_[i] = maxSuper;
  }
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType>::
initialize ()
{
  //note: in BlockRelaxation, Container_->setParameters() immediately
  //precedes Container_->initialize(). So no matter what, either all
  //the block bandwidths were set (in setParameters()) or none of them were.
  //If none were they must be computed individually.
  if(kl_[0] == -1)
    computeBandwidth();
  GO totalScalars = 0;
  for(LO b = 0; b < this->numBlocks_; b++)
  {
    LO stride = 2 * kl_[b] + ku_[b] + 1;
    scalarOffsets_[b] = totalScalars;
    totalScalars += stride * this->blockSizes_[b] * this->scalarsPerRow_;
  }
  scalars_.resize(totalScalars);
  for(int b = 0; b < this->numBlocks_; b++)
  {
    //NOTE: the stride and upper bandwidth used to construct the SerialBandDenseMatrix looks
    //too large, but the extra kl_ in upper band space is needed by the LAPACK factorization routine.
    LO nrows = this->blockSizes_[b] * this->scalarsPerRow_;
    diagBlocks_.emplace_back(Teuchos::View, scalars_.data() + scalarOffsets_[b], 2 * kl_[b] + ku_[b] + 1, nrows, nrows, kl_[b], kl_[b] + ku_[b]);
    diagBlocks_[b].putScalar(Teuchos::ScalarTraits<LSC>::zero());
  }
  std::fill (ipiv_.begin (), ipiv_.end (), 0);
  // We assume that if you called this method, you intend to recompute
  // everything.
  this->IsComputed_ = false;
  this->IsInitialized_ = true;
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType>::
compute ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ipiv_.size () != this->blockRows_.size() * this->scalarsPerRow_, std::logic_error,
    "Ifpack2::BandedContainer::compute: ipiv_ array has the wrong size.  "
    "Please report this bug to the Ifpack2 developers.");

  this->IsComputed_ = false;
  if (!this->isInitialized ()) {
    this->initialize ();
  }

  extract (); // Extract the submatrices
  factor ();  // Factor them

  this->IsComputed_ = true;
}

template<class MatrixType, class LocalScalarType>
void BandedContainer<MatrixType, LocalScalarType>::extract()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using STS = Teuchos::ScalarTraits<SC>;
  SC zero = STS::zero();
  const LO INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  //To extract diagonal blocks, need to translate local rows to local columns.
  //Strategy: make a lookup table that translates local cols in the matrix to offsets in blockRows_:
  //blockOffsets_[b] <= offset < blockOffsets_[b+1]: tests whether the column is in block b.
  //offset - blockOffsets_[b]: gives the column within block b.
  //
  //This provides the block and col within a block in O(1).
  if(this->scalarsPerRow_ > 1)
  {
    Array<LO> colToBlockOffset(this->inputBlockMatrix_->getLocalNumCols(), INVALID);
    for(int i = 0; i < this->numBlocks_; i++)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockEnd = blockStart + this->blockSizes_[i];
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(size_t j = 0; j < size_t(blockRows.size()); j++)
      {
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      using h_inds_type = typename block_crs_matrix_type::local_inds_host_view_type;
      using h_vals_type = typename block_crs_matrix_type::values_host_view_type;
      for(LO blockRow = 0; blockRow < LO(blockRows.size()); blockRow++)
      {
        //get a raw view of the whole block row
        h_inds_type indices;
        h_vals_type values;
        LO inputRow = this->blockRows_[blockStart + blockRow];
        this->inputBlockMatrix_->getLocalRowView(inputRow, indices, values);
        LO numEntries = (LO) indices.size();        
        for(LO k = 0; k < numEntries; k++)
        {
          LO colOffset = colToBlockOffset[indices[k]];
          if(blockStart <= colOffset && colOffset < blockEnd)
          {
            //This entry does appear in the diagonal block.
            //(br, bc) identifies the scalar's position in the BlockCrs block.
            //Convert this to (r, c) which is its position in the container block.
            LO blockCol = colOffset - blockStart;
            for(LO bc = 0; bc < this->bcrsBlockSize_; bc++)
            {
              for(LO br = 0; br < this->bcrsBlockSize_; br++)
              {
                LO r = this->bcrsBlockSize_ * blockRow + br;
                LO c = this->bcrsBlockSize_ * blockCol + bc;
                auto val = values[k * (this->bcrsBlockSize_ * this->bcrsBlockSize_) + (br + this->bcrsBlockSize_ * bc)];
                if(val != zero)
                  diagBlocks_[i](r, c) = val;
              }
            }
          }
        }
      }
    }
  }
  else
  {
    //get the mapping from point-indexed matrix columns to offsets in blockRows_
    //(this includes regular CrsMatrix columns, in which case bcrsBlockSize_ == 1)
    Array<LO> colToBlockOffset(this->inputMatrix_->getLocalNumCols() * this->bcrsBlockSize_, INVALID);
    for(int i = 0; i < this->numBlocks_; i++)
    {
      //Get the interval where block i is defined in blockRows_
      LO blockStart = this->blockOffsets_[i];
      LO blockEnd = blockStart + this->blockSizes_[i];
      ArrayView<const LO> blockRows = this->getBlockRows(i);
      //Set the lookup table entries for the columns appearing in block i.
      //If OverlapLevel_ > 0, then this may overwrite values for previous blocks, but
      //this is OK. The values updated here are only needed to process block i's entries.
      for(size_t j = 0; j < size_t(blockRows.size()); j++)
      {
        //translateRowToCol will return the corresponding split column
        LO localCol = this->translateRowToCol(blockRows[j]);
        colToBlockOffset[localCol] = blockStart + j;
      }
      for(size_t blockRow = 0; blockRow < size_t(blockRows.size()); blockRow++)
      {
        //get a view of the split row
        LO inputPointRow = this->blockRows_[blockStart + blockRow];
        auto rowView = this->getInputRowView(inputPointRow);
        for(size_t k = 0; k < rowView.size(); k++)
        {
          LO colOffset = colToBlockOffset[rowView.ind(k)];
          if(blockStart <= colOffset && colOffset < blockEnd)
          {
            LO blockCol = colOffset - blockStart;
            auto val = rowView.val(k);
            if(val != zero)
              diagBlocks_[i](blockRow, blockCol) = rowView.val(k);
          }
        }
      }
    }
  }
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType>::
clearBlocks ()
{
  diagBlocks_.clear();
  scalars_.clear();
  ContainerImpl<MatrixType, LocalScalarType>::clearBlocks();
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType>::
factor ()
{
  Teuchos::LAPACK<int, LSC> lapack;
  int INFO = 0;

  for(int i = 0; i < this->numBlocks_; i++)
  {
    // Plausibility checks for matrix
    TEUCHOS_TEST_FOR_EXCEPTION(diagBlocks_[i].values()==0, std::invalid_argument,
                       "BandedContainer<T>::factor: Diagonal block is an empty SerialBandDenseMatrix<T>!");
    TEUCHOS_TEST_FOR_EXCEPTION(diagBlocks_[i].upperBandwidth() < diagBlocks_[i].lowerBandwidth(), std::invalid_argument,
                       "BandedContainer<T>::factor: Diagonal block needs kl additional superdiagonals for factorization! However, the number of superdiagonals is smaller than the number of subdiagonals!");
    int* blockIpiv = &ipiv_[this->blockOffsets_[i] * this->scalarsPerRow_];
    lapack.GBTRF (diagBlocks_[i].numRows(),
        diagBlocks_[i].numCols(),
        diagBlocks_[i].lowerBandwidth(),
        diagBlocks_[i].upperBandwidth() - diagBlocks_[i].lowerBandwidth(), /* enter the real number of superdiagonals (see Teuchos_SerialBandDenseSolver)*/
        diagBlocks_[i].values(),
        diagBlocks_[i].stride(),
        blockIpiv,
        &INFO);

    // INFO < 0 is a bug.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO < 0, std::logic_error, "Ifpack2::BandedContainer::factor: "
      "LAPACK's _GBTRF (LU factorization with partial pivoting) was called "
      "incorrectly.  INFO = " << INFO << " < 0.  "
      "Please report this bug to the Ifpack2 developers.");
    // INFO > 0 means the matrix is singular.  This is probably an issue
    // either with the choice of rows the rows we extracted, or with the
    // input matrix itself.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO > 0, std::runtime_error, "Ifpack2::BandedContainer::factor: "
      "LAPACK's _GBTRF (LU factorization with partial pivoting) reports that the "
      "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
      "(one-based index i) is exactly zero.  This probably means that the input "
      "matrix has a singular diagonal block.");
  }
}

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType>::
solveBlock(ConstHostSubviewLocal X,
           HostSubviewLocal Y,
           int blockIndex,
           Teuchos::ETransp mode,
           const LSC alpha,
           const LSC beta) const
{
  #ifdef HAVE_IFPACK2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != Y.extent (0),
    std::logic_error, "Ifpack2::BandedContainer::solveBlock: X and Y have "
    "incompatible dimensions (" << X.extent (0) << " resp. "
    << Y.extent (0) << ").  Please report this bug to "
    "the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != static_cast<size_t> (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numCols() : diagBlocks_[blockIndex].numRows()),
    std::logic_error, "Ifpack2::BandedContainer::solveBlock: The input "
    "multivector X has incompatible dimensions from those of the "
    "inverse operator (" << X.extent (0) << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numCols() : diagBlocks_[blockIndex].numRows())
    << ").  Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y.extent (0) != static_cast<size_t> (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numRows() : diagBlocks_[blockIndex].numCols()),
    std::logic_error, "Ifpack2::BandedContainer::solveBlock: The output "
    "multivector Y has incompatible dimensions from those of the "
    "inverse operator (" << Y.extent (0) << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numRows() : diagBlocks_[blockIndex].numCols())
    << ").  Please report this bug to the Ifpack2 developers.");
  #endif

  size_t numRows = X.extent (0);
  size_t numVecs = X.extent (1);

  LSC zero = Teuchos::ScalarTraits<LSC>::zero ();
  if (alpha == zero) { // don't need to solve the linear system
    if (beta == zero) {
      // Use BLAS AXPY semantics for beta == 0: overwrite, clobbering
      // any Inf or NaN values in Y (rather than multiplying them by
      // zero, resulting in NaN values).
      for(size_t j = 0; j < numVecs; j++)
        for(size_t i = 0; i < numRows; i++)
          Y(i, j) = zero;
    }
    else { // beta != 0
      for(size_t j = 0; j < numVecs; j++)
        for(size_t i = 0; i < numRows; i++)
          Y(i, j) = beta * Y(i, j);
    }
  }
  else { // alpha != 0; must solve the linear system
    Teuchos::LAPACK<int, LSC> lapack;
    // If beta is nonzero or Y is not constant stride, we have to use
    // a temporary output multivector.  It gets a copy of X, since
    // GBTRS overwrites its (multi)vector input with its output.
    std::vector<LSC> yTemp(numVecs * numRows);
    for(size_t j = 0; j < numVecs; j++)
    {
      for(size_t i = 0; i < numRows; i++)
        yTemp[j * numRows + i] = X(i, j);
    }

    int INFO = 0;
    const char trans = (mode == Teuchos::CONJ_TRANS ? 'C' : (mode == Teuchos::TRANS ? 'T' : 'N'));

    const int* blockIpiv = &ipiv_[this->blockOffsets_[blockIndex] * this->scalarsPerRow_];

    lapack.GBTRS(trans,
        diagBlocks_[blockIndex].numCols(),
        diagBlocks_[blockIndex].lowerBandwidth(),
        /* enter the real number of superdiagonals (see Teuchos_SerialBandDenseSolver)*/
        diagBlocks_[blockIndex].upperBandwidth() - diagBlocks_[blockIndex].lowerBandwidth(),
        numVecs,
        diagBlocks_[blockIndex].values(),
        diagBlocks_[blockIndex].stride(),
        blockIpiv,
        yTemp.data(),
        numRows,
        &INFO);

    if (beta != zero) {
      for(size_t j = 0; j < numVecs; j++)
      {
        for(size_t i = 0; i < numRows; i++)
        {
          Y(i, j) *= ISC(beta);
          Y(i, j) += ISC(alpha * yTemp[j * numRows + i]);
        }
      }
    }
    else {
      for(size_t j = 0; j < numVecs; j++)
      {
        for(size_t i = 0; i < numRows; i++)
          Y(i, j) = ISC(yTemp[j * numRows + i]);
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO != 0, std::runtime_error, "Ifpack2::BandedContainer::solveBlock: "
      "LAPACK's _GBTRS (solve using LU factorization with partial pivoting) "
      "failed with INFO = " << INFO << " != 0.");
  }
}

template<class MatrixType, class LocalScalarType>
std::ostream&
BandedContainer<MatrixType, LocalScalarType>::
print (std::ostream& os) const
{
  Teuchos::FancyOStream fos (Teuchos::rcpFromRef (os));
  fos.setOutputToRootOnly (0);
  describe (fos);
  return os;
}

template<class MatrixType, class LocalScalarType>
std::string
BandedContainer<MatrixType, LocalScalarType>::
description () const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (this->isInitialized()) {
    if (this->isComputed()) {
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

template<class MatrixType, class LocalScalarType>
void
BandedContainer<MatrixType, LocalScalarType>::
describe (Teuchos::FancyOStream& os,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << std::endl;
  os << "Ifpack2::BandedContainer" << std::endl;
  for(int i = 0; i < this->numBlocks_; i++)
  {
    os << "Block " << i << ": Number of rows           = " << this->blockSizes_[i] << std::endl;
    os << "Block " << i << ": Number of subdiagonals   = " << diagBlocks_[i].lowerBandwidth() << std::endl;
    os << "Block " << i << ": Number of superdiagonals = " << diagBlocks_[i].upperBandwidth() << std::endl;
  }
  os << "isInitialized()          = " << this->IsInitialized_ << std::endl;
  os << "isComputed()             = " << this->IsComputed_ << std::endl;
  os << "================================================================================" << std::endl;
  os << std::endl;
}

template<class MatrixType, class LocalScalarType>
std::string BandedContainer<MatrixType, LocalScalarType>::getName()
{
  return "Banded";
}

} // namespace Ifpack2

#define IFPACK2_BANDEDCONTAINER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::BandedContainer< Tpetra::RowMatrix<S, LO, GO, N>, S >;

#endif // IFPACK2_BANDEDCONTAINER_HPP
