// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DENSECONTAINER_DEF_HPP
#define IFPACK2_DENSECONTAINER_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Tpetra_BlockMultiVector.hpp"

#ifdef HAVE_MPI
#  include <mpi.h>
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI

namespace Ifpack2 {

template<class MatrixType, class LocalScalarType>
DenseContainer<MatrixType, LocalScalarType>::
DenseContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                const Teuchos::Array<Teuchos::Array<LO> >& partitions,
                const Teuchos::RCP<const import_type>&,
                bool pointIndexed) :
  ContainerImpl<MatrixType, LocalScalarType> (matrix, partitions, pointIndexed),
  scalarOffsets_ (this->numBlocks_)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    !matrix->hasColMap(), std::invalid_argument, "Ifpack2::DenseContainer: "
    "The constructor's input matrix must have a column Map.");

  //compute scalarOffsets_
  GO totalScalars = 0;
  for(LO i = 0; i < this->numBlocks_; i++)
  {
    scalarOffsets_[i] = totalScalars;
    totalScalars += this->blockSizes_[i] * this->scalarsPerRow_ *
                    this->blockSizes_[i] * this->scalarsPerRow_;
  }
  scalars_.resize(totalScalars);
  for(int i = 0; i < this->numBlocks_; i++)
  {
    LO denseRows = this->blockSizes_[i] * this->scalarsPerRow_;
    //create square dense matrix (stride is same as rows and cols)
    diagBlocks_.emplace_back(Teuchos::View, scalars_.data() + scalarOffsets_[i], denseRows, denseRows, denseRows);
  }

  ipiv_.resize(this->blockRows_.size() * this->scalarsPerRow_);
}

template<class MatrixType, class LocalScalarType>
DenseContainer<MatrixType, LocalScalarType>::
~DenseContainer () {}

template<class MatrixType, class LocalScalarType>
void
DenseContainer<MatrixType, LocalScalarType>::
initialize ()
{
  // Fill the diagonal block and LU permutation array with zeros.
  for(int i = 0; i < this->numBlocks_; i++)
    diagBlocks_[i].putScalar(Teuchos::ScalarTraits<LSC>::zero());
  std::fill (ipiv_.begin (), ipiv_.end (), 0);

  this->IsInitialized_ = true;
  // We assume that if you called this method, you intend to recompute
  // everything.
  this->IsComputed_ = false;
}

template<class MatrixType, class LocalScalarType>
void
DenseContainer<MatrixType, LocalScalarType>::
compute ()
{
// FIXME: I am commenting this out because it breaks block CRS support
//  TEUCHOS_TEST_FOR_EXCEPTION(
//    static_cast<size_t> (ipiv_.size ()) != numRows_, std::logic_error,
//    "Ifpack2::DenseContainer::compute: ipiv_ array has the wrong size.  "
//    "Please report this bug to the Ifpack2 developers.");

  this->IsComputed_ = false;
  if (!this->isInitialized ()) {
    this->initialize();
  }

  extract (); // Extract the submatrices
  factor ();  // Factor them
  this->IsComputed_ = true;
}

template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType, LocalScalarType>::extract()
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
DenseContainer<MatrixType, LocalScalarType>::
factor ()
{
  Teuchos::LAPACK<int, LSC> lapack;
  for(int i = 0; i < this->numBlocks_; i++)
  {
    int INFO = 0;
    int* blockIpiv = &ipiv_[this->blockOffsets_[i] * this->scalarsPerRow_];
    lapack.GETRF(diagBlocks_[i].numRows(),
                 diagBlocks_[i].numCols(),
                 diagBlocks_[i].values(),
                 diagBlocks_[i].stride(),
                 blockIpiv, &INFO);
    // INFO < 0 is a bug.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO < 0, std::logic_error, "Ifpack2::DenseContainer::factor: "
      "LAPACK's _GETRF (LU factorization with partial pivoting) was called "
      "incorrectly.  INFO = " << INFO << " < 0.  "
      "Please report this bug to the Ifpack2 developers.");
    // INFO > 0 means the matrix is singular.  This is probably an issue
    // either with the choice of rows the rows we extracted, or with the
    // input matrix itself.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO > 0, std::runtime_error, "Ifpack2::DenseContainer::factor: "
      "LAPACK's _GETRF (LU factorization with partial pivoting) reports that the "
      "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
      "(one-based index i) is exactly zero.  This probably means that the input "
      "matrix has a singular diagonal block.");
  }
}

template<class MatrixType, class LocalScalarType>
void
DenseContainer<MatrixType, LocalScalarType>::
solveBlock(ConstHostSubviewLocal X,
           HostSubviewLocal Y,
           int blockIndex,
           Teuchos::ETransp mode,
           LSC alpha,
           LSC beta) const
{
  #ifdef HAVE_IFPACK2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != Y.extent (0),
    std::logic_error, "Ifpack2::DenseContainer::solveBlock: X and Y have "
    "incompatible dimensions (" << X.extent (0) << " resp. "
    << Y.extent (0) << ").  Please report this bug to "
    "the Ifpack2 developers.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (1) != Y.extent(1),
    std::logic_error, "Ifpack2::DenseContainer::solveBlock: X and Y have "
    "incompatible numbers of vectors (" << X.extent (1) << " resp. "
    << Y.extent (1) << ").  Please report this bug to "
    "the Ifpack2 developers.");
  #endif

  typedef Teuchos::ScalarTraits<LSC> STS;
  size_t numRows = X.extent(0);
  size_t numVecs = X.extent(1);
  if(alpha == STS::zero()) { // don't need to solve the linear system
    if(beta == STS::zero()) {
      // Use BLAS AXPY semantics for beta == 0: overwrite, clobbering
      // any Inf or NaN values in Y (rather than multiplying them by
      // zero, resulting in NaN values).
      for(size_t j = 0; j < numVecs; j++)
      {
        for(size_t i = 0; i < numRows; i++)
          Y(i, j) = STS::zero();
      }
    }
    else // beta != 0
      for(size_t j = 0; j < numVecs; j++)
      {
        for(size_t i = 0; i < numRows; i++)
          Y(i, j) = beta * Y(i, j);
      }
  }
  else { // alpha != 0; must solve the linear system
    Teuchos::LAPACK<int, LSC> lapack;
    // If beta is nonzero or Y is not constant stride, we have to use
    // a temporary output multivector.  It gets a (deep) copy of X,
    // since GETRS overwrites its (multi)vector input with its output.
    std::vector<LSC> yTemp(numVecs * numRows);
    for(size_t j = 0; j < numVecs; j++)
    {
      for(size_t i = 0; i < numRows; i++)
        yTemp[j * numRows + i] = X(i, j);
    }

    int INFO = 0;
    int* blockIpiv = &ipiv_[this->blockOffsets_[blockIndex] * this->scalarsPerRow_];
    const char trans =
      (mode == Teuchos::CONJ_TRANS ? 'C' : (mode == Teuchos::TRANS ? 'T' : 'N'));
    lapack.GETRS (trans,
                  diagBlocks_[blockIndex].numRows (),
                  numVecs,
                  diagBlocks_[blockIndex].values (),
                  diagBlocks_[blockIndex].stride (),
                  blockIpiv,
                  yTemp.data(),
                  numRows,
                  &INFO);

    if (beta != STS::zero ()) {
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
      INFO != 0, std::runtime_error, "Ifpack2::DenseContainer::solveBlock: "
      "LAPACK's _GETRS (solve using LU factorization with partial pivoting) "
      "failed with INFO = " << INFO << " != 0.");
  }
}

template<class MatrixType, class LocalScalarType>
std::ostream&
DenseContainer<MatrixType, LocalScalarType>::
print (std::ostream& os) const
{
  Teuchos::FancyOStream fos (Teuchos::rcpFromRef (os));
  fos.setOutputToRootOnly (0);
  this->describe (fos);
  return os;
}

template<class MatrixType, class LocalScalarType>
std::string
DenseContainer<MatrixType, LocalScalarType>::
description () const
{
  std::ostringstream oss;
  oss << "Ifpack::DenseContainer: ";
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
DenseContainer<MatrixType, LocalScalarType>::
describe (Teuchos::FancyOStream& os,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << endl;
  os << "Ifpack2::DenseContainer" << endl;
  for(int i = 0; i < this->numBlocks_; i++)
  {
    os << "Block " << i << " number of rows          = " << this->blockSizes_[i] << endl;
  }
  os << "isInitialized()         = " << this->IsInitialized_ << endl;
  os << "isComputed()            = " << this->IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

template<class MatrixType, class LocalScalarType>
void DenseContainer<MatrixType, LocalScalarType>::clearBlocks()
{
  diagBlocks_.clear();
  scalars_.clear();
  ContainerImpl<MatrixType, LocalScalarType>::clearBlocks();
}

template<class MatrixType, class LocalScalarType>
std::string DenseContainer<MatrixType, LocalScalarType>::getName()
{
  return "Dense";
}

} // namespace Ifpack2

// There's no need to instantiate for CrsMatrix too.  All Ifpack2
// preconditioners can and should do dynamic casts if they need a type
// more specific than RowMatrix.

#define IFPACK2_DENSECONTAINER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::DenseContainer< Tpetra::RowMatrix<S, LO, GO, N>, S >;

#endif // IFPACK2_DENSECONTAINER_DEF_HPP
