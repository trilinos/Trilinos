// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_TRIDICONTAINER_DEF_HPP
#define IFPACK2_TRIDICONTAINER_DEF_HPP

#include "Teuchos_LAPACK.hpp"

#ifdef HAVE_MPI
#  include <mpi.h>
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI


namespace Ifpack2 {

template<class MatrixType, class LocalScalarType>
TriDiContainer<MatrixType, LocalScalarType>::
TriDiContainer (const Teuchos::RCP<const row_matrix_type>& matrix,
                const Teuchos::Array<Teuchos::Array<LO> >& partitions,
                const Teuchos::RCP<const import_type>&,
                bool pointIndexed) :
  ContainerImpl<MatrixType, LocalScalarType>(matrix, partitions, pointIndexed),
  ipiv_(this->blockRows_.size() * this->scalarsPerRow_),
  scalarOffsets_ (this->numBlocks_)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! matrix->hasColMap (), std::invalid_argument, "Ifpack2::TriDiContainer: "
    "The constructor's input matrix must have a column Map.");
  // FIXME (mfh 25 Aug 2013) What if the matrix's row Map has a
  // different index base than zero?
  //compute scalar array offsets (probably different from blockOffsets_)
  LO scalarTotal = 0;
  for(LO i = 0; i < this->numBlocks_; i++)
  {
    scalarOffsets_[i] = scalarTotal;
    //actualBlockRows is how many scalars it takes to store the diagonal for block i
    LO actualBlockRows = this->scalarsPerRow_ * this->blockSizes_[i];
    if(actualBlockRows == 1)
    {
      scalarTotal += 1;
    }
    else
    {
      //this formula is exact for any block size of 1 or more
      //includes 1 subdiagonal and 2 superdiagonals.
      scalarTotal += 4 * (actualBlockRows  - 1);
    }
  }
  //Allocate scalar arrays
  scalars_.resize(scalarTotal);
  diagBlocks_.reserve(this->numBlocks_);
  for(int i = 0; i < this->numBlocks_; i++)
  {
    diagBlocks_.emplace_back(Teuchos::View, scalars_.data() + scalarOffsets_[i], this->blockSizes_[i] * this->scalarsPerRow_);
  }
}

template<class MatrixType, class LocalScalarType>
TriDiContainer<MatrixType, LocalScalarType>::
~TriDiContainer () {}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType>::initialize ()
{
  for(int i = 0; i < this->numBlocks_; i++)
    diagBlocks_[i].putScalar(Teuchos::ScalarTraits<LSC>::zero());
  this->IsInitialized_ = true;
  // We assume that if you called this method, you intend to recompute
  // everything.
  this->IsComputed_ = false;
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType>::compute ()
{
  #ifdef HAVE_IFPACK2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    ipiv_.size () != this->blockRows_.size() * this->scalarsPerRow_, std::logic_error,
    "Ifpack2::TriDiContainer::compute: ipiv_ array has the wrong size.  "
    "Please report this bug to the Ifpack2 developers.");
  #endif

  this->IsComputed_ = false;
  if (! this->isInitialized ()) {
    this->initialize ();
  }
  extract (); // extract the submatrices
  factor (); // factor them
  this->IsComputed_ = true;
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType>::extract()
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  SC zero = Teuchos::ScalarTraits<SC>::zero();
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
void TriDiContainer<MatrixType, LocalScalarType>::clearBlocks ()
{
  diagBlocks_.clear();
  scalars_.clear();
  ContainerImpl<MatrixType, LocalScalarType>::clearBlocks();
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType>::factor ()
{
  for(int i = 0; i < this->numBlocks_; i++)
  {
    Teuchos::LAPACK<int, LSC> lapack;
    int INFO = 0;
    int* blockIpiv = &ipiv_[this->blockOffsets_[i] * this->scalarsPerRow_];
    lapack.GTTRF (diagBlocks_[i].numRowsCols (),
                  diagBlocks_[i].DL(),
                  diagBlocks_[i].D(),
                  diagBlocks_[i].DU(),
                  diagBlocks_[i].DU2(),
                  blockIpiv, &INFO);
    // INFO < 0 is a bug.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO < 0, std::logic_error, "Ifpack2::TriDiContainer::factor: "
      "LAPACK's _GTTRF (LU factorization with partial pivoting) was called "
      "incorrectly.  INFO = " << INFO << " < 0.  "
      "Please report this bug to the Ifpack2 developers.");
    // INFO > 0 means the matrix is singular.  This is probably an issue
    // either with the choice of rows the rows we extracted, or with the
    // input matrix itself.
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO > 0, std::runtime_error, "Ifpack2::TriDiContainer::factor: "
      "LAPACK's _GTTRF (LU factorization with partial pivoting) reports that the "
      "computed U factor is exactly singular.  U(" << INFO << "," << INFO << ") "
      "(one-based index i) is exactly zero.  This probably means that the input "
      "matrix has a singular diagonal block.");
  }
}

template<class MatrixType, class LocalScalarType>
void TriDiContainer<MatrixType, LocalScalarType>::
solveBlock(ConstHostSubviewLocal X,
           HostSubviewLocal Y,
           int blockIndex,
           Teuchos::ETransp mode,
           LSC alpha,
           LSC beta) const
{
  LSC zero = Teuchos::ScalarTraits<LSC>::zero();
  size_t numVecs = X.extent(1);
  size_t numRows = X.extent(0);

  #ifdef HAVE_IFPACK2_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != Y.extent (0),
    std::logic_error, "Ifpack2::TriDiContainer::solveBlock: X and Y have "
    "incompatible dimensions (" << X.extent (0) << " resp. "
    << Y.extent (0) << ").  Please report this bug to "
    "the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.extent (0) != static_cast<size_t> (diagBlocks_[blockIndex].numRowsCols()),
    std::logic_error, "Ifpack2::TriDiContainer::solveBlock: The input "
    "multivector X has incompatible dimensions from those of the "
    "inverse operator (" << X.extent (0) << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numRowsCols () : diagBlocks_[blockIndex].numRowsCols())
    << ").  Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    Y.extent (0) != static_cast<size_t> (diagBlocks_[blockIndex].numRowsCols()),
    std::logic_error, "Ifpack2::TriDiContainer::solveBlock: The output "
    "multivector Y has incompatible dimensions from those of the "
    "inverse operator (" << Y.extent (0) << " vs. "
    << (mode == Teuchos::NO_TRANS ? diagBlocks_[blockIndex].numRowsCols() : diagBlocks_[blockIndex].numRowsCols ())
    << ").  Please report this bug to the Ifpack2 developers.");
  #endif

  if(alpha == zero) { // don't need to solve the linear system
    if(beta == zero) {
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
          Y(i, j) *= ISC(beta);
    }
  }
  else { // alpha != 0; must solve the linear system
    Teuchos::LAPACK<int, LSC> lapack;
    // If beta is nonzero or Y is not constant stride, we have to use
    // a temporary output multivector.  It gets a copy of X, since
    // GETRS overwrites its (multi)vector input with its output.
    
    std::vector<LSC> yTemp(numVecs * numRows);
    for(size_t j = 0; j < numVecs; j++)
    {
      for(size_t i = 0; i < numRows; i++)
        yTemp[j * numRows + i] = X(i, j);
    }

    int INFO = 0;
    const char trans =
      (mode == Teuchos::CONJ_TRANS ? 'C' : (mode == Teuchos::TRANS ? 'T' : 'N'));
    int* blockIpiv = &ipiv_[this->blockOffsets_[blockIndex] * this->scalarsPerRow_];
    lapack.GTTRS (trans,
                  diagBlocks_[blockIndex].numRowsCols(),
                  numVecs,
                  diagBlocks_[blockIndex].DL(),
                  diagBlocks_[blockIndex].D(),
                  diagBlocks_[blockIndex].DU(),
                  diagBlocks_[blockIndex].DU2(),
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
      INFO != 0, std::runtime_error, "Ifpack2::TriDiContainer::solveBlock: "
      "LAPACK's _GETRS (solve using LU factorization with partial pivoting) "
      "failed with INFO = " << INFO << " != 0.");
  }
}

template<class MatrixType, class LocalScalarType>
std::ostream& TriDiContainer<MatrixType, LocalScalarType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}

template<class MatrixType, class LocalScalarType>
std::string TriDiContainer<MatrixType, LocalScalarType>::description() const
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
TriDiContainer<MatrixType, LocalScalarType>::
describe (Teuchos::FancyOStream& os,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;
  os << "================================================================================" << endl;
  os << "Ifpack2::TriDiContainer" << endl;
  os << "Number of blocks        = " << this->numBlocks_ << endl;
  os << "isInitialized()         = " << this->IsInitialized_ << endl;
  os << "isComputed()            = " << this->IsComputed_ << endl;
  os << "================================================================================" << endl;
  os << endl;
}

template<class MatrixType, class LocalScalarType>
std::string TriDiContainer<MatrixType, LocalScalarType>::getName()
{
  return "TriDi";
}

#define IFPACK2_TRIDICONTAINER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::TriDiContainer< Tpetra::RowMatrix<S, LO, GO, N>, S >;

} // namespace Ifpack2

#endif
