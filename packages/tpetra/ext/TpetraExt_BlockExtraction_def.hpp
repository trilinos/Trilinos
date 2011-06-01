//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef TPETRAEXT_BLOCKEXTRACTION_DEF_HPP
#define TPETRAEXT_BLOCKEXTRACTION_DEF_HPP

#include <algorithm>
#include <map>

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
Tpetra::Ext::extractBlockDiagonals(
          const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & matrix, 
          const Teuchos::ArrayView<const LocalOrdinal>            & first_points,
          Teuchos::ArrayRCP<Scalar>                               & out_diags,
          Teuchos::ArrayRCP<LocalOrdinal>                         & out_offsets)
{
  // block meta-data
  const int numBlocks = (int)first_points.size()-1;
  const size_t numRows = matrix.getNodeNumRows();
  TEST_FOR_EXCEPTION(numRows != 0 && numBlocks <= 0, std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): specified zero blocks for a matrix with more than zero rows.");
  TEST_FOR_EXCEPTION(numBlocks > 0 && first_points[0] != 0, std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): first point of first block must be zero.");
  TEST_FOR_EXCEPTION(matrix.isFillComplete() == false, std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): matrix must be fill-complete.");
  // INVARIANT: for all i=0,...,numBlock-1: first_points[i] - first_points[i-1] = block_size[i]
  int alloc_size      = 0,
      sum_block_sizes = 0;
  if (numBlocks) {
    out_offsets = arcp<LocalOrdinal>(numBlocks);
  }
  for (int b=0; b < numBlocks; ++b) 
  {
    const int block_size_b = (int)first_points[b+1] - (int)first_points[b];
    TEST_FOR_EXCEPTION(block_size_b < 0, std::runtime_error,
        "Tpetra::Ext::extractBlockDiagonals(): first points are not coherent.");
    sum_block_sizes += block_size_b;
    out_offsets[b]   = alloc_size;
    alloc_size      += block_size_b*block_size_b;
  }
  TEST_FOR_EXCEPTION( sum_block_sizes != (int)matrix.getNodeNumRows(), std::runtime_error, 
      "Tpetra::Ext::extractBlockDiagonals(): specified blocks are not compatible with specified matrix (blocks are too large or too small or the last offset was missing).");
  if (alloc_size) {
    out_diags = arcp<Scalar>(alloc_size);
    // must explicitly fill with zeros, because we will only insert non-zeros below
    std::fill( out_diags.begin(), out_diags.end(), ScalarTraits<Scalar>::zero() );
  }
  // extract blocks from matrix
  if (alloc_size) {
    const LocalOrdinal first_row = matrix.getRowMap()->getMinLocalIndex(),
                        last_row = matrix.getRowMap()->getMaxLocalIndex();
    ArrayView<const Scalar> rowvals;
    ArrayView<const LocalOrdinal> rowinds;
    // b is ready to be incremented to zero
    // block is invalid
    // subrow and block_size_b are prepared to trigger the while loop and properly initialize the others
    int b            = -1, // zero minus one
        subrow       =  0, 
        block_size_b =  subrow;
    typename ArrayRCP<Scalar>::iterator block = out_diags.end();
    // loop over all local rows
    for (LocalOrdinal lrow = first_row; lrow <= last_row; ++lrow)
    {
      // the while loop accounts for blocks of size zero
      while (subrow == block_size_b) 
      {
        // we busted the block, move to the next
        b += 1;
        block_size_b = first_points[b+1] - first_points[b];
        // an iterator to the beginning of this particular space will ensure bounds in a debug build
        // in a release build, it will simply be a pointer
        if (block_size_b) {
          block = out_diags.persistingView( out_offsets[b], block_size_b*block_size_b ).begin();
        }
        subrow = 0;
      }
      // extract the row, put the members of the block diagonal into the current block
      matrix.getLocalRowView(lrow, rowinds, rowvals);
      for (int k=0; k < (int)rowinds.size(); ++k) {
        const int subcol = rowinds[k] - first_points[b];
        if (subcol >= 0 && subcol < block_size_b) {
          block[subcol*block_size_b + subrow] += rowvals[k];
        }
      }
      ++subrow;
    }
    // this should have simultaneously finished matrix and the last block
    TEST_FOR_EXCEPTION( subrow != block_size_b, std::logic_error,
        "Tpetra::Ext::extractBlockDiagonals(): internal logic error. Please contact Tpetra team.");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
Tpetra::Ext::extractBlockDiagonals(
          const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & matrix, 
          const BlockMap<LocalOrdinal,GlobalOrdinal,Node>         & block_map,
          Teuchos::ArrayRCP<Scalar>                               & out_diags,
          Teuchos::ArrayRCP<LocalOrdinal>                         & out_offsets)
{
  ArrayRCP<const LocalOrdinal> block_firsts = block_map.getNodeFirstPointInBlocks();
  try {
    extractBlockDiagonals( matrix, block_firsts(), out_diags, out_offsets );
  }
  catch (std::exception &e) {
    TEST_FOR_EXCEPTION(true, std::runtime_error, 
        "Tpetra::Ext::extractBlockDiagonals(RowMatrix,BlockMap,...) caught exception:\n\n" << e.what());
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
Tpetra::Ext::extractBlockRow(LocalOrdinal localBlockRow,
                             const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &matrix, 
                             const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>         &block_row_map,
                             const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>         &block_col_map,
                             Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> >                   &out_block_entries,
                             Teuchos::ArrayRCP<LocalOrdinal>                                 &out_block_indices)
{
  std::string errstr("Tpetra::Ext::extractBlockRow(): ");
  TEST_FOR_EXCEPTION( matrix.isFillComplete() == false, std::runtime_error, 
      errstr << "matrix not fill-complete.");
  TEST_FOR_EXCEPTION( !( localBlockRow >= 0 && (size_t)localBlockRow < block_row_map.getNodeNumBlocks() ), std::runtime_error,
      errstr << "requested block row is not valid:\n" 
            << "localBlockRow == " << localBlockRow << "\n"
            << "block_row_map.getNodeNumBlocks() == " << block_row_map.getNodeNumBlocks() << "\n");
  TEST_FOR_EXCEPTION( ! block_row_map.getPointMap()->isCompatible( *matrix.getRowMap() ), std::runtime_error,
      errstr << "specified block row map is not compatible with the matix row map.");
  // 
  out_block_entries = null;
  out_block_indices = null;
  size_t numAllocatedBlocks = 0;
  //
  typedef typename ArrayView<const LocalOrdinal>::iterator AVLOIter;
  typedef typename ArrayView<const Scalar      >::iterator AVSIter;
  // 
  const size_t blockSize = block_row_map.getLocalBlockSize(localBlockRow);
  const ArrayRCP<const LocalOrdinal> firstPoints = block_col_map.getNodeFirstPointInBlocks();
  const LocalOrdinal rowOffset = block_row_map.getFirstLocalPointInLocalBlock( localBlockRow );
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> &colMap = *matrix.getColMap();
  std::map<LocalOrdinal,ArrayRCP<Scalar> > blocks;
  for (size_t subRow=0; subRow < blockSize; ++subRow) 
  {
    ArrayView<const LocalOrdinal> rowIndices;
    ArrayView<const Scalar>       rowValues;
    matrix.getLocalRowView( rowOffset+subRow, rowIndices, rowValues );
    const AVLOIter rowindend = rowIndices.end();
    AVLOIter       rowind    = rowIndices.begin();
    AVSIter        rowval    = rowValues.begin();
    for (; rowind != rowindend; ++rowind, ++rowval)
    {
      // what is the global id for this entry?
      const GlobalOrdinal mGind = colMap.getGlobalElement( *rowind );
      // what is the local id of that global id, under the block column map's point map?
      const LocalOrdinal  bLind = block_col_map.getPointMap()->getLocalElement( mGind );
      // did it even exist?
      if (bLind != OrdinalTraits<LocalOrdinal>::invalid()) {
        // which block does that element belong to, according to the partitioning in the block column map?
        typename ArrayRCP<const LocalOrdinal>::iterator nfp;
        nfp = std::upper_bound( firstPoints.begin(), firstPoints.end(), bLind );
#ifdef HAVE_TPETRA_DEBUG
        // firstPoints[0] == 0, so search for any bLind >= 0 should return nfp > firstPoints.end()
        TEST_FOR_EXCEPTION( nfp == firstPoints.begin(), std::logic_error, errstr << "internal Tpetra logic error. Please contact Tpetra team.\n");
#endif
        const size_t bj     = nfp - firstPoints.begin() - 1;
        // which sub-column in that block?
        const size_t subCol = bLind - *(nfp-1);
        const size_t blockWidth = block_col_map.getLocalBlockSize( bj );
#ifdef HAVE_TPETRA_DEBUG
        // if this fails, it is probably an error with blockWidth from BlockMap
        TEST_FOR_EXCEPTION( subCol >= blockWidth, std::logic_error, errstr << "internal Tpetra logic error. Please contact Tpetra team.\n");
#endif
        // find the block by looking it up in our std::map object
        ArrayRCP<Scalar> &block = blocks[bj];
        if (block == null) {
          // wasn't yet allocated; do so
          block = arcp<Scalar>(blockWidth*blockSize);
          std::fill(block.begin(), block.end(), ScalarTraits<Scalar>::zero());
          ++numAllocatedBlocks;
        }
        else {
#ifdef HAVE_TPETRA_DEBUG
          TEST_FOR_EXCEPTION((size_t)block.size() != blockWidth*blockSize, std::logic_error, errstr << "internal Tpetra logic error. Please contact Tpetra team.\n");
#endif
        }
        block[subCol*blockSize + subRow] = *rowval;
      }
    }
  }
#ifdef HAVE_TPETRA_DEBUG
  TEST_FOR_EXCEPTION( numAllocatedBlocks != blocks.size(), std::logic_error, errstr << "internal Tpetra logic error. Please contact Tpetra team.\n");
#endif
  out_block_entries = arcp<ArrayRCP<Scalar> >( numAllocatedBlocks );
  out_block_indices = arcp<LocalOrdinal>     ( numAllocatedBlocks );
  typedef typename std::map<LocalOrdinal,ArrayRCP<Scalar> >::iterator MapIter;
  MapIter block = blocks.begin();
  const MapIter blockend = blocks.end();
  for (size_t b=0; block != blockend; ++block, ++b) 
  {
    out_block_indices[b] = (*block).first;
    out_block_entries[b] = (*block).second;
  }
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra::Ext namespace!
//

#define TPETRAEXT_BLOCKEXTRACTION_INSTANT(SCALAR,LO,GO,NODE)   \
  template void                                                \
  extractBlockDiagonals(const RowMatrix<SCALAR,LO,GO,NODE> &,  \
                        const ArrayView<const LO> &,  \
                              ArrayRCP<SCALAR> &,     \
                              ArrayRCP<LO> &);        \
  \
  template void                                                \
  extractBlockDiagonals(const RowMatrix<SCALAR,LO,GO,NODE> &,  \
                        const BlockMap<LO,GO,NODE> &,          \
                              ArrayRCP<SCALAR> &,     \
                              ArrayRCP<LO> &);        \
  \
  template void                                         \
  extractBlockRow(LO localBlockRow,                     \
                  const RowMatrix<SCALAR,LO,GO,NODE> &, \
                  const BlockMap<LO,GO,NODE>         &, \
                  const BlockMap<LO,GO,NODE>         &, \
                  ArrayRCP<ArrayRCP<SCALAR> >        &, \
                  ArrayRCP<LO>                       &);


#endif // TPETRAEXT_BLOCKEXTRACTION_DEF_HPP
