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

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
Tpetra::Ext::extractBlockDiagonals(
          const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & matrix, 
          const Teuchos::ArrayView<const LocalOrdinal>            & first_points,
          Teuchos::ArrayRCP<Scalar>                               & out_diags,
          Teuchos::ArrayRCP<LocalOrdinal>                         & out_offsets)
{
  // block meta-data
  const int numBlocks = (int)first_points.size();
  const size_t numRows = matrix.getNodeNumRows();
  TEST_FOR_EXCEPTION(numRows != 0 && numBlocks == 0, std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): specified zero blocks for a matrix with more than zero rows.");
  TEST_FOR_EXCEPTION(numBlocks > 0 && first_points[0] != 0, std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): first point of first block must be zero.");
  TEST_FOR_EXCEPTION(matrix.isFillComplete() == false, std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): matrix must be fill-complete.");
  // INVARIANT: first_points[i] - first_points[i-1] = block_size[i]
  // therefore, block_size[numBlocks-1] = first_points[numBlocks] - first_points[numBlocks-1]
  // however, first_points[numBlocks] does not exist. therefore, the size of the last block must be inferred by 
  // subtracting the number of rows from the sum of the size of all other blocks
  int alloc_size      = 0, 
      last_block_size = numRows;
  if (numBlocks) {
    out_offsets = arcp<LocalOrdinal>(numBlocks);
  }
  for (int b=0; b < numBlocks-1; ++b) 
  {
    const int block_size_b = (int)first_points[b+1] - (int)first_points[b];
    TEST_FOR_EXCEPTION(block_size_b < 0, std::runtime_error,
        "Tpetra::Ext::extractBlockDiagonals(): first points are not coherent.");
    last_block_size -= block_size_b;
    out_offsets[b]   = alloc_size;
    alloc_size      += block_size_b*block_size_b;
  }
  out_offsets[numBlocks-1] = alloc_size;
  alloc_size += last_block_size*last_block_size;
  TEST_FOR_EXCEPTION( last_block_size < 0, std::runtime_error, 
      "Tpetra::Ext::extractBlockDiagonals(): specified blocks are not compatible with specified matrix (blocks are too large).");
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
    const int first_block = 0,
               last_block = numBlocks-1;
    int                 b = first_block-1, 
                   subrow = 0, 
             block_size_b = 0;
    typename ArrayRCP<Scalar>::iterator block = out_diags.end();
    // loop over all local rows
    for (LocalOrdinal lrow = first_row; lrow <= last_row; ++lrow)
    {
      // the while loop accounts for blocks of size zero
      while (subrow == block_size_b) 
      {
        // we busted the block, move to the next
        b += 1;
        block_size_b = (b == last_block) ? last_block_size 
                                         : (first_points[b+1] - first_points[b]);
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
  ArrayRCP<const LocalOrdinal> block_firsts_arcp_is_too_long = block_map.getNodeFirstPointInBlocks();
  ArrayView<const LocalOrdinal> block_firsts = block_firsts_arcp_is_too_long(0,block_firsts_arcp_is_too_long.size()-1);
  try {
    extractBlockDiagonals( matrix, block_firsts, out_diags, out_offsets );
  }
  catch (std::exception &e) {
    TEST_FOR_EXCEPTION(true, std::runtime_error, 
        "Tpetra::Ext::extractBlockDiagonals(RowMatrix,BlockMap,...) caught exception:\n\n" << e.what());
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
                        const Teuchos::ArrayView<const LO> &,  \
                              Teuchos::ArrayRCP<SCALAR> &,     \
                              Teuchos::ArrayRCP<LO> &);        \
  \
  template void                                                \
  extractBlockDiagonals(const RowMatrix<SCALAR,LO,GO,NODE> &,  \
                        const BlockMap<LO,GO,NODE> &,          \
                              Teuchos::ArrayRCP<SCALAR> &,     \
                              Teuchos::ArrayRCP<LO> &);


#endif // TPETRAEXT_BLOCKEXTRACTION_DEF_HPP
