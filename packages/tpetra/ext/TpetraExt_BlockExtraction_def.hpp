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
Teuchos::ArrayRCP<Scalar>
Tpetra::Ext::extractBlockDiagonals(
          const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &matrix, 
          const Teuchos::ArrayView<const LocalOrdinal> &block_offsets)
{
  // block meta-data
  const int numBlocks = (int)block_offsets.size();
  TEST_FOR_EXCEPTION(numBlocks > (int)matrix.getNodeNumRows(), std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): specified block offsets are not compatible with specified matrix.");
  TEST_FOR_EXCEPTION(block_offsets[0] != 0, std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): offset to first block must be zero.");
  // INVARIANT: block_offsets[i] - block_offsets[i-1] = block_size[i-1]**2
  // therefore, block_size[numBlocks-1]**2 = block_offsets[numBlocks] - block_offsets[numBlocks-1]
  // however, block_offsets[numBlocks] does not exist. therefore, the size of the last block must be inferred by 
  // subtracting the number of rows from the sum of the size of all other blocks
  const size_t numRows = matrix.getNodeNumRows();
  int alloc_size      = 0, 
      last_block_size = numRows;
  for (int b=0; b < numBlocks-1; ++b) {
    const int block_size_b = ScalarTraits<int>::squareroot( block_offsets[b+1] - block_offsets[b] );
    last_block_size -= block_size_b;
    alloc_size      += block_size_b*block_size_b;
  }
  alloc_size += last_block_size*last_block_size;
  TEST_FOR_EXCEPTION( last_block_size < 0, std::runtime_error, 
      "Tpetra::Ext::extractBlockDiagonals(): specified block offsets are not compatible with specified matrix.");
  ArrayRCP<Scalar> out_block_diagonals;
  if (alloc_size) {
    out_block_diagonals = arcp<Scalar>(alloc_size);
    std::fill( out_block_diagonals.begin(), out_block_diagonals.end(), ScalarTraits<Scalar>::zero() );
  }
  // extract blocks
  if (alloc_size) {
    const LocalOrdinal first_row = matrix.getRowMap()->getMinLocalIndex(),
                        last_row = matrix.getRowMap()->getMaxLocalIndex();
    ArrayView<const Scalar> rowvals;
    ArrayView<const LocalOrdinal> rowinds;
    int b = 0, subrow = 0, coloffset = 0, block_size_b;
    block_size_b = ScalarTraits<int>::squareroot( block_offsets[b+1] - block_offsets[b] );
    typename ArrayRCP<Scalar>::iterator block;
    block = out_block_diagonals.persistingView( block_offsets[b], block_size_b*block_size_b ).begin();
    // loop over all local rows
    for (LocalOrdinal lrow = first_row; lrow <= last_row; ++lrow)
    {
      if (subrow == block_size_b) {
        // we busted the block, move to the next
        coloffset += block_size_b;
        ++b;
        block_size_b = (b == numBlocks-1) ? last_block_size : ScalarTraits<int>::squareroot( block_offsets[b+1] - block_offsets[b] );
        block = out_block_diagonals.persistingView( block_offsets[b], block_size_b*block_size_b ).begin();
        subrow = 0;
      }
      // extract the row, put the members of the block diagonal into the current block
      matrix.getLocalRowView(lrow, rowinds, rowvals);
      for (int k=0; k < (int)rowinds.size(); ++k) {
        const int subcol = rowinds[k] - coloffset;
        if (subcol >= 0 && subcol < block_size_b) {
          block[subcol*block_size_b + subrow] += rowvals[k];
        }
      }
      ++subrow;
    }
  }
  return out_block_diagonals;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra::Ext namespace!
//

#define TPETRAEXT_BLOCKEXTRACTION_INSTANT(SCALAR,LO,GO,NODE)   \
  template Teuchos::ArrayRCP<SCALAR>                           \
  extractBlockDiagonals(const RowMatrix<SCALAR,LO,GO,NODE> &,  \
                        const Teuchos::ArrayView<const LO> &);

#endif // TPETRAEXT_BLOCKEXTRACTION_DEF_HPP
