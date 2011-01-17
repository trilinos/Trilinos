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
void Tpetra::Ext::extractBlockDiagonals(
          const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &matrix, 
          const Teuchos::ArrayView<const LocalOrdinal>            &block_sizes,
          Teuchos::ArrayRCP<Scalar>                               &out_block_diagonals,
          Teuchos::ArrayRCP<LocalOrdinal>                         &out_block_offsets)
{
  // block meta-data
  const int numBlocks = (int)block_sizes.size();
  if (numBlocks) {
    out_block_offsets = arcp<LocalOrdinal>(numBlocks);
  }
  int blockSum = 0;
  int allocSize = 0;
  for (int b=0; b != (int)numBlocks; ++b) {
    out_block_offsets[b] = allocSize;
    allocSize += block_sizes[b]*block_sizes[b];
    blockSum  += block_sizes[b];
  }
  TEST_FOR_EXCEPTION( (size_t)blockSum != matrix.getNodeNumRows(), std::runtime_error,
      "Tpetra::Ext::extractBlockDiagonals(): block_sizes must sum to number of local matrix rows.");
  if (allocSize) {
    out_block_diagonals = arcp<Scalar>(allocSize);
    std::fill( out_block_diagonals.begin(), out_block_diagonals.end(), ScalarTraits<Scalar>::zero() );
  }
  // extract blocks
  if (allocSize) {
    const LocalOrdinal first_row = matrix.getRowMap()->getMinLocalIndex(),
                        last_row = matrix.getRowMap()->getMaxLocalIndex();
    ArrayView<const Scalar> rowvals;
    ArrayView<const LocalOrdinal> rowinds;
    int b = 0, subrow = 0, coloffset = 0;
    typename ArrayRCP<Scalar>::iterator block;
    block = out_block_diagonals.persistingView( out_block_offsets[b], block_sizes[b]*block_sizes[b] ).begin();
    // loop over all local rows
    for (LocalOrdinal lrow = first_row; lrow <= last_row; ++lrow)
    {
      if (subrow == block_sizes[b]) {
        // we busted the block, move to the next
        coloffset += block_sizes[b];
        ++b;
        block = out_block_diagonals.persistingView( out_block_offsets[b], block_sizes[b]*block_sizes[b] ).begin();
        subrow = 0;
      }
      // extract the row, put the members of the block diagonal into the current block
      matrix.getLocalRowView(lrow, rowinds, rowvals);
      for (int k=0; k < (int)rowinds.size(); ++k) {
        const int subcol = rowinds[k] - coloffset;
        if (subcol >= 0 && subcol < block_sizes[b]) {
          block[subcol*block_sizes[b] + subrow] += rowvals[k];
        }
      }
      ++subrow;
    }
  }
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra::Ext namespace!
//

#define TPETRAEXT_BLOCKEXTRACTION_INSTANT(SCALAR,LO,GO,NODE) \
  template void extractBlockDiagonals(const RowMatrix<SCALAR,LO,GO,NODE> &,  \
                                      const Teuchos::ArrayView<const LO> &,  \
                                      Teuchos::ArrayRCP<SCALAR>          &,  \
                                      Teuchos::ArrayRCP<LO>              &);

#endif // TPETRAEXT_BLOCKEXTRACTION_DEF_HPP
