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

#ifndef TPETRAEXT_BLOCKEXTRACTION_DECL_HPP
#define TPETRAEXT_BLOCKEXTRACTION_DECL_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_BlockMap.hpp"

/*! \file TpetraExt_BlockExtraction_decl.hpp 

    Methods for block extraction of data from Tpetra objects.
 */

namespace Tpetra {
  namespace Ext {

    /** \brief Extracts the block diagonals from a RowMatrix into contiguous, host storage.

        This method does not initiate any communication, and therefore can be called safely on a single node.

        \param in matrix - The sparse matrix source for the diagonals.
        \param in first_points - A list of ordinals, where <tt>first_points[b]</tt> indicates the first local element in the <tt>b</tt>-th block.
        \param out out_diags - A reference-counted array, containing the block diagonals in contiguous storage.
        \param out out_offsets - A reference-counted array, indicating the offset to reach each block in \c out_diags.

        \pre - <tt>first_points[0] == 0</tt>
        \pre - elements in \c first_points are non-decreasing
        \pre - <tt>matrix.isFillComplete() == true</tt>

        \post - <tt>out_offsets.size() == block_sizes.size()</tt>
        \post - the non-trivial <tt>b</tt>-th block is stored contiguously (column-major) in <tt>out_diags( out_offsets[b], block_size_b )</tt>, where \f$block\_size\_b = ( first\_points\[b+1\] - first\_points\[b\] )^2 \f$.

        \relatesalso Tpetra::CrsMatrix
      */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void
    extractBlockDiagonals(const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &matrix, 
                          const Teuchos::ArrayView<const LocalOrdinal> &first_points,
                          Teuchos::ArrayRCP<Scalar>       &out_diags,
                          Teuchos::ArrayRCP<LocalOrdinal> &out_offsets);

    /** \brief Extracts the block diagonals from a RowMatrix into contiguous, host storage using a BlockMap.

        This method does not initiate any communication, and therefore can be called safely on a single node.

        \param in matrix - The sparse matrix source for the diagonals.
        \param in block_map - A BlockMap describing the blocks
        \param out out_diags - A reference-counted array, containing the block diagonals in contiguous storage.
        \param out out_offsets - A reference-counted array, indicating the offset to reach each block in \c out_diags.

        \pre - <tt>block_map.getPointMap().isCompatible( matrix.getRowMap() )</tt>
        \pre - <tt>matrix.isFillComplete() == true</tt>

        \post - <tt>out_offsets.size() == block_sizes.size()</tt>
        \post - the non-trivial <tt>b</tt>-th block is stored contiguously (column-major) in <tt>out_diags( out_offsets[b], block_map.getLocalBlockSize(b) )</tt>
      
        Calls extractBlockDiagonals().

        \relatesalso Tpetra::CrsMatrix
      */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void
    extractBlockDiagonals(const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &matrix, 
                          const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>         &block_map,
                          Teuchos::ArrayRCP<Scalar>       &out_diags,
                          Teuchos::ArrayRCP<LocalOrdinal> &out_offsets);

    /** \brief Extracts block elements from a RowMatrix into contiguous, host storage.

        This method does not initiate any communication, and therefore can be called safely on a single node.

        \param in block_row     - The block row to be extracted
        \param in block_row_map - A BlockMap describing the row blocks
        \param in block_col_map - A BlockMap describing the column blocks
        \param in matrix - The sparse matrix source for the diagonals.
        \param out block_entries - The block entries
        \param out block_indices - The indices for the block entries

        \pre - <tt>block_row_map.getPointMap().isCompatible( matrix.getRowMap() )</tt>
        \pre - <tt>block_row_map.getPointMap().isCompatible( matrix.getRowMap() )</tt>
        \pre - <tt>block_col_map.getPointMap().isCompatible( matrix.getColMap() )</tt>
        \pre - <tt>block_col_map.getPointMap().isCompatible( matrix.getColMap() )</tt>
        \pre - <tt>matrix.isFillComplete() == true</tt>

        \relatesalso Tpetra::CrsMatrix
      */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void
    extractBlockRow(LocalOrdinal localBlockRow,
                    const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &matrix, 
                    const BlockMap<LocalOrdinal,GlobalOrdinal,Node> &block_row_map,
                    const BlockMap<LocalOrdinal,GlobalOrdinal,Node> &block_col_map,
                    ArrayRCP<ArrayRCP<Scalar> >      &out_block_entries,
                    ArrayRCP<LocalOrdinal>           &out_block_indices);

  }
}

#endif // TPETRAEXT_BLOCKEXTRACTION_DECL_HPP
