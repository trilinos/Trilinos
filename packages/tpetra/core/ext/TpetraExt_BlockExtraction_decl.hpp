// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRAEXT_BLOCKEXTRACTION_DECL_HPP
#define TPETRAEXT_BLOCKEXTRACTION_DECL_HPP

/// \file TpetraExt_BlockExtraction_decl.hpp
/// \brief Methods for block extraction of data from Tpetra objects.

#include <Tpetra_ConfigDefs.hpp>

#ifndef HAVE_TPETRA_CLASSIC_VBR
#  error "It is an error to include this file if VBR (variable-block-size) sparse matrix support is disabled in Tpetra.  If you would like to enable VBR support, please reconfigure Trilinos with the CMake option Tpetra_ENABLE_VBR set to ON, and rebuild Trilinos."
#else

#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_BlockMap.hpp"

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

/**
  \example CrsMatrix_BlockExtraction.cpp
  An example for using the block extraction methods Tpetra::Ext::extractBlockDiagonals() and Tpetra::Ext::extractBlockRow().
 */

#endif // ! HAVE_TPETRA_CLASSIC_VBR
#endif // ! TPETRAEXT_BLOCKEXTRACTION_DECL_HPP
