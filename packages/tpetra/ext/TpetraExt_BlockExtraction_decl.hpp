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

namespace Tpetra {
  namespace Ext {

    /** \brief Extracts the block diagonals from a RowMatrix into contiguous, host storage.

        \param in matrix - The sparse matrix
        \param in block_sizes - A list of ordinals indicating the sizes defining the block diagonals
        \param out out_block_diagonals - a host-allocated buffer for the extracted diagonal blocks
        \param out out_block_offsets  - a convenience array of offsets into \c out_block_diagonals

        \pre The sum of the block sizes must be equal to the number of local rows in the matrix
        \post - The sparse entries in the <tt>i</tt>-th block diagonal are stored in <tt>out_block_diagonals( out_block_offsets[i] , block_sizes[i]*block_sizes[i] )</tt>
        \post - <tt>out_block_offsets.size() == block_sizes.size()</tt>
      */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
    void extractBlockDiagonals(const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &matrix, 
                               const Teuchos::ArrayView<const LocalOrdinal> &block_sizes,
                               Teuchos::ArrayRCP<Scalar>       &out_block_diagonals,
                               Teuchos::ArrayRCP<LocalOrdinal> &out_block_offsets);
  }
}

#endif // TPETRAEXT_BLOCKEXTRACTION_DECL_HPP
