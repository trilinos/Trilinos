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

#include "Tpetra_CrsMatrix.hpp"

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void Tpetra::Ext::extractBlockDiagonals(
          const RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> &matrix, 
          const Teuchos::ArrayView<const LocalOrdinal> &block_sizes,
          Teuchos::ArrayRCP<Scalar>       &out_block_diagonals,
          Teuchos::ArrayRCP<LocalOrdinal> &out_block_offsets)
{
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

