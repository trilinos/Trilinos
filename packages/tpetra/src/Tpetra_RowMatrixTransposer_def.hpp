//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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

#ifndef TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP
#define TPETRA_ROWMATRIXTRANSPOSER_DEF_HPP

#include "Tpetra_Map.hpp"
#include "Teuchos_DefaultSerialComm.hpp"
#ifdef DOXYGEN_USE_ONLY
  // #include "Tpetra_RowMatrixtransposer_decl.hpp"
#endif

namespace Tpetra {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatOps>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::RowMatrixTransposer(const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& origMatrix)
  : origMatrix_(origMatrix), comm_(origMatrix.getComm()), indexBase_(origMatrix_.getIndexBase()) {}

template<class Scalar, 
  class LocalOrdinal, 
  class GlobalOrdinal, 
  class Node, 
  class SpMatOps>
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::~RowMatrixTransposer() {}

template<class Scalar, 
  class LocalOrdinal,
  class GlobalOrdinal, 
  class Node, 
  class SpMatOps>
RCP<CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps> >
RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>::createTranspose (
  const OptimizeOption optimizeTranspose,
  Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > transposeRowMap)
{
  optimizeTranspose_ = optimizeTranspose;
  transposeMatrix_ = rcp(new CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, SpMatOps>(
    transposeRowMap == null? origMatrix_.getDomainMap() : transposeRowMap, 0));
  ArrayView<const LocalOrdinal> localIndicies;
  ArrayView<const Scalar> localValues;
  ArrayView<const GlobalOrdinal> myGlobalRows = origMatrix_.getRowMap()->getNodeElementList();
  size_t numEntriesInRow;
  Array<GlobalOrdinal> rowNum(1);
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > origColMap = origMatrix_.getColMap();
  for(
    size_t i = Teuchos::OrdinalTraits<size_t>::zero();
    i < origMatrix_.getNodeNumRows();
    ++i
  )
  {
    rowNum[0] = origMatrix_.getRowMap()->getGlobalElement(i);
    numEntriesInRow = origMatrix_.getNumEntriesInLocalRow(i);
    origMatrix_.getLocalRowView(i, localIndicies, localValues);
    for(size_t j=0; j<numEntriesInRow; ++j){
      transposeMatrix_->insertGlobalValues(origColMap->getGlobalElement(localIndicies[j]), rowNum(0,1), localValues(j,1));
    }
  }

  transposeMatrix_->fillComplete(origMatrix_.getRangeMap(), origMatrix_.getDomainMap(), optimizeTranspose_);

  return transposeMatrix_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_ROWMATRIXTRANSPOSE_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class RowMatrixTransposer< SCALAR, LO , GO , NODE >;


}

#endif
