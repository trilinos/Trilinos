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

#ifndef TPETRA_VBRMATRIX_DEF_HPP
#define TPETRA_VBRMATRIX_DEF_HPP

#include <Kokkos_NodeHelpers.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#ifdef DOXYGEN_USE_ONLY
#include "Tpetra_VbrMatrix_decl.hpp"
#endif

namespace Tpetra {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::VbrMatrix(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype)
 : blkRowMap_(blkRowMap),
   blkColMap_(Teuchos::null),
   blkDomainMap_(Teuchos::null),
   blkRangeMap_(Teuchos::null),
   blkGraph_(),
   lclMatrix_(blkRowMap->getNodeNumBlocks(), blkRowMap->getPointMap()->getNode()),
   pbuf_values1D_(),
   lclMatVec_(blkRowMap->getPointMap()->getNode()),
   col_ind_2D_global_(),
   col_ind_2D_local_(),
   pbuf_values2D_()
{
  //The graph of this VBR matrix will be a CrsGraph where each entry in the graph
  //corresponds to a block-entry in the matrix.
  //That is, you can think of a VBR matrix as a Crs matrix of dense submatrices...

  global_size_t numGlobalElems = Teuchos::OrdinalTraits<global_size_t>::invalid();
  GlobalOrdinal indexBase = blkRowMap->getPointMap()->getIndexBase();
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm = blkRowMap->getPointMap()->getComm();
  const Teuchos::RCP<Node>& node = blkRowMap->getPointMap()->getNode();

  //Create a point-entry map where each point
  //corresponds to a block in our block-row-map:
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > blkPointMap = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElems, blkRowMap->getBlockIDs(), indexBase, comm, node));
  blkGraph_ = Teuchos::rcp(new CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(blkPointMap, maxNumEntriesPerRow, pftype));
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~VbrMatrix()
{
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getDomainMap() const
{
  return blkDomainMap_->getPointMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRangeMap() const
{
  return blkRangeMap_->getPointMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(
         const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
               Teuchos::ETransp mode,
               Scalar alpha,
               Scalar beta) const
{
  throw std::runtime_error("Tpetra::VbrMatrix::apply not yet implemented!");
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
bool
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApply() const
{
  return false;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getBlockRowMap() const
{
  return blkRowMap_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getPointRowMap() const
{
  return blkRowMap_->getPointMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getBlockColMap() const
{
  return blkColMap_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getPointColMap() const
{
  return blkColMap_->getPointMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
Teuchos::ArrayRCP<Scalar>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalBlockEntryView(
    GlobalOrdinal globalBlockRow,
    GlobalOrdinal globalBlockCol,
    size_t rowsPerBlock,
    size_t colsPerBlock)
{
  //This private method returns a block-entry (as an ArrayRCP),
  //creating/allocating the block-entry if it doesn't already exist.

  Teuchos::ArrayRCP<Scalar> blockEntryView;

  LocalOrdinal numBlockRows = blkGraph_->getNodeNumRows();
  Teuchos::RCP<Node> node = getNode();

  if (pbuf_values1D_ != Teuchos::null) {
    //still need to implement the packed (storage-optimized) stuff...
  }
  else {
    if (pbuf_values2D_.size() == 0) {
      pbuf_values2D_.resize(numBlockRows);
      col_ind_2D_global_.resize(numBlockRows);
    }
  }

  LocalOrdinal localBlockRow = blkRowMap_->getLocalBlockID(globalBlockRow);

  //this is essentially a range-check for globalBlockRow:
  TEST_FOR_EXCEPTION( localBlockRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
     std::runtime_error,
     "Tpetra::VbrMatrix::getGlobalBlockEntryView, globalBlockRow not local.");

  MapGlobalArrayRCP& blkrow = col_ind_2D_global_[localBlockRow];
  typename MapGlobalArrayRCP::iterator col_iter = blkrow.find(globalBlockCol);

  if (col_iter != blkrow.end()) {
    blockEntryView = col_iter->second;
  }
  else {
    //block-entry doesn't already exist, so we will create it.

    //make sure block-size is specified:
    TEST_FOR_EXCEPTION(rowsPerBlock==0 || colsPerBlock==0, std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockEntryView ERROR: creating block-entry, but rowsPerBlock and/or colsPerBlock is 0.");

    size_t blockSize = rowsPerBlock*colsPerBlock;
    blockEntryView = node->template allocBuffer<Scalar>(blockSize);
    pbuf_values2D_[globalBlockRow].push_back(blockEntryView);
    blkrow.insert(std::make_pair(globalBlockCol, blockEntryView));
    blkGraph_->insertGlobalIndices(globalBlockRow, Teuchos::ArrayView<GlobalOrdinal>(&globalBlockCol, 1));
  }

  return blockEntryView;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
Teuchos::RCP<Node>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNode()
{
  return blkRowMap_->getPointMap()->getNode();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry = getGlobalBlockEntryView(globalBlockRow,globalBlockCol, blockEntry.numRows(), blockEntry.numCols());

  //now copy the incoming block-entry into internal storage:
  size_t offset = 0;
  for(GlobalOrdinal col=0; col<blockEntry.numCols(); ++col) {
    for(GlobalOrdinal row=0; row<blockEntry.numRows(); ++row) {
      internalBlockEntry[offset++] = blockEntry[col][row];
    }
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol,
                           LocalOrdinal& numPtRows, LocalOrdinal& numPtCols,
                           Teuchos::ArrayRCP<const Scalar>& blockEntry)
{
  //obtain pointer to internal storage.
  //(This will throw if it doesn't already exist, since we're not specifying
  //the dimensions of the block-entry.)
  blockEntry = getGlobalBlockEntryView(globalBlockRow, globalBlockCol);

  LocalOrdinal localBlockID = blkRowMap_->getLocalBlockID(globalBlockRow);
  numPtRows = blkRowMap_->getBlockSize(localBlockID);

  TEST_FOR_EXCEPTION(numPtRows == 0, std::runtime_error,
    "Tpetra::VbrMatrix::getGlobalBlockEntry ERROR, numPtRows == 0.");

  numPtCols = blockEntry.size() / numPtRows;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap, OptimizeOption opt)
{
  blkGraph_->fillComplete(opt);
  throw std::runtime_error("Tpetra::VbrMatrix::fillComplete(domainMap,rangeMap) not yet implemented!");
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillComplete(OptimizeOption opt)
{
  fillComplete(getBlockRowMap(), getBlockRowMap(), opt);
}

}//namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_VBRMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class VbrMatrix< SCALAR , LO , GO , NODE >;

#endif //TPETRA_VBRMATRIX_DEF_HPP

