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
#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Tpetra_Vector.hpp>

#ifdef DOXYGEN_USE_ONLY
#include "Tpetra_VbrMatrix_decl.hpp"
#endif

namespace Tpetra {

//-----------------------------------------------------------------
template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
convertBlockMapToPointMap(const Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockMap)
{
  global_size_t numGlobalElems = Teuchos::OrdinalTraits<global_size_t>::invalid();
  GlobalOrdinal indexBase = blockMap->getPointMap()->getIndexBase();
  const Teuchos::RCP<const Teuchos::Comm<int> >& comm = blockMap->getPointMap()->getComm();
  const Teuchos::RCP<Node>& node = blockMap->getPointMap()->getNode();

  //Create a point-entry map where each point
  //corresponds to a block in the block-map:
  return Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElems, blockMap->getBlockIDs(), indexBase, comm, node));
}

//-----------------------------------------------------------------
template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >
makeBlockColumnMap(const Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockMap, const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& ptColMap)
{
  global_size_t numGlobalBlocks = Teuchos::OrdinalTraits<global_size_t>::invalid();
  Teuchos::ArrayView<const GlobalOrdinal> blockIDs = ptColMap->getNodeElementList();
  Teuchos::ArrayView<const LocalOrdinal> blockSizes = blockMap->getBlockSizes();

  //We will create a block-map where each block corresponds to a point in
  //the input-map 'ptColMap', and where each block's size is obtained from
  //the input-block-map 'blockMap'.

  //if we are running on a single processor, then it's easy because
  //we know that blockMap is distributed the same as ptColMap:
  int numProcessors = ptColMap->getComm()->getSize();
  if (numProcessors == 1) {
    if (blockMap->isBlockSizeConstant()) {
    return Teuchos::rcp(new Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>(
        numGlobalBlocks, blockIDs, blockSizes,
        ptColMap->getIndexBase(), ptColMap->getComm(), ptColMap->getNode()
    ));
    }
    else {
    return Teuchos::rcp(new Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>(
        numGlobalBlocks, blockIDs, blockSizes,
        ptColMap->getIndexBase(), ptColMap->getComm(), ptColMap->getNode()
    ));
    }
  }

  //if we get to here, we're running on multiple processors, which means
  //blockMap is possibly (probably) not distributed the same as ptColMap,
  //so we have to do some communication to get the block-sizes from blockMap
  //corresponding to the blockIDs we got from ptColMap.
  //
  //I think the simplest way to do this is to create vectors where the values
  //are block-sizes, and import one to the other.
  typedef Tpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> OrdVec;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > blkPtMap =
      convertBlockMapToPointMap(blockMap);

  OrdVec source(blkPtMap, blockSizes);
  OrdVec target(ptColMap);
  Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> importer(blkPtMap, ptColMap);
  target.doImport(source, importer, REPLACE);
  blockSizes = target.get1dView()();

  //now we can create our block-column-map:
  return Teuchos::rcp(new Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>(
    numGlobalBlocks, blockIDs, blockSizes,
    ptColMap->getIndexBase(), ptColMap->getComm(), ptColMap->getNode()
  ));
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::VbrMatrix(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype)
 : blkRowMap_(blkRowMap),
   blkColMap_(Teuchos::null),
   blkDomainMap_(Teuchos::null),
   blkRangeMap_(Teuchos::null),
   blkGraph_(),
   lclMatrix_(blkRowMap->getNodeNumBlocks(), blkRowMap->getPointMap()->getNode()),
   pbuf_values1D_(),
   pbuf_rptr_(),
   pbuf_cptr_(),
   pbuf_bptr_(),
   pbuf_bindx_(),
   pbuf_indx_(),
   lclMatVec_(blkRowMap->getPointMap()->getNode()),
   col_ind_2D_global_(Teuchos::rcp(new Teuchos::Array<std::map<GlobalOrdinal,Teuchos::ArrayRCP<Scalar> > >)),
   col_ind_2D_local_(Teuchos::rcp(new Teuchos::Array<std::map<LocalOrdinal,Teuchos::ArrayRCP<Scalar> > >)),
   pbuf_values2D_(Teuchos::rcp(new Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > >)),
   is_fill_completed_(false),
   is_storage_optimized_(false)
{
  //The graph of this VBR matrix will be a CrsGraph where each entry in the
  //graph corresponds to a block-entry in the matrix.
  //That is, you can think of a VBR matrix as a Crs matrix of dense
  //submatrices...

  //First create a point-entry map where each point corresponds to
  //a block in our block-row-map:
  Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >
      blkPointMap = convertBlockMapToPointMap(blkRowMap);

  //now create the graph:
  blkGraph_ = Teuchos::rcp(new CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(blkPointMap, maxNumEntriesPerRow, pftype));
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::~VbrMatrix()
{
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getDomainMap() const
{
  return blkDomainMap_->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getRangeMap() const
{
  return blkRangeMap_->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::apply(
         const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
               Teuchos::ETransp mode,
               Scalar alpha,
               Scalar beta) const
{
  const Kokkos::MultiVector<Scalar,Node> *lclX = &X.getLocalMV();
  Kokkos::MultiVector<Scalar,Node>        *lclY = &Y.getLocalMVNonConst();

  lclMatVec_.template multiply<Scalar,Scalar>(mode,alpha,*lclX,beta,*lclY);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
bool
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::hasTransposeApply() const
{
  return false;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getBlockRowMap() const
{
  return blkRowMap_;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getPointRowMap() const
{
  return blkRowMap_->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getBlockColMap() const
{
  return blkColMap_;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getPointColMap() const
{
  return blkColMap_->getPointMap();
}

//-------------------------------------------------------------------
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
    throw std::runtime_error("Tpetra::VbrMatrix ERROR, packed storage not yet implemented.");
  }
  else {
    if (pbuf_values2D_->size() == 0) {
      pbuf_values2D_->resize(numBlockRows);
      col_ind_2D_global_->resize(numBlockRows);
    }
  }

  LocalOrdinal localBlockRow = blkRowMap_->getLocalBlockID(globalBlockRow);

  //this is essentially a range-check for globalBlockRow:
  TEST_FOR_EXCEPTION( localBlockRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
     std::runtime_error,
     "Tpetra::VbrMatrix::getGlobalBlockEntryView, globalBlockRow not local.");

  MapGlobalArrayRCP& blkrow = (*col_ind_2D_global_)[localBlockRow];
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
    (*pbuf_values2D_)[localBlockRow].push_back(blockEntryView);
    blkrow.insert(std::make_pair(globalBlockCol, blockEntryView));
    blkGraph_->insertGlobalIndices(globalBlockRow, Teuchos::ArrayView<GlobalOrdinal>(&globalBlockCol, 1));
  }

  return blockEntryView;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
Teuchos::RCP<Node>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getNode()
{
  return blkRowMap_->getPointMap()->getNode();
}

//-------------------------------------------------------------------
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

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry = getGlobalBlockEntryView(globalBlockRow,globalBlockCol, blkRowSize, blkColSize);

  LocalOrdinal blk_size = blockEntry.size();
  TEST_FOR_EXCEPTION(blkColSize*LDA > blk_size, std::runtime_error,
    "Tpetra::VbrMatrix::setGlobalBlockEntry ERROR, inconsistent block-entry sizes.");

  //copy the incoming block-entry into internal storage:
  size_t offset = 0;
  size_t src_offset = 0;
  for(GlobalOrdinal col=0; col<blkColSize; ++col) {
    for(GlobalOrdinal row=0; row<blkRowSize; ++row) {
      internalBlockEntry[offset++] = blockEntry[src_offset + row];
    }
    src_offset += LDA;
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::getGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol,
                           LocalOrdinal& numPtRows, LocalOrdinal& numPtCols,
                           Teuchos::ArrayRCP<const Scalar>& blockEntry)
{
  //obtain pointer to internal storage.
  //(This will throw if it doesn't already exist, since we're not specifying
  //the default arguments that are the dimensions of the block-entry.)
  blockEntry = getGlobalBlockEntryView(globalBlockRow, globalBlockCol);

  LocalOrdinal localBlockID = blkRowMap_->getLocalBlockID(globalBlockRow);
  numPtRows = blkRowMap_->getBlockSize(localBlockID);

  TEST_FOR_EXCEPTION(numPtRows == 0, std::runtime_error,
    "Tpetra::VbrMatrix::getGlobalBlockEntry ERROR, numPtRows == 0.");

  numPtCols = blockEntry.size() / numPtRows;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::optimizeStorage()
{
  if (is_storage_optimized_ == true) return;

  size_t num_block_rows = blkGraph_->getNodeNumRows();
  size_t num_block_cols = blkGraph_->getNodeNumCols();
  size_t num_block_nonzeros = blkGraph_->getNodeNumEntries();

  //need to count the number of point-entries:
  size_t num_point_entries = 0;
  typedef typename Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > >::size_type Tsize_t;
  for(Tsize_t r=0; r<pbuf_values2D_->size(); ++r) {
    Tsize_t rlen = (*pbuf_values2D_)[r].size();
    for(Tsize_t c=0; c<rlen; ++c) {
      num_point_entries += (*pbuf_values2D_)[r][c].size();
    }
  }

  const Teuchos::RCP<Node>& node = blkRowMap_->getPointMap()->getNode();

  //Don't change these allocation sizes unless you know what you're doing!
  //The '+1's are in the right place. (Trust me, I'm a trained professional.)
  pbuf_rptr_ = node->template allocBuffer<LocalOrdinal>(num_block_rows+1);
  pbuf_cptr_ = node->template allocBuffer<LocalOrdinal>(num_block_cols+1);
  pbuf_bptr_ = node->template allocBuffer<LocalOrdinal>(num_block_rows+1);
  pbuf_bindx_= node->template allocBuffer<LocalOrdinal>(num_block_nonzeros);
  pbuf_indx_ = node->template allocBuffer<LocalOrdinal>(num_block_nonzeros+1);
  pbuf_values1D_ = node->template allocBuffer<Scalar>(num_point_entries);

  size_t roffset = 0;
  size_t pt_r_offset = 0;
  size_t ioffset = 0;
  size_t offset = 0;
  Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > >& pbuf_vals2D = *pbuf_values2D_;
  for(Tsize_t r=0; r<pbuf_vals2D.size(); ++r) {
    pbuf_rptr_[r] = pt_r_offset;
    LocalOrdinal rsize = blkRowMap_->getBlockSize(r);
    pt_r_offset += rsize;

    pbuf_bptr_[r] = roffset;
    Tsize_t rlen = pbuf_vals2D[r].size();
    roffset += rlen;

    Teuchos::ArrayRCP<const LocalOrdinal> blk_row_inds = blkGraph_->getLocalRowView(r);

    for(Tsize_t c=0; c<rlen; ++c) {
      pbuf_bindx_[ioffset] = blk_row_inds[c];
      pbuf_indx_[ioffset++] = offset;

      Tsize_t blkSize = pbuf_vals2D[r][c].size(); 
      //Here we're putting blk-col-size in cptr.
      //Later we will convert cptr to offsets.
      pbuf_cptr_[blk_row_inds[c]] = blkSize/rsize;

      for(Tsize_t n=0; n<blkSize; ++n) {
        pbuf_values1D_[offset++] = pbuf_vals2D[r][c][n];
      }
    }
  }
  pbuf_rptr_[num_block_rows] = pt_r_offset;
  pbuf_bptr_[num_block_rows] = roffset;
  pbuf_indx_[ioffset] = offset;

  //now convert cptr from sizes to offsets;
  LocalOrdinal coffset = 0;
  //we know that pbuf_cptr_.size() is always strictly greater than 0:
  for(Tsize_t c=0; c<pbuf_cptr_.size()-1; ++c) {
    LocalOrdinal csz = pbuf_cptr_[c];
    pbuf_cptr_[c] = coffset;
    coffset += csz; 
  }
  pbuf_cptr_[num_block_cols] = coffset;

  //Final step: release memory for the "2D" storage:
  col_ind_2D_global_ = Teuchos::null;
  col_ind_2D_local_ = Teuchos::null;
  pbuf_values2D_ = Teuchos::null;

  is_storage_optimized_ = true;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillLocalMatrix()
{
  //We insist that optimzeStorage has already been called.
  //We don't care whether this function (fillLocalMatris()) is being
  //called for the first time or not.
  TEST_FOR_EXCEPTION(is_storage_optimized_ != true, std::runtime_error,
    "Tpetra::VbrMatrix::fillLocalMatrix ERROR, optimizeStorage is required to have already been called.");

  lclMatrix_.setPackedValues(pbuf_values1D_,
                             pbuf_rptr_,
                             pbuf_cptr_,
                             pbuf_bptr_,
                             pbuf_bindx_,
                             pbuf_indx_);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillLocalMatVec()
{
  //We insist that optimzeStorage has already been called.
  //We don't care whether this function (fillLocalMatVec()) is being
  //called for the first time or not.
  TEST_FOR_EXCEPTION(is_storage_optimized_ != true, std::runtime_error,
    "Tpetra::VbrMatrix::fillLocalMatrix ERROR, optimizeStorage is required to have already been called.");

  lclMatVec_.initializeValues(lclMatrix_);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap)
{
  //should we return now if is_fill_completed_ is already true ??

  blkDomainMap_ = blockDomainMap;
  blkRangeMap_  = blockRangeMap;

  blkGraph_->fillComplete(convertBlockMapToPointMap(blockDomainMap),
                          convertBlockMapToPointMap(blockRangeMap));

  //Now we need to take the column-map from blkGraph_ and create a
  //corresponding block-column-map.
  //Our block-column-map will have the same distribution as the
  //blkGraph_->getColMap, and we'll get block-size info from the
  //blkDomainMap_. This will require some communication in cases
  //where blkDomainMap is distributed differently.
  blkColMap_ = makeBlockColumnMap(blkDomainMap_, blkGraph_->getColMap());

  optimizeStorage();

  fillLocalMatrix();
  fillLocalMatVec();

  is_fill_completed_ = true;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatVec, class LocalMatSolve>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve>::fillComplete()
{
  //In this case, our block-row-map will also be our domain-map and
  //range-map.

  fillComplete(getBlockRowMap(), getBlockRowMap());
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

