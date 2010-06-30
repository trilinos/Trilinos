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

#include <algorithm>

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
  return Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalElems, blockMap->getNodeBlockIDs(), indexBase, comm, node));
}

//-----------------------------------------------------------------
template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >
makeBlockColumnMap(const Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockMap, const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& ptColMap)
{
  global_size_t numGlobalBlocks = Teuchos::OrdinalTraits<global_size_t>::invalid();
  Teuchos::ArrayView<const GlobalOrdinal> blockIDs = ptColMap->getNodeElementList();
  Teuchos::ArrayView<const LocalOrdinal> blockSizes = blockMap->getNodeBlockSizes();
  Teuchos::Array<GlobalOrdinal> points(blockIDs.size());

  typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type Tsize_t;
  for(Tsize_t i=0; i<blockSizes.size(); ++i) {
    points[i] = blockMap->getFirstGlobalPointInBlock(i);
  }

  //We will create a block-map where each block corresponds to a point in
  //the input-map 'ptColMap', and where each block's size is obtained from
  //the input-block-map 'blockMap'.

  //if we are running on a single processor, then it's easy because
  //we know that blockMap is distributed the same as ptColMap:
  int numProcessors = ptColMap->getComm()->getSize();
  if (numProcessors == 1) {
    return Teuchos::rcp(new Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>(
        numGlobalBlocks, blockIDs, points(), blockSizes,
        ptColMap->getIndexBase(), ptColMap->getComm(), ptColMap->getNode()
    ));
  }

  //If we get to here, we're running on multiple processors, and blockMap
  //is probably not distributed the same as ptColMap, so we have to do
  //some communication to get the block-sizes from blockMap corresponding
  //to the blockIDs we got from ptColMap.
  //We also have to do communication to get the global first-points-in-block
  //for the blocks in the new block-column-map.
  //
  //I think the simplest way to do this is to create vectors where the values
  //are block-sizes (or first-points-in-block), and import one to the other.
  typedef Tpetra::Vector<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Node> GOVec;
  typedef Tpetra::Vector<LocalOrdinal,LocalOrdinal,GlobalOrdinal,Node> LOVec;

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > blkPtMap =
      convertBlockMapToPointMap(blockMap);

  LOVec source_sizes(blkPtMap, blockSizes);
  GOVec source_points(blkPtMap, points());
  LOVec target_sizes(ptColMap);
  GOVec target_points(ptColMap);

  Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> importer(blkPtMap, ptColMap);
  target_sizes.doImport(source_sizes, importer, REPLACE);
  target_points.doImport(source_points, importer, REPLACE);

  blockSizes = target_sizes.get1dView()();

  //now we can create our block-column-map:
  return Teuchos::rcp(new Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>(
    numGlobalBlocks, blockIDs, target_points.get1dView()(), blockSizes,
    ptColMap->getIndexBase(), ptColMap->getComm(), ptColMap->getNode()
  ));
}

template<class Scalar,class Node>
void fill_device_ArrayRCP(Teuchos::RCP<Node>& node, Teuchos::ArrayRCP<Scalar>& ptr, Scalar value)
{
  Kokkos::ReadyBufferHelper<Node> rbh(node);
  Kokkos::InitOp<Scalar> wdp;
  wdp.alpha = value;
  rbh.begin();
  wdp.x = rbh.addNonConstBuffer(ptr);
  rbh.end();
  node->template parallel_for<Kokkos::InitOp<Scalar> >(0, ptr.size(), wdp);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::VbrMatrix(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRowMap, size_t maxNumEntriesPerRow, ProfileType pftype)
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
   importer_(),
   exporter_(),
   importedVec_(),
   exportedVec_(),
   col_ind_2D_global_(Teuchos::rcp(new Teuchos::Array<std::map<GlobalOrdinal,Teuchos::ArrayRCP<Scalar> > >)),
   values2D_(Teuchos::rcp(new Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > >)),
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
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::~VbrMatrix()
{
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getDomainMap() const
{
  return blkDomainMap_->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getRangeMap() const
{
  return blkRangeMap_->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
bool
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::isFillComplete() const
{
  return is_fill_completed_;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
template <class DomainScalar, class RangeScalar>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::multiply(const MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, RangeScalar alpha, RangeScalar beta) const
{
  TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
    "Tpetra::VbrMatrix::multiply ERROR, multiply may only be called after fillComplete has been called.");

  const Kokkos::MultiVector<Scalar,Node> *lclX = &X.getLocalMV();
  Kokkos::MultiVector<Scalar,Node>        *lclY = &Y.getLocalMVNonConst();

  lclMatVec_.template multiply<Scalar,Scalar>(trans,alpha,*lclX,beta,*lclY);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::updateImport(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X) const
{
  if (importer_ != Teuchos::null) {
    if (importedVec_ == Teuchos::null || importedVec_->getNumVectors() != X.getNumVectors()) {
      importedVec_ = Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(blkColMap_->getPointMap(), X.getNumVectors()));
    }

    importedVec_->doImport(X, *importer_, REPLACE);
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::updateExport(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
  if (exporter_ != Teuchos::null) {
    if (exportedVec_ == Teuchos::null || exportedVec_->getNumVectors() != Y.getNumVectors()) {
      exportedVec_ = Teuchos::rcp(new MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(blkColMap_->getPointMap(), Y.getNumVectors()));
    }

    exportedVec_->doExport(Y, *exporter_, REPLACE);
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::apply(
         const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
               MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
               Teuchos::ETransp trans,
               Scalar alpha,
               Scalar beta) const
{
  if (trans == Teuchos::NO_TRANS) {
    updateImport(X);
    const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Xref = importedVec_ == Teuchos::null ? X : *importedVec_;
    this->template multiply<Scalar,Scalar>(Xref, Y, trans, alpha, beta);
  }
  else if (trans == Teuchos::TRANS || trans == Teuchos::CONJ_TRANS) {
    updateImport(Y);
    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Yref = importedVec_ == Teuchos::null ? Y : *importedVec_;
    this->template multiply<Scalar,Scalar>(X, Yref, trans, alpha, beta);
    if (importedVec_ != Teuchos::null) {
      Y.doExport(*importedVec_, *importer_, ADD);
    }
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
bool
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::hasTransposeApply() const
{
  return false;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getBlockRowMap() const
{
  return blkRowMap_;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getPointRowMap() const
{
  return blkRowMap_->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getBlockColMap() const
{
  return blkColMap_;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getPointColMap() const
{
  return blkColMap_->getPointMap();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalBlockEntryViewNonConst(
    GlobalOrdinal globalBlockRow,
    GlobalOrdinal globalBlockCol,
    LocalOrdinal& numPtRows,
    LocalOrdinal& numPtCols,
    Teuchos::ArrayRCP<Scalar>& blockEntry)
{
  //Return a non-const, read-write view of a block-entry (as an ArrayRCP),
  //creating/allocating the block-entry if it doesn't already exist, (but only
  //if fillComplete hasn't been called yet, and if the arguments numPtRows
  //and numPtCols have been set on input).

  if (is_storage_optimized_) {
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockEntryViewNonConst internal ERROR, storage is optimized but isFillComplete() is false.");

    LocalOrdinal localBlockRow = blkRowMap_->getLocalBlockID(globalBlockRow);
    LocalOrdinal localBlockCol = blkColMap_->getLocalBlockID(globalBlockCol);

    getLocalBlockEntryViewNonConst(localBlockRow, localBlockCol,
                                   numPtRows, numPtCols, blockEntry);
    return;
  }

  //If we get to here, fillComplete hasn't been called yet, and the matrix data
  //is stored in un-packed '2D' form.

  if (values2D_->size() == 0) {
    LocalOrdinal numBlockRows = blkGraph_->getNodeNumRows();
    values2D_->resize(numBlockRows);
    col_ind_2D_global_->resize(numBlockRows);
  }

  //this acts as a range-check for globalBlockRow:
  LocalOrdinal localBlockRow = blkRowMap_->getLocalBlockID(globalBlockRow);
  TEST_FOR_EXCEPTION( localBlockRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
     std::runtime_error,
     "Tpetra::VbrMatrix::getGlobalBlockEntryViewNonConst, globalBlockRow not on the local processor.");

  MapGlobalArrayRCP& blkrow = (*col_ind_2D_global_)[localBlockRow];
  typename MapGlobalArrayRCP::iterator col_iter = blkrow.find(globalBlockCol);

  if (col_iter != blkrow.end()) {
    blockEntry = col_iter->second;
  }
  else {
    //blockEntry doesn't already exist, so we will create it.

    //make sure block-size is specified:
    TEST_FOR_EXCEPTION(numPtRows==0 || numPtCols==0, std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockEntryViewNonConst ERROR: creating block-entry, but numPtRows and/or numPtCols is 0.");

    Teuchos::RCP<Node> node = getNode();
    size_t blockSize = numPtRows*numPtCols;
    blockEntry = Teuchos::arcp(new Scalar[blockSize], 0, blockSize);
    std::fill(blockEntry.begin(), blockEntry.end(), 0);
    (*values2D_)[localBlockRow].push_back(blockEntry);
    blkrow.insert(std::make_pair(globalBlockCol, blockEntry));
    blkGraph_->insertGlobalIndices(globalBlockRow, Teuchos::ArrayView<GlobalOrdinal>(&globalBlockCol, 1));
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getGlobalBlockEntryView(
      GlobalOrdinal globalBlockRow,
      GlobalOrdinal globalBlockCol,
      LocalOrdinal& numPtRows,
      LocalOrdinal& numPtCols,
      Teuchos::ArrayRCP<const Scalar>& blockEntry) const
{
  //This method returns a const-view of a block-entry (as an ArrayRCP).
  //Throws an exception if the block-entry doesn't already exist.

  if (is_storage_optimized_) {
    TEST_FOR_EXCEPTION(!isFillComplete(), std::runtime_error,
      "Tpetra::VbrMatrix::getGlobalBlockEntryView internal ERROR, storage is optimized but isFillComplete() is false.");

    LocalOrdinal localBlockRow = blkRowMap_->getLocalBlockID(globalBlockRow);
    LocalOrdinal localBlockCol = blkColMap_->getLocalBlockID(globalBlockCol);
    getLocalBlockEntryView(localBlockRow, localBlockCol, numPtRows, numPtCols, blockEntry);
    return;
  }

  TEST_FOR_EXCEPTION(values2D_->size() == 0, std::runtime_error,
    "Tpetra::VbrMatrix::getGlobalBlockEntryView ERROR, matrix storage not yet allocated, can't return a const view.");

  //this acts as a range-check for globalBlockRow:
  LocalOrdinal localBlockRow = blkRowMap_->getLocalBlockID(globalBlockRow);
  TEST_FOR_EXCEPTION( localBlockRow == Teuchos::OrdinalTraits<LocalOrdinal>::invalid(),
     std::runtime_error,
     "Tpetra::VbrMatrix::getGlobalBlockEntryView, globalBlockRow not on the local processor.");

  MapGlobalArrayRCP& blkrow = (*col_ind_2D_global_)[localBlockRow];
  typename MapGlobalArrayRCP::iterator col_iter = blkrow.find(globalBlockCol);

  if (col_iter == blkrow.end()) {
    throw std::runtime_error("Tpetra::VbrMatrix::getGlobalBlockEntryView ERROR, specified block-entry doesn't exist.");
  }

  numPtRows = blkRowMap_->getBlockSize(localBlockRow);

  TEST_FOR_EXCEPTION(numPtRows == 0, std::runtime_error,
    "Tpetra::VbrMatrix::getGlobalBlockEntryView ERROR, numPtRows == 0.");

  blockEntry = col_iter->second;
  numPtCols = blockEntry.size() / numPtRows;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalBlockEntryViewNonConst(
    LocalOrdinal localBlockRow,
    LocalOrdinal localBlockCol,
    LocalOrdinal& numPtRows,                     
    LocalOrdinal& numPtCols,
    Teuchos::ArrayRCP<Scalar>& blockEntry)
{
  typedef Teuchos::ArrayRCP<const LocalOrdinal> Host_View;
  typedef typename Host_View::iterator ITER;
  //This method returns a non-constant view of a block-entry (as an ArrayRCP).

  TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
   "Tpetra::VbrMatrix::getLocalBlockEntryViewNonConst ERROR, this method can only be called after fillComplete() has been called.");

  TEST_FOR_EXCEPTION(is_storage_optimized_ == false, std::runtime_error,
   "Tpetra::VbrMatrix::getLocalBlockEntryViewNonConst ERROR, this method can only be called if storage is optimized.");

  Teuchos::RCP<Node> node = getNode();

  Host_View bptr = node->template viewBuffer<LocalOrdinal>(2,pbuf_bptr_+localBlockRow);
  LocalOrdinal bindx_offset = bptr[0];
  LocalOrdinal length = bptr[1] - bindx_offset;

  TEST_FOR_EXCEPTION( length < 1, std::runtime_error,
    "Tpetra::VbrMatrix::getLocalBlockEntryViewNonConst ERROR, specified localBlockCol not found in localBlockRow.");

  Host_View bindx = node->template viewBuffer<LocalOrdinal>(length, pbuf_bindx_+bindx_offset);
  ITER bindx_beg = bindx.begin(), bindx_end = bindx.end();
  ITER it = std::lower_bound(bindx_beg, bindx_end, localBlockCol);

  TEST_FOR_EXCEPTION(it == bindx_end || *it != localBlockCol, std::runtime_error,
    "Tpetra::VbrMatrix::getLocalBlockEntryViewNonConst ERROR, specified localBlockCol not found.");

  Host_View rptr = node->template viewBuffer<LocalOrdinal>(2,pbuf_rptr_+localBlockRow),
            cptr = node->template viewBuffer<LocalOrdinal>(2,pbuf_cptr_+localBlockCol);
  numPtRows = rptr[1] - rptr[0];
  numPtCols = cptr[1] - cptr[0];

  const LocalOrdinal blkSize = numPtRows*numPtCols;
  Host_View indx = node->template viewBuffer<LocalOrdinal>(1,pbuf_indx_+bptr[0]+(it-bindx_beg));
  const LocalOrdinal offset = indx[0];
  blockEntry = node->template viewBufferNonConst<Scalar>(Kokkos::ReadWrite, blkSize, pbuf_values1D_ + offset);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getLocalBlockEntryView(
      LocalOrdinal localBlockRow,
      LocalOrdinal localBlockCol,
      LocalOrdinal& numPtRows,
      LocalOrdinal& numPtCols,
      Teuchos::ArrayRCP<const Scalar>& blockEntry) const
{
  typedef Teuchos::ArrayRCP<const LocalOrdinal> Host_View;
  typedef typename Host_View::iterator ITER;
  //This method returns a constant view of a block-entry (as an ArrayRCP).

  TEST_FOR_EXCEPTION(isFillComplete() == false, std::runtime_error,
   "Tpetra::VbrMatrix::getLocalBlockEntryView ERROR, this method can only be called after fillComplete() has been called.");

  TEST_FOR_EXCEPTION(is_storage_optimized_ == false, std::runtime_error,
   "Tpetra::VbrMatrix::getLocalBlockEntryView ERROR, this method can only be called if storage is optimized.");

  Teuchos::RCP<Node> node = getNode();

  Host_View bptr = node->template viewBuffer<LocalOrdinal>(2, pbuf_bptr_+localBlockRow);
  LocalOrdinal bindx_offset = bptr[0];
  LocalOrdinal length = bptr[1] - bindx_offset;

  Host_View bindx = node->template viewBuffer<LocalOrdinal>(length, pbuf_bindx_+bindx_offset);
  ITER bindx_beg = bindx.begin(), bindx_end = bindx.end();
  ITER it = std::lower_bound(bindx_beg, bindx_end, localBlockCol);

  TEST_FOR_EXCEPTION(it == bindx_end || *it != localBlockCol, std::runtime_error,
    "Tpetra::VbrMatrix::getLocalBlockEntryView ERROR, specified localBlockCol not found.");

  Host_View rptr = node->template viewBuffer<LocalOrdinal>(2, pbuf_rptr_+localBlockRow),
            cptr = node->template viewBuffer<LocalOrdinal>(2, pbuf_cptr_+localBlockCol);
  numPtRows = rptr[1] - rptr[0];
  numPtCols = cptr[1] - cptr[0];

  const LocalOrdinal blkSize = numPtRows*numPtCols;
  Host_View indx = node->template viewBuffer<LocalOrdinal>(1, pbuf_indx_+bptr[0]+(it-bindx_beg));
  const LocalOrdinal offset = indx[0];
  blockEntry = node->template viewBuffer<Scalar>(blkSize, pbuf_values1D_ + offset);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
Teuchos::RCP<Node>
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::getNode() const
{
  return blkRowMap_->getPointMap()->getNode();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::putScalar(Scalar s)
{
  Teuchos::RCP<Node> node = getNode();

  if (isFillComplete()) {
    fill_device_ArrayRCP(node, pbuf_values1D_, s);
  }
  else {
    typedef typename Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > >::size_type Tsize_t;
    Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > >& pbuf_vals2D = *values2D_;
    LocalOrdinal numBlockRows = blkGraph_->getNodeNumRows();
    for(LocalOrdinal r=0; r<numBlockRows; ++r) {
      for(Tsize_t c=0; c<pbuf_vals2D[r].size(); ++c) {
        std::fill(pbuf_vals2D[r][c].begin(), pbuf_vals2D[r][c].end(), s);
      }
    }
  }
}

//-------------------------------------------------------------------
template<class ArrayType1,class ArrayType2,class Ordinal>
void set_array_values(ArrayType1& dest, ArrayType2& src, Ordinal rows, Ordinal cols, Ordinal stride, Tpetra::CombineMode mode)
{
  size_t offset = 0;
  size_t src_offset = 0;

  if (mode == ADD) {
    for(Ordinal col=0; col<cols; ++col) {
      for(Ordinal row=0; row<rows; ++row) {
        dest[offset++] += src[src_offset + row];
      }
      src_offset += stride;
    }
  }
  else {
    for(Ordinal col=0; col<cols; ++col) {
      for(Ordinal row=0; row<rows; ++row) {
        dest[offset++] = src[src_offset + row];
      }
      src_offset += stride;
    }
  }
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  LocalOrdinal blkRowSize = blockEntry.numRows();
  LocalOrdinal blkColSize = blockEntry.numCols();
  getGlobalBlockEntryViewNonConst(globalBlockRow,globalBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  //now copy the incoming block-entry into internal storage:
  const Scalar* inputvalues = blockEntry.values();
  set_array_values(internalBlockEntry, inputvalues, blkRowSize, blkColSize, blockEntry.stride(), REPLACE);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, const Teuchos::SerialDenseMatrix<GlobalOrdinal,Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  LocalOrdinal blkRowSize = blockEntry.numRows();
  LocalOrdinal blkColSize = blockEntry.numCols();
  getGlobalBlockEntryViewNonConst(globalBlockRow,globalBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  //now sum (add) the incoming block-entry into internal storage:
  const Scalar* inputvalues = blockEntry.values();
  set_array_values(internalBlockEntry, inputvalues, blkRowSize, blkColSize, blockEntry.stride(), ADD);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::setGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  getGlobalBlockEntryViewNonConst(globalBlockRow,globalBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  LocalOrdinal blk_size = blockEntry.size();
  TEST_FOR_EXCEPTION(blkColSize*LDA > blk_size, std::runtime_error,
    "Tpetra::VbrMatrix::setGlobalBlockEntry ERROR, inconsistent block-entry sizes.");

  //copy the incoming block-entry into internal storage:
  set_array_values(internalBlockEntry, blockEntry, blkRowSize, blkColSize, LDA, REPLACE);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::sumIntoGlobalBlockEntry(GlobalOrdinal globalBlockRow, GlobalOrdinal globalBlockCol, LocalOrdinal blkRowSize, LocalOrdinal blkColSize, LocalOrdinal LDA, const Teuchos::ArrayView<const Scalar>& blockEntry)
{
  //first get an ArrayRCP for the internal storage for this block-entry:
  Teuchos::ArrayRCP<Scalar> internalBlockEntry;
  getGlobalBlockEntryViewNonConst(globalBlockRow,globalBlockCol, blkRowSize, blkColSize, internalBlockEntry);

  LocalOrdinal blk_size = blockEntry.size();
  TEST_FOR_EXCEPTION(blkColSize*LDA > blk_size, std::runtime_error,
    "Tpetra::VbrMatrix::setGlobalBlockEntry ERROR, inconsistent block-entry sizes.");

  //copy the incoming block-entry into internal storage:
  set_array_values(internalBlockEntry, blockEntry, blkRowSize, blkColSize, LDA, ADD);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::optimizeStorage()
{
  if (is_storage_optimized_ == true) return;

  size_t num_block_rows = blkGraph_->getNodeNumRows();
  size_t num_block_cols = blkGraph_->getNodeNumCols();
  size_t num_block_nonzeros = blkGraph_->getNodeNumEntries();

  //need to count the number of point-entries:
  size_t num_point_entries = 0;
  typedef typename Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > >::size_type Tsize_t;
  for(Tsize_t r=0; r<values2D_->size(); ++r) {
    Tsize_t rlen = (*values2D_)[r].size();
    for(Tsize_t c=0; c<rlen; ++c) {
      num_point_entries += (*values2D_)[r][c].size();
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

  Teuchos::ArrayRCP<LocalOrdinal> v_rptr = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, num_block_rows+1, pbuf_rptr_);
  Teuchos::ArrayRCP<LocalOrdinal> v_cptr = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, num_block_cols+1, pbuf_cptr_);
  Teuchos::ArrayRCP<LocalOrdinal> v_bptr = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, num_block_rows+1, pbuf_bptr_);
  Teuchos::ArrayRCP<LocalOrdinal> v_bindx = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, num_block_nonzeros, pbuf_bindx_);
  Teuchos::ArrayRCP<LocalOrdinal> v_indx = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, num_block_nonzeros+1, pbuf_indx_);
  Teuchos::ArrayRCP<Scalar> v_values1D = node->template viewBufferNonConst<Scalar>(Kokkos::WriteOnly, num_point_entries, pbuf_values1D_);

  size_t roffset = 0;
  size_t pt_r_offset = 0;
  size_t ioffset = 0;
  size_t offset = 0;
  Teuchos::Array<Teuchos::Array<Teuchos::ArrayRCP<Scalar> > >& pbuf_vals2D = *values2D_;
  for(Tsize_t r=0; r<pbuf_vals2D.size(); ++r) {
    v_rptr[r] = pt_r_offset;
    LocalOrdinal rsize = blkRowMap_->getBlockSize(r);
    pt_r_offset += rsize;

    v_bptr[r] = roffset;
    Tsize_t rlen = pbuf_vals2D[r].size();
    roffset += rlen;

    Teuchos::ArrayRCP<const LocalOrdinal> blk_row_inds = blkGraph_->getLocalRowView(r);
    MapGlobalArrayRCP& blk_row = (*col_ind_2D_global_)[r];

    for(Tsize_t c=0; c<rlen; ++c) {
      v_bindx[ioffset] = blk_row_inds[c];
      v_indx[ioffset++] = offset;

      Tsize_t blkSize = pbuf_vals2D[r][c].size(); 
      //Here we're putting blk-col-size in cptr.
      //Later we will convert cptr to offsets.
      v_cptr[blk_row_inds[c]] = blkSize/rsize;

      GlobalOrdinal global_col = blkGraph_->getColMap()->getGlobalElement(blk_row_inds[c]);
      typename MapGlobalArrayRCP::iterator iter = blk_row.find(global_col);
      TEST_FOR_EXCEPTION(iter == blk_row.end(), std::runtime_error,
        "Tpetra::VbrMatrix::optimizeStorage ERROR, global_col not found in row.");

      Teuchos::ArrayRCP<Scalar> vals = iter->second;
      for(Tsize_t n=0; n<blkSize; ++n) {
        v_values1D[offset++] = vals[n];
      }
    }
  }
  v_rptr[num_block_rows] = pt_r_offset;
  v_bptr[num_block_rows] = roffset;
  v_indx[ioffset] = offset;

  //now convert cptr from sizes to offsets;
  LocalOrdinal coffset = 0;
  //we know that v_cptr.size() is always strictly greater than 0:
  Tsize_t cptr_len = v_cptr.size()-1;
  for(Tsize_t c=0; c<cptr_len; ++c) {
    LocalOrdinal csz = v_cptr[c];
    v_cptr[c] = coffset;
    coffset += csz; 
  }
  v_cptr[num_block_cols] = coffset;

  //Final step: release memory for the "2D" storage:
  col_ind_2D_global_ = Teuchos::null;
  values2D_ = Teuchos::null;

  is_storage_optimized_ = true;
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillLocalMatrix()
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
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillLocalMatVec()
{
  //We insist that optimzeStorage has already been called.
  //We don't care whether this function (fillLocalMatVec()) is being
  //called for the first time or not.
  TEST_FOR_EXCEPTION(is_storage_optimized_ != true, std::runtime_error,
    "Tpetra::VbrMatrix::fillLocalMatrix ERROR, optimizeStorage is required to have already been called.");

  lclMatVec_.initializeValues(lclMatrix_);
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void
VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockRangeMap)
{
  if (isFillComplete()) return;

  blkDomainMap_ = blockDomainMap;
  blkRangeMap_  = blockRangeMap;

  blkGraph_->fillComplete(convertBlockMapToPointMap(blockDomainMap),
                          convertBlockMapToPointMap(blockRangeMap));

  //Now we need to take the point-column-map from blkGraph_ and create a
  //corresponding block-column-map.
  //Our block-column-map will have the same distribution as
  //blkGraph_->getColMap, and we'll get block-size info from the
  //blkDomainMap_. This will require some communication in cases
  //where blkDomainMap is distributed differently.
  blkColMap_ = makeBlockColumnMap(blkDomainMap_, blkGraph_->getColMap());

  if (!blkDomainMap_->getPointMap()->isSameAs(*blkColMap_->getPointMap())) {
    importer_ = Teuchos::rcp(new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(blkDomainMap_->getPointMap(), blkColMap_->getPointMap()));
  }
  if (!blkRangeMap_->getPointMap()->isSameAs(*blkRowMap_->getPointMap())) {
    exporter_ = Teuchos::rcp(new Tpetra::Export<LocalOrdinal,GlobalOrdinal,Node>(blkRowMap_->getPointMap(), blkRangeMap_->getPointMap()));
  }

  is_fill_completed_ = true;

  optimizeStorage();

  fillLocalMatrix();
  fillLocalMatVec();
}

//-------------------------------------------------------------------
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void VbrMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::fillComplete()
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

