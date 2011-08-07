// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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
// ***********************************************************************
// @HEADER

#ifndef TPETRA_BLOCKMAP_DEF_HPP
#define TPETRA_BLOCKMAP_DEF_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Distributor.hpp"

#ifdef DOXYGEN_USE_ONLY
#include "Tpetra_BlockMap_decl.hpp"
#endif

namespace Tpetra {

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::BlockMap(
    global_size_t numGlobalBlocks,
    LocalOrdinal blockSize,
    GlobalOrdinal indexBase,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
    const Teuchos::RCP<Node> &node)
 : pointMap_(),
   globalNumBlocks_(numGlobalBlocks),
   myGlobalBlockIDs_(),
   pbuf_firstPointInBlock_(),
   view_firstPointInBlock_(),
   blockIDsAreContiguous_(true),
   constantBlockSize_(blockSize)
{
  TEST_FOR_EXCEPTION( blockSize <= 0, std::runtime_error,
       "Tpetra::BlockMap::BlockMap ERROR: blockSize must be greater than 0.");

  global_size_t numGlobalPoints = numGlobalBlocks*blockSize;

  //we compute numLocalPoints to make sure that Tpetra::Map doesn't split
  //numGlobalPoints in a way that would separate the points within a block
  //onto different mpi processors.

  size_t numLocalPoints = numGlobalBlocks/comm->getSize();
  int remainder = numGlobalBlocks%comm->getSize();
  int localProc = comm->getRank();
  if (localProc < remainder) ++numLocalPoints;
  numLocalPoints *= blockSize;

  //now create the point-map specifying both numGlobalPoints and numLocalPoints:
  pointMap_ = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalPoints, numLocalPoints, indexBase, comm, node));

  size_t numLocalBlocks = pointMap_->getNodeNumElements()/blockSize;
  size_t checkLocalBlocks = numLocalPoints/blockSize;
  //can there be an inconsistency here???
  TEST_FOR_EXCEPTION(numLocalBlocks != checkLocalBlocks, std::runtime_error,
       "Tpetra::BlockMap::BlockMap ERROR: internal failure, numLocalBlocks not consistent with point-map.");
  
  myGlobalBlockIDs_.resize(numLocalBlocks);
  pbuf_firstPointInBlock_ = node->template allocBuffer<LocalOrdinal>(numLocalBlocks+1);
  Teuchos::ArrayRCP<LocalOrdinal> v_firstPoints = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, numLocalBlocks+1, pbuf_firstPointInBlock_);

  LocalOrdinal firstPoint = pointMap_->getMinLocalIndex();
  GlobalOrdinal blockID = pointMap_->getMinGlobalIndex()/blockSize;
  for(size_t i=0; i<numLocalBlocks; ++i) {
    myGlobalBlockIDs_[i] = blockID++;
    v_firstPoints[i] = firstPoint;
    firstPoint += blockSize;
  }
  v_firstPoints[numLocalBlocks] = firstPoint;
  v_firstPoints = Teuchos::null;
  view_firstPointInBlock_ = node->template viewBuffer<LocalOrdinal>(numLocalBlocks+1, pbuf_firstPointInBlock_);
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::BlockMap(
    global_size_t numGlobalBlocks,
    size_t numLocalBlocks,
    LocalOrdinal blockSize,
    GlobalOrdinal indexBase,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
    const Teuchos::RCP<Node> &node)
 : pointMap_(),
   globalNumBlocks_(numGlobalBlocks),
   myGlobalBlockIDs_(),
   pbuf_firstPointInBlock_(),
   view_firstPointInBlock_(),
   blockIDsAreContiguous_(true),
   constantBlockSize_(blockSize)
{
  TEST_FOR_EXCEPTION( blockSize <= 0, std::runtime_error,
       "Tpetra::BlockMap::BlockMap ERROR: blockSize must be greater than 0.");

  global_size_t numGlobalPoints = numGlobalBlocks*blockSize;
  if (numGlobalBlocks == Teuchos::OrdinalTraits<global_size_t>::invalid()) {
    numGlobalPoints = Teuchos::OrdinalTraits<global_size_t>::invalid();
  }

  size_t numLocalPoints = numLocalBlocks*blockSize;

  //now create the point-map specifying both numGlobalPoints and numLocalPoints:
  pointMap_ = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(numGlobalPoints, numLocalPoints, indexBase, comm, node));

  if (numGlobalBlocks == Teuchos::OrdinalTraits<global_size_t>::invalid()) {
    Teuchos::reduceAll<int, global_size_t>(*comm, Teuchos::REDUCE_SUM,
             1, &numLocalBlocks, &globalNumBlocks_);
  }

  myGlobalBlockIDs_.resize(numLocalBlocks);
  pbuf_firstPointInBlock_ = node->template allocBuffer<LocalOrdinal>(numLocalBlocks+1);
  Teuchos::ArrayRCP<LocalOrdinal> v_firstPoints = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, numLocalBlocks+1, pbuf_firstPointInBlock_);

  LocalOrdinal firstPoint = pointMap_->getMinLocalIndex();
  GlobalOrdinal blockID = pointMap_->getMinGlobalIndex()/blockSize;
  for(size_t i=0; i<numLocalBlocks; ++i) {
    myGlobalBlockIDs_[i] = blockID++;
    v_firstPoints[i] = firstPoint;
    firstPoint += blockSize;
  }
  v_firstPoints[numLocalBlocks] = firstPoint;
  v_firstPoints = Teuchos::null;
  view_firstPointInBlock_ = node->template viewBuffer<LocalOrdinal>(numLocalBlocks+1, pbuf_firstPointInBlock_);
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::BlockMap(
    global_size_t numGlobalBlocks,
    const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs,
    const Teuchos::ArrayView<const GlobalOrdinal>& myFirstGlobalPointInBlocks,
    const Teuchos::ArrayView<const LocalOrdinal>& blockSizes,
    GlobalOrdinal indexBase,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
    const Teuchos::RCP<Node> &node)
 : pointMap_(),
   globalNumBlocks_(numGlobalBlocks),
   myGlobalBlockIDs_(myGlobalBlockIDs),
   pbuf_firstPointInBlock_(),
   view_firstPointInBlock_(),
   blockIDsAreContiguous_(false),
   constantBlockSize_(0)
{
  TEST_FOR_EXCEPTION(myGlobalBlockIDs_.size()!=blockSizes.size(), std::runtime_error,
             "Tpetra::BlockMap::BlockMap ERROR: input myGlobalBlockIDs and blockSizes arrays must have the same length.");

  size_t sum_blockSizes = 0;
  Teuchos::Array<GlobalOrdinal> myGlobalPoints;
  typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator
    iter = blockSizes.begin(), iend = blockSizes.end();
  size_t i = 0;
  for(; iter!=iend; ++iter) {
    LocalOrdinal bsize = *iter;
    sum_blockSizes += bsize;
    GlobalOrdinal firstPoint = myFirstGlobalPointInBlocks[i++];
    for(LocalOrdinal j=0; j<bsize; ++j) {
      myGlobalPoints.push_back(firstPoint+j);
    }
  }

  pointMap_ = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(Teuchos::OrdinalTraits<global_size_t>::invalid(), myGlobalPoints(), indexBase, comm, node));

  global_size_t global_sum;
  Teuchos::reduceAll<int,global_size_t>(*comm, Teuchos::REDUCE_SUM,
                Teuchos::as<global_size_t>(myGlobalBlockIDs.size()), Teuchos::outArg(global_sum));
  globalNumBlocks_ = global_sum;

  iter = blockSizes.begin();
  LocalOrdinal firstBlockSize = iter==iend ? 0 : *iter;
  LocalOrdinal firstPoint = pointMap_->getMinLocalIndex();
  LocalOrdinal numLocalBlocks = myGlobalBlockIDs.size();
  pbuf_firstPointInBlock_ = node->template allocBuffer<LocalOrdinal>(numLocalBlocks+1);
  Teuchos::ArrayRCP<LocalOrdinal> v_firstPoints = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, numLocalBlocks+1, pbuf_firstPointInBlock_);

  bool blockSizesAreConstant = true;
  i=0;
  for(; iter!=iend; ++iter, ++i) {
    v_firstPoints[i] = firstPoint;
    firstPoint += *iter;
    if (*iter != firstBlockSize) {
      blockSizesAreConstant = false;
    }
  }
  v_firstPoints[i] = firstPoint;
  v_firstPoints = Teuchos::null;
  if (blockSizesAreConstant) constantBlockSize_ = firstBlockSize;

  size_t num_points = pointMap_->getNodeNumElements();
  TEST_FOR_EXCEPTION(sum_blockSizes != num_points, std::runtime_error,
            "Tpetra::BlockMap::BlockMap ERROR: internal failure, sum of block-sizes must equal pointMap->getNodeNumElements().");

  typename Teuchos::Array<GlobalOrdinal>::const_iterator
    b_iter = myGlobalBlockIDs_.begin(), b_end = myGlobalBlockIDs_.end();
  GlobalOrdinal id = b_iter==b_end ? 0 : *b_iter;
  if (b_iter != b_end) ++b_iter;
  blockIDsAreContiguous_ = true;
  for(; b_iter != b_end; ++b_iter) {
    if (*b_iter != id+1) {
      blockIDsAreContiguous_ = false;
      break;
    }
    ++id;
  }
  if (blockIDsAreContiguous_ == false) {
    setup_noncontig_mapping();
  }

  view_firstPointInBlock_ = node->template viewBuffer<LocalOrdinal>(numLocalBlocks+1, pbuf_firstPointInBlock_);
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::BlockMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& pointMap, const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs, const Teuchos::ArrayView<const LocalOrdinal>& blockSizes, const Teuchos::RCP<Node>& node)
 : pointMap_(pointMap),
   globalNumBlocks_(Teuchos::OrdinalTraits<LocalOrdinal>::invalid()),
   myGlobalBlockIDs_(myGlobalBlockIDs),
   pbuf_firstPointInBlock_(),
   view_firstPointInBlock_(),
   blockIDsAreContiguous_(false),
   constantBlockSize_(0)
{
  TEST_FOR_EXCEPTION(myGlobalBlockIDs_.size()!=blockSizes.size(), std::runtime_error,
             "Tpetra::BlockMap::BlockMap ERROR: input myGlobalBlockIDs and blockSizes arrays must have the same length.");

  global_size_t numLocalBlocks = myGlobalBlockIDs.size();
  Teuchos::reduceAll<int,global_size_t>(*pointMap->getComm(), Teuchos::REDUCE_SUM, 1, &numLocalBlocks, &globalNumBlocks_);

  LocalOrdinal firstPoint = pointMap->getMinLocalIndex();
  size_t sum_blockSizes = 0;
  typename Teuchos::ArrayView<const LocalOrdinal>::const_iterator
    iter = blockSizes.begin(), iend = blockSizes.end();
  LocalOrdinal firstBlockSize = *iter;

  pbuf_firstPointInBlock_ = node->template allocBuffer<LocalOrdinal>(myGlobalBlockIDs.size()+1);
  Teuchos::ArrayRCP<LocalOrdinal> v_firstPoints = node->template viewBufferNonConst<LocalOrdinal>(Kokkos::WriteOnly, numLocalBlocks+1, pbuf_firstPointInBlock_);

  bool blockSizesAreConstant = true;
  size_t i=0;
  for(; iter!=iend; ++iter, ++i) {
    sum_blockSizes += *iter;
    v_firstPoints[i] = firstPoint;
    firstPoint += *iter;
    if (*iter != firstBlockSize) {
      blockSizesAreConstant = false;
    }
  }
  v_firstPoints[i] = firstPoint;
  v_firstPoints = Teuchos::null;
  if (blockSizesAreConstant) constantBlockSize_ = firstBlockSize;

  size_t num_points = pointMap->getNodeNumElements();
  TEST_FOR_EXCEPTION(sum_blockSizes != num_points, std::runtime_error,
            "Tpetra::BlockMap::BlockMap ERROR: sum of block-sizes must equal pointMap->getNodeNumElements().");

  typename Teuchos::Array<GlobalOrdinal>::const_iterator
    b_iter = myGlobalBlockIDs_.begin(), b_end = myGlobalBlockIDs_.end();
  GlobalOrdinal id = *b_iter;
  ++b_iter;
  blockIDsAreContiguous_ = true;
  for(; b_iter != b_end; ++b_iter) {
    if (*b_iter != id+1) {
      blockIDsAreContiguous_ = false;
      break;
    }
    ++id;
  }
  if (blockIDsAreContiguous_ == false) {
    setup_noncontig_mapping();
  }

  view_firstPointInBlock_ = node->template viewBuffer<LocalOrdinal>(numLocalBlocks+1, pbuf_firstPointInBlock_);
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
void
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getRemoteBlockInfo(
    const Teuchos::ArrayView<const GlobalOrdinal>& GBIDs,
    const Teuchos::ArrayView<GlobalOrdinal>& firstGlobalPointInBlocks,
    const Teuchos::ArrayView<LocalOrdinal>& blockSizes) const
{
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > pbmap =
      convertBlockMapToPointMap(*this);

  Teuchos::Array<int> remoteProcs(GBIDs.size());
  pbmap->getRemoteIndexList(GBIDs(), remoteProcs());

  Tpetra::Distributor distor(pbmap->getComm());

  Tpetra::ArrayRCP<GlobalOrdinal> exportGBIDs;
  Tpetra::ArrayRCP<int> exportProcs;

  distor.createFromRecvs(GBIDs(), remoteProcs(), exportGBIDs, exportProcs);

  Tpetra::Array<GlobalOrdinal> exportFirstPointInBlocks(exportGBIDs.size());
  Tpetra::Array<LocalOrdinal> exportSizes(exportGBIDs.size());

  typename Teuchos::ArrayRCP<GlobalOrdinal>::const_iterator
    iter = exportGBIDs.begin(),
    iter_end = exportGBIDs.end();
  for(int i=0; iter!=iter_end; ++iter, ++i) {
    GlobalOrdinal gbid = *iter;
    LocalOrdinal lbid = getLocalBlockID(gbid);
    exportFirstPointInBlocks[i] = getFirstGlobalPointInLocalBlock(lbid);
    exportSizes[i] = getLocalBlockSize(lbid);
  }

  size_t numPackets = 1;
  distor.doPostsAndWaits(exportFirstPointInBlocks().getConst(),
                         numPackets, firstGlobalPointInBlocks);

  distor.doPostsAndWaits(exportSizes().getConst(), numPackets, blockSizes);
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
void
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::setup_noncontig_mapping()
{
  //possibly need to use a hash (unordered_map) here...
  typedef typename Teuchos::Array<GlobalOrdinal>::size_type Tsize_t;
  for(Tsize_t i=0; i<myGlobalBlockIDs_.size(); ++i) {
    LocalOrdinal li = i;
    map_global_to_local_.insert(std::make_pair(myGlobalBlockIDs_[i],li));
  }
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
global_size_t
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumBlocks() const
{
  return globalNumBlocks_;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
size_t
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumBlocks() const
{
  return myGlobalBlockIDs_.size();
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::ArrayView<const GlobalOrdinal>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getNodeBlockIDs() const
{
  return myGlobalBlockIDs_();
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
bool
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::isBlockSizeConstant() const
{ return constantBlockSize_ != 0; }

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::ArrayRCP<const LocalOrdinal>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getNodeFirstPointInBlocks() const
{
  return view_firstPointInBlock_;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::ArrayRCP<const LocalOrdinal>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getNodeFirstPointInBlocks_Device() const
{
  return pbuf_firstPointInBlock_;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
GlobalOrdinal
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getGlobalBlockID(LocalOrdinal localBlockID) const
{
  LocalOrdinal invalid = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  if (localBlockID < 0 || localBlockID >= myGlobalBlockIDs_.size()) {
    return invalid;
  }

  return myGlobalBlockIDs_[localBlockID];
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
LocalOrdinal
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getLocalBlockID(GlobalOrdinal globalBlockID) const
{
  LocalOrdinal invalid = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
  if (myGlobalBlockIDs_.size() == 0) {
    return invalid;
  }

  if (blockIDsAreContiguous_ == false) {
    //need to use a hash (unordered_map) here instead of a map...
    typename std::map<GlobalOrdinal,LocalOrdinal>::const_iterator iter =
      map_global_to_local_.find(globalBlockID);
    if (iter == map_global_to_local_.end()) return invalid;
    return iter->second;
  }

  LocalOrdinal localBlockID = globalBlockID - myGlobalBlockIDs_[0];
  LocalOrdinal numLocalBlocks = myGlobalBlockIDs_.size();

  if (localBlockID < 0 || localBlockID >= numLocalBlocks) {
    return invalid;
  }

  return localBlockID;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
LocalOrdinal
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getLocalBlockSize(LocalOrdinal localBlockID) const
{
  if (constantBlockSize_ != 0) {
    return constantBlockSize_;
  }

  //should this be a debug-mode-only range check?
  if (localBlockID < 0 || localBlockID >= view_firstPointInBlock_.size()) {
    throw std::runtime_error("Tpetra::BlockMap::getLocalBlockSize ERROR: localBlockID out of range.");
  }

  return view_firstPointInBlock_[localBlockID+1]-view_firstPointInBlock_[localBlockID];
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
LocalOrdinal
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getFirstLocalPointInLocalBlock(LocalOrdinal localBlockID) const
{
  //should this be a debug-mode-only range check?
  if (localBlockID < 0 || localBlockID >= view_firstPointInBlock_.size()) {
    throw std::runtime_error("Tpetra::BlockMap::getFirstLocalPointInLocalBlock ERROR: localBlockID out of range.");
  }

  return view_firstPointInBlock_[localBlockID];
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
GlobalOrdinal
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getFirstGlobalPointInLocalBlock(LocalOrdinal localBlockID) const
{
  //should this be a debug-mode-only range check?
  if (localBlockID < 0 || localBlockID >= view_firstPointInBlock_.size()) {
    throw std::runtime_error("Tpetra::BlockMap::getFirstGlobalPointInLocalBlock ERROR: localBlockID out of range.");
  }

  return pointMap_->getGlobalElement(view_firstPointInBlock_[localBlockID]);
}

}//namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_BLOCKMAP_INSTANT(LO,GO,NODE) \
  \
  template class BlockMap< LO , GO , NODE >;

#endif
