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
#include "Teuchos_Array.hpp"

#ifdef DOXYGEN_USE_ONLY
#include "Tpetra_BlockMap_decl.hpp"
#endif

namespace Tpetra {

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::BlockMap(global_size_t numGlobalBlocks, LocalOrdinal blockSize,
      GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
      const Teuchos::RCP<Node> &node)
 : pointMap_(),
   myGlobalBlockIDs_(),
   blockSizes_(),
   firstPointInBlock_(),
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
  blockSizes_.resize(numLocalBlocks, blockSize);
  firstPointInBlock_.resize(numLocalBlocks);
  LocalOrdinal firstPoint = pointMap_->getMinLocalIndex();
  GlobalOrdinal blockID = pointMap_->getMinGlobalIndex()/blockSize;
  for(size_t i=0; i<numLocalBlocks; ++i) {
    myGlobalBlockIDs_[i] = blockID++;
    firstPointInBlock_[i] = firstPoint;
    firstPoint += blockSize;
  }
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::BlockMap(global_size_t numGlobalBlocks, size_t numLocalBlocks, LocalOrdinal blockSize,
      GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
      const Teuchos::RCP<Node> &node)
 : pointMap_(),
   myGlobalBlockIDs_(),
   blockSizes_(numLocalBlocks, blockSize),
   firstPointInBlock_(),
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

  //we don't need to allocate or fill the blockSizes_ array since blockSize is
  //constant.

  myGlobalBlockIDs_.resize(numLocalBlocks);
  firstPointInBlock_.resize(numLocalBlocks);
  LocalOrdinal firstPoint = pointMap_->getMinLocalIndex();
  GlobalOrdinal blockID = pointMap_->getMinGlobalIndex()/blockSize;
  for(size_t i=0; i<numLocalBlocks; ++i) {
    myGlobalBlockIDs_[i] = blockID++;
    firstPointInBlock_[i] = firstPoint;
    firstPoint += blockSize;
  }
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::BlockMap(global_size_t numGlobalBlocks, const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs, const Teuchos::ArrayView<const LocalOrdinal>& blockSizes, GlobalOrdinal indexBase, const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
      const Teuchos::RCP<Node> &node)
 : pointMap_(),
   myGlobalBlockIDs_(myGlobalBlockIDs),
   blockSizes_(blockSizes),
   firstPointInBlock_(),
   blockIDsAreContiguous_(false),
   constantBlockSize_(0)
{
  TEST_FOR_EXCEPTION(myGlobalBlockIDs_.size()!=blockSizes_.size(), std::runtime_error,
             "Tpetra::BlockMap::BlockMap ERROR: input myGlobalBlockIDs and blockSizes arrays must have the same length.");

  size_t sum_blockSizes = 0;
  typename Teuchos::Array<LocalOrdinal>::const_iterator
    iter = blockSizes_.begin(), iend = blockSizes_.end();
  for(; iter!=iend; ++iter) {
    sum_blockSizes += *iter;
  }

  pointMap_ = Teuchos::rcp(new Map<LocalOrdinal,GlobalOrdinal,Node>(Teuchos::OrdinalTraits<global_size_t>::invalid(), sum_blockSizes, indexBase, comm, node));

  iter = blockSizes_.begin();
  LocalOrdinal firstBlockSize = *iter;
  LocalOrdinal firstPoint = pointMap_->getMinLocalIndex();
  bool blockSizesAreConstant = true;
  for(; iter!=iend; ++iter) {
    firstPointInBlock_.push_back(firstPoint);
    firstPoint += *iter;
    if (*iter != firstBlockSize) {
      blockSizesAreConstant = false;
    }
  }
  if (blockSizesAreConstant) constantBlockSize_ = firstBlockSize;

  size_t num_points = pointMap_->getNodeNumElements();
  TEST_FOR_EXCEPTION(sum_blockSizes != num_points, std::runtime_error,
            "Tpetra::BlockMap::BlockMap ERROR: internal failure, sum of block-sizes must equal pointMap->getNodeNumElements().");

  typename Teuchos::Array<GlobalOrdinal>::const_iterator
    b_iter = myGlobalBlockIDs_.begin(), b_end = myGlobalBlockIDs_.end();
  GlobalOrdinal id = *b_iter;
  ++b_iter;
  for(; b_iter != b_end; ++b_iter) {
    if (*b_iter != id+1) break;
    ++id;
  }
  if (b_iter == b_end) blockIDsAreContiguous_ = true;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::BlockMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& pointMap, const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs, const Teuchos::ArrayView<const LocalOrdinal>& blockSizes)
 : pointMap_(pointMap),
   myGlobalBlockIDs_(myGlobalBlockIDs),
   blockSizes_(blockSizes),
   firstPointInBlock_(),
   blockIDsAreContiguous_(false),
   constantBlockSize_(0)
{
  TEST_FOR_EXCEPTION(myGlobalBlockIDs_.size()!=blockSizes_.size(), std::runtime_error,
             "Tpetra::BlockMap::BlockMap ERROR: input myGlobalBlockIDs and blockSizes arrays must have the same length.");

  LocalOrdinal firstPoint = pointMap->getMinLocalIndex();
  size_t sum_blockSizes = 0;
  typename Teuchos::Array<LocalOrdinal>::const_iterator
    iter = blockSizes_.begin(), iend = blockSizes_.end();
  LocalOrdinal firstBlockSize = *iter;
  bool blockSizesAreConstant = true;
  for(; iter!=iend; ++iter) {
    sum_blockSizes += *iter;
    firstPointInBlock_.push_back(firstPoint);
    firstPoint += *iter;
    if (*iter != firstBlockSize) {
      blockSizesAreConstant = false;
    }
  }
  if (blockSizesAreConstant) constantBlockSize_ = firstBlockSize;

  size_t num_points = pointMap->getNodeNumElements();
  TEST_FOR_EXCEPTION(sum_blockSizes != num_points, std::runtime_error,
            "Tpetra::BlockMap::BlockMap ERROR: sum of block-sizes must equal pointMap->getNodeNumElements().");

  typename Teuchos::Array<GlobalOrdinal>::const_iterator
    b_iter = myGlobalBlockIDs_.begin(), b_end = myGlobalBlockIDs_.end();
  GlobalOrdinal id = *b_iter;
  ++b_iter;
  for(; b_iter != b_end; ++b_iter) {
    if (*b_iter != id+1) break;
    ++id;
  }
  if (b_iter == b_end) blockIDsAreContiguous_ = true;
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
size_t
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumBlocks() const
{
  return myGlobalBlockIDs_.size();
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::ArrayView<const GlobalOrdinal>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getBlockIDs() const
{
  return myGlobalBlockIDs_();
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
bool
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::isBlockSizeConstant() const
{ return constantBlockSize_ != 0; }

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::ArrayView<const LocalOrdinal>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getBlockSizes() const
{
  return blockSizes_();
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::ArrayView<const LocalOrdinal>
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getFirstPointInBlocks() const
{
  return firstPointInBlock_();
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

  if (blockIDsAreContiguous_ != true) {
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
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getBlockSize(LocalOrdinal localBlockID) const
{
  if (constantBlockSize_ != 0) {
    return constantBlockSize_;
  }

  //should this be a debug-mode-only range check?
  if (localBlockID < 0 || localBlockID >= blockSizes_.size()) {
    throw std::runtime_error("Tpetra::BlockMap::getBlockSize ERROR: localBlockID out of range.");
  }

  return blockSizes_[localBlockID];
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
LocalOrdinal
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getFirstLocalPointInBlock(LocalOrdinal localBlockID) const
{
  //should this be a debug-mode-only range check?
  if (localBlockID < 0 || localBlockID >= firstPointInBlock_.size()) {
    throw std::runtime_error("Tpetra::BlockMap::getFirstLocalPointInBlock ERROR: localBlockID out of range.");
  }

  return firstPointInBlock_[localBlockID];
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
GlobalOrdinal
BlockMap<LocalOrdinal,GlobalOrdinal,Node>::getFirstGlobalPointInBlock(LocalOrdinal localBlockID) const
{
  //should this be a debug-mode-only range check?
  if (localBlockID < 0 || localBlockID >= firstPointInBlock_.size()) {
    throw std::runtime_error("Tpetra::BlockMap::getFirstGlobalPointInBlock ERROR: localBlockID out of range.");
  }

  return pointMap_->getGlobalElement(firstPointInBlock_[localBlockID]);
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_BLOCKMAP_INSTANT(LO,GO,NODE) \
  \
  template class BlockMap< LO , GO , NODE >;


}//namespace Tpetra

#endif

