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

#ifndef TPETRA_BLOCKCRSGRAPH_DEF_HPP
#define TPETRA_BLOCKCRSGRAPH_DEF_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Vector.hpp"

#ifdef DOXYGEN_USE_ONLY
#include "Tpetra_BlockCrsGraph_decl.hpp"
#endif

namespace Tpetra {

//-----------------------------------------------------------------
template<class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >
makeBlockColumnMap(
    const Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockMap,
    const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >& ptColMap)
{
  global_size_t numGlobalBlocks = Teuchos::OrdinalTraits<global_size_t>::invalid();
  Teuchos::ArrayView<const GlobalOrdinal> blockIDs = ptColMap->getNodeElementList();
  Teuchos::ArrayRCP<const LocalOrdinal> firstPoints = blockMap->getNodeFirstPointInBlocks();
  Teuchos::Array<GlobalOrdinal> points(firstPoints.size()-1);
  Teuchos::Array<LocalOrdinal> blockSizes(firstPoints.size()-1);

  typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type Tsize_t;
  for(Tsize_t i=0; i<blockSizes.size(); ++i) {
    points[i] = blockMap->getFirstGlobalPointInLocalBlock(i);
    blockSizes[i] = firstPoints[i+1]-firstPoints[i];
  }

  //We will create a block-map where each block corresponds to a point in
  //the input-map 'ptColMap', and where each block's size is obtained from
  //the input-block-map 'blockMap'.

  //if we are running on a single processor, then it's easy because
  //we know that blockMap is distributed the same as ptColMap:
  int numProcessors = ptColMap->getComm()->getSize();
  if (numProcessors == 1) {
    return Teuchos::rcp(new Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>(
        numGlobalBlocks, blockIDs, points(), blockSizes(),
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
      convertBlockMapToPointMap(*blockMap);

  LOVec source_sizes(blkPtMap, blockSizes());
  GOVec source_points(blkPtMap, points());
  LOVec target_sizes(ptColMap);
  GOVec target_points(ptColMap);

  Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> importer(blkPtMap, ptColMap);
  target_sizes.doImport(source_sizes, importer, REPLACE);
  target_points.doImport(source_points, importer, REPLACE);

  //now we can create our block-column-map:
  return Teuchos::rcp(new Tpetra::BlockMap<LocalOrdinal,GlobalOrdinal,Node>(
    numGlobalBlocks, blockIDs, target_points.get1dView()(), target_sizes.get1dView()(),
    ptColMap->getIndexBase(), ptColMap->getComm(), ptColMap->getNode()
  ));
}

template<class LocalOrdinal,class GlobalOrdinal,class Node>
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::BlockCrsGraph(
   const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blkRowMap,
   size_t maxNumEntriesPerRow, ProfileType pftype
    )
 : ptGraph_(),
   blkRowMap_(blkRowMap),
   blkColMap_(),
   blkDomainMap_(),
   blkRangeMap_()
{
  ptGraph_ = Teuchos::rcp(new CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(
          convertBlockMapToPointMap(*blkRowMap), maxNumEntriesPerRow, pftype));
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::insertGlobalIndices(GlobalOrdinal row, const Teuchos::ArrayView<const GlobalOrdinal> &indices)
{
  ptGraph_->insertGlobalIndices(row, indices);
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const size_t>
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeRowOffsets() const
{
  return ptGraph_->getNodeRowPtrs();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::ArrayRCP<const LocalOrdinal>
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodePackedIndices() const
{
  return ptGraph_->getNodePackedIndices();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::globalAssemble()
{
  ptGraph_->globalAssemble();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::fillComplete(const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkDomainMap, const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &blkRangeMap, OptimizeOption os)
{
  blkDomainMap_ = blkDomainMap;
  blkRangeMap_  = blkRangeMap;

  RCP<ParameterList> params = parameterList();
  if (os == DoOptimizeStorage) params->set("Optimize Storage",true);
  else                         params->set("Optimize Storage",false);
  ptGraph_->fillComplete(convertBlockMapToPointMap(*blkDomainMap),
                         convertBlockMapToPointMap(*blkRangeMap), params);

  //Now we need to take the point-column-map from ptGraph_ and create a
  //corresponding block-column-map.
  //Our block-column-map will have the same distribution as
  //blkGraph_->getColMap, and we'll get block-size info from the
  //blkDomainMap_. This will require some communication in cases
  //where blkDomainMap is distributed differently.
  blkColMap_ = makeBlockColumnMap(blkDomainMap_, ptGraph_->getColMap());
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::fillComplete(OptimizeOption os)
{
  fillComplete(getBlockRowMap(), getBlockRowMap(), os);
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::optimizeStorage()
{
  if (ptGraph_->isFillComplete()) {
    ptGraph_->resumeFill();
  }
  RCP<ParameterList> params = parameterList();
  params->set("Optimize Storage",true);
  ptGraph_->fillComplete(params);
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isFillComplete() const
{
  return ptGraph_->isFillComplete();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isUpperTriangular() const
{
  return ptGraph_->isUpperTriangular();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isLowerTriangular() const
{
  return ptGraph_->isLowerTriangular();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
bool
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::isLocallyIndexed() const
{
  return ptGraph_->isLocallyIndexed();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumBlockEntries() const
{
  return ptGraph_->getNodeNumEntries();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalBlockRowLength(GlobalOrdinal row) const
{
  return ptGraph_->getNumEntriesInGlobalRow(row);
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalBlockRowView(GlobalOrdinal row,
                             Teuchos::ArrayView<const GlobalOrdinal>& blockCols) const
{
  ptGraph_->getGlobalRowView(row, blockCols);
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
void
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getLocalBlockRowView(LocalOrdinal row,
                             Teuchos::ArrayView<const LocalOrdinal>& blockCols) const
{
  ptGraph_->getLocalRowView(row, blockCols);
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumBlockRows() const
{
  return ptGraph_->getNodeNumRows();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getGlobalNumBlockRows() const
{
  return ptGraph_->getGlobalNumRows();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
size_t
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getNodeNumBlockDiags() const
{
  return ptGraph_->getNodeNumDiags();
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getBlockRowMap() const
{
  return blkRowMap_;
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getBlockColMap() const
{
  return blkColMap_;
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getBlockDomainMap() const
{
  return blkDomainMap_;
}

//-------------------------------------------------------------------
template<class LocalOrdinal, class GlobalOrdinal, class Node>
const Teuchos::RCP<const BlockMap<LocalOrdinal,GlobalOrdinal,Node> > &
BlockCrsGraph<LocalOrdinal,GlobalOrdinal,Node>::getBlockRangeMap() const
{
  return blkRangeMap_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_BLOCKCRSGRAPH_INSTANT(LO,GO,NODE) \
  \
  template class BlockCrsGraph< LO , GO , NODE >;


}//namespace Tpetra

#endif

