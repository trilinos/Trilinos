// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "CartesianConnManager.hpp"

#include "Panzer_FieldPattern.hpp"
#include "Shards_CellTopology.hpp"
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {
namespace unit_test {

void CartesianConnManager::
initialize(const Teuchos::MpiComm<int> & comm,GlobalOrdinal nx, GlobalOrdinal ny,
                                              int px, int py,
                                              int bx, int by,
                                              const shards::CellTopology elemTopo)
{
  // preconditions
  TEUCHOS_ASSERT(0<nx*ny);
  TEUCHOS_TEST_FOR_EXCEPTION(comm.getSize()!=px*py,std::logic_error,
                                "Processor assignment does not equal processor count: " << comm.getSize()
                             << " != " << px << " * " << py );
  TEUCHOS_ASSERT(0<bx*by);
  TEUCHOS_ASSERT(nx/px>=1.0);
  TEUCHOS_ASSERT(ny/py>=1.0);

  numProc_ = comm.getSize();
  myRank_ = comm.getRank();

  dim_ = 2;

  brickElements_.x = nx;
  brickElements_.y = ny;
  brickElements_.z = 1;

  processors_.x = px;
  processors_.y = py;
  processors_.z = 1;  // even though its 2D do this

  blocks_.x = bx;
  blocks_.y = by;
  blocks_.z = 1;  // even though its 2D do this

  totalBrickElements_.x = nx*bx;
  totalBrickElements_.y = ny*by;
  totalBrickElements_.z = 1*1;

  totalNodes_ = (bx*nx+1)*(by*ny+1);
  totalEdges_ = (bx*nx+1)*(by*ny)+(bx*nx)*(by*ny+1);
  totalFaces_ = (bx*nx)*(by*ny);
  totalElements_ = totalFaces_;

  switch (elemTopo.getBaseKey()) {
  case shards::Quadrilateral<4>::key:
  {
    numSubElemsPerBrickElement_ =1;
    subElemToBrickElementNodesMap_.resize(numSubElemsPerBrickElement_);
    subElemToBrickElementNodesMap_[0].resize(4);
    for (int i=0; i<4; ++i) {
      subElemToBrickElementNodesMap_[0][i] = i;
    }

    subElemToBrickElementEdgesMap_.resize(numSubElemsPerBrickElement_);
    subElemToBrickElementEdgesMap_[0].resize(4);
    for (int i=0; i<4; ++i)
      subElemToBrickElementEdgesMap_[0][i] = i;
  } break;
  case shards::Triangle<3>::key: //split Quad into 2 tets
  {
    totalEdges_ = totalEdges_ + totalFaces_;
    totalFaces_ *= 2;
    totalElements_ *= totalFaces_;
    numSubElemsPerBrickElement_ = 2;

    int trianglesInQuad[2][3] = {
        {0, 1, 3},
        {1, 2, 3}
    };

    // edges from 0 to 3 are the Quad edges,
    // edge 4 is a diagonal edge internal to the Quad
    int triEdgesInQuad[2][3] = {
        {0, 4, 3},
        {1, 2, 4}
    };

    subElemToBrickElementNodesMap_.resize(numSubElemsPerBrickElement_);
    for (int i=0; i<numSubElemsPerBrickElement_; ++i) {
      subElemToBrickElementNodesMap_[i].resize(3);
      for (int j=0; j<3; ++j)
        subElemToBrickElementNodesMap_[i][j]=trianglesInQuad[i][j];
    }

    subElemToBrickElementEdgesMap_.resize(numSubElemsPerBrickElement_);
    for (int i=0; i<numSubElemsPerBrickElement_; ++i) {
      subElemToBrickElementEdgesMap_[i].resize(3);
      for (int j=0; j<3; ++j)
        subElemToBrickElementEdgesMap_[i][j] = triEdgesInQuad[i][j];
    }
  } break;
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true,
            std::logic_error,
            "CartesianConnManager: element not supported" << std::endl);
  }

  elemTopology_ = elemTopo;

  buildLocalElements();
}

void CartesianConnManager::
initialize(const Teuchos::MpiComm<int> & comm,GlobalOrdinal nx, GlobalOrdinal ny, GlobalOrdinal nz,
                                              int px, int py,int pz,
                                              int bx, int by,int bz,
                                              const shards::CellTopology elemTopo)
{
  // preconditions
  TEUCHOS_ASSERT(0<nx*ny*nz);
  TEUCHOS_TEST_FOR_EXCEPTION(comm.getSize()!=px*py*pz,std::logic_error,
                                "Processor assignment does not equal processor count: " << comm.getSize()
                             << " != " << px << " * " << py << " * " << pz);
  TEUCHOS_ASSERT(0<bx*by*bz);
  TEUCHOS_ASSERT(nx/px>=1.0);
  TEUCHOS_ASSERT(ny/py>=1.0);
  TEUCHOS_ASSERT(nz/pz>=1.0);

  numProc_ = comm.getSize();
  myRank_ = comm.getRank();

  dim_ = 3;

  brickElements_.x = nx;
  brickElements_.y = ny;
  brickElements_.z = nz;

  processors_.x = px;
  processors_.y = py;
  processors_.z = pz;

  blocks_.x = bx;
  blocks_.y = by;
  blocks_.z = bz;

  totalBrickElements_.x = nx*bx;
  totalBrickElements_.y = ny*by;
  totalBrickElements_.z = nz*bz;

  totalNodes_ = (bx*nx+1)*(by*ny+1)*(bz*nz+1);
  totalEdges_ = (bx*nx+1)*(by*ny)*(bz*nz+1)+(bx*nx)*(by*ny+1)*(bz*nz+1)+(bx*nx+1)*(by*ny+1)*(bz*nz);
  totalFaces_ = (bx*nx)*(by*ny)*(bz*nz+1) + (bx*nx)*(by*ny+1)*(bz*nz) + (bx*nx+1)*(by*ny)*(bz*nz);
  totalElements_ = (bx*nx)*(by*ny)*(bz*nz);

  switch (elemTopo.getBaseKey()) {
  case shards::Hexahedron<8>::key:
  {
    numSubElemsPerBrickElement_ =1;
    subElemToBrickElementNodesMap_.resize(numSubElemsPerBrickElement_);
    subElemToBrickElementNodesMap_[0].resize(8);
    for (int i=0; i<8; ++i) {
      subElemToBrickElementNodesMap_[0][i] = i;
    }
    subElemToBrickElementFacesMap_.resize(numSubElemsPerBrickElement_);
    subElemToBrickElementFacesMap_[0].resize(6);
    for (int i=0; i<6; ++i)
      subElemToBrickElementFacesMap_[0][i] = i;
    subElemToBrickElementEdgesMap_.resize(numSubElemsPerBrickElement_);
    subElemToBrickElementEdgesMap_[0].resize(12);
    for (int i=0; i<12; ++i)
      subElemToBrickElementEdgesMap_[0][i] = i;
  } break;
  case shards::Tetrahedron<4>::key: //split Hex into 6 tets
  {
    totalEdges_ = totalEdges_ + totalFaces_ + totalElements_;
    totalFaces_ = 2*totalFaces_ + 6*totalElements_;
    totalElements_ *= 6;
    numSubElemsPerBrickElement_ = 6;

    int tetsInHex[6][4] = {
        {0,2,3,6},
        {0,3,7,6},
        {0,7,4,6},
        {0,5,6,4},
        {1,5,6,0},
        {1,6,2,0}
    };

    // edges from 0 to 11 are the Hex edges,
    // edges from 12 to 17, are diagonal edges associated to the Hex faces
    // edge 18 is a diagonal edge internal to the Hex
    int tetsEdgesInHex[6][6] = {
        {16,  2,  3, 18, 10 ,14},
        { 3, 11, 15, 18, 14,  6},
        {15,  7,  8, 18,  6, 17},
        {12,  5, 18,  8,  4, 17},
        { 9,  5, 13,  0, 12, 18},
        {13, 10,  1,  0, 18, 16}
    };

    // faces from 0 to 11 are associated to the Hex faces (each of the 4 Hex faces is split in two triangles),
    // faces from 12 to 17, are faces internal to the Hex
    int tetsFacesInHex[6][4] = {
        {12, 2, 13, 10},
        {13, 8, 14, 3},
        {14, 11, 15, 9},
        {6, 5, 15, 16},
        {0, 16, 17, 7},
        {17, 12, 4, 1}
    };


    subElemToBrickElementNodesMap_.resize(numSubElemsPerBrickElement_);
    for (int i=0; i<numSubElemsPerBrickElement_; ++i) {
      subElemToBrickElementNodesMap_[i].resize(4);
      for (int j=0; j<4; ++j)
        subElemToBrickElementNodesMap_[i][j]=tetsInHex[i][j];
    }


    subElemToBrickElementFacesMap_.resize(numSubElemsPerBrickElement_);
    for (int i=0; i<numSubElemsPerBrickElement_; ++i) {
      subElemToBrickElementFacesMap_[i].resize(4);
      for (int j=0; j<4; ++j)
        subElemToBrickElementFacesMap_[i][j] = tetsFacesInHex[i][j];
    }

    subElemToBrickElementEdgesMap_.resize(numSubElemsPerBrickElement_);
    for (int i=0; i<numSubElemsPerBrickElement_; ++i) {
      subElemToBrickElementEdgesMap_[i].resize(6);
      for (int j=0; j<6; ++j)
        subElemToBrickElementEdgesMap_[i][j] = tetsEdgesInHex[i][j];
    }
  } break;
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true,
            std::logic_error,
            "CartesianConnManager: element not supported" << std::endl);
  }

  elemTopology_ = elemTopo;

  buildLocalElements();
}

void CartesianConnManager::buildConnectivity(const panzer::FieldPattern & fp)
{
  int numIds = fp.numberIds();

  connectivity_.clear();
  connectivity_.resize(totalElements_);

  std::vector<std::string> elementBlockIds;
  getElementBlockIds(elementBlockIds);

  // loop over all element blocks
  for(std::size_t i=0;i<elementBlockIds.size();i++) {

    // loop over each element block
    const std::vector<LocalOrdinal> & elmts = getElementBlock(elementBlockIds[i]);
    for(std::size_t e=0;e<elmts.size();e++) {
      LocalOrdinal element = elmts[e];
      std::vector<GlobalOrdinal> & conn = connectivity_[element];

      int ordered_dim[4] = { 0, 1, 2, 3 }; // this is the order of creation
      for(int d=0;d<4;d++) {
        int od = ordered_dim[d];
        if(dim_==2) {
          if(od!=3) // ordered dimension in 2D is incorrect
            updateConnectivity_2d(fp,od,element,conn);
        }
        else
          updateConnectivity_3d(fp,od,element,conn);
      }

      TEUCHOS_ASSERT_EQUALITY(numIds,Teuchos::as<int>(conn.size()));
    }
  }
}


Teuchos::RCP<panzer::ConnManager>
CartesianConnManager::noConnectivityClone() const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<CartesianConnManager> clone = rcp(new CartesianConnManager);

  clone->numProc_ = numProc_;
  clone->myRank_ = myRank_;

  clone->dim_ = dim_;

  clone->brickElements_.x = brickElements_.x;
  clone->brickElements_.y = brickElements_.y;
  clone->brickElements_.z = brickElements_.z;

  clone->processors_.x = processors_.x;
  clone->processors_.y = processors_.y;
  clone->processors_.z = processors_.z;

  clone->blocks_.x = blocks_.x;
  clone->blocks_.y = blocks_.y;
  clone->blocks_.z = blocks_.z;

  clone->totalBrickElements_.x = totalBrickElements_.x;
  clone->totalBrickElements_.y = totalBrickElements_.y;
  clone->totalBrickElements_.z = totalBrickElements_.z;

  clone->totalNodes_ = totalNodes_;
  clone->totalEdges_ = totalEdges_;
  clone->totalFaces_ = totalFaces_;
  clone->totalElements_ = totalElements_;
  clone->numSubElemsPerBrickElement_ = numSubElemsPerBrickElement_;
  clone->subElemToBrickElementNodesMap_ = subElemToBrickElementNodesMap_;
  clone->subElemToBrickElementEdgesMap_ = subElemToBrickElementEdgesMap_;
  clone->subElemToBrickElementFacesMap_ = subElemToBrickElementFacesMap_;
  clone->elemTopology_ = elemTopology_;

  clone->buildLocalElements();

  return clone;
}

std::size_t CartesianConnManager::numElementBlocks() const
{
  return Teuchos::as<std::size_t>(blocks_.x * blocks_.y * blocks_.z);
}

void CartesianConnManager::
getElementBlockIds(std::vector<std::string> & elementBlockIds) const
{
  for(int i=0;i<blocks_.x;i++) {
    for(int j=0;j<blocks_.y;j++) {
      for(int k=0;k<blocks_.z;k++) {
        // build up block name
        std::stringstream ss;
        ss << "eblock-" << i << "_" << j;
        if(dim_==3)
          ss << "_" << k;

        elementBlockIds.push_back(ss.str());
      }
    }
  }
}

void
CartesianConnManager::getElementBlockTopologies(std::vector<shards::CellTopology> & elementBlockTopologies) const
{

  if ( dim_ == 2 ) {
    int nblocks = blocks_.x*blocks_.y;
    const CellTopologyData & myCellData = *elemTopology_.getCellTopologyData();
    struct shards::CellTopology my_topo(&myCellData);
    elementBlockTopologies = std::vector<shards::CellTopology>(nblocks, my_topo);
  }
  if ( dim_ == 3 ) {
    int nblocks = blocks_.x*blocks_.y*blocks_.z;
    const CellTopologyData & myCellData = *elemTopology_.getCellTopologyData();
    struct shards::CellTopology my_topo(&myCellData);
    elementBlockTopologies = std::vector<shards::CellTopology>(nblocks, my_topo);
  }

}

  const std::vector<ConnManager::LocalOrdinal> &
CartesianConnManager::getElementBlock(const std::string & blockId) const
{
  // find the element block
  auto found = localElements_.find(blockId);
  if(found!=localElements_.end())
    return found->second;

  // if not found return an empty vector
  return emptyVector_;
}

std::string
CartesianConnManager::getBlockId(ConnManager::LocalOrdinal localBrickElmentId) const
{
  auto globalTriplet = computeLocalBrickElementGlobalTriplet(localBrickElmentId,myBrickElements_,myBrickOffset_);

  int i = Teuchos::as<int>(globalTriplet.x/brickElements_.x);
  int j = Teuchos::as<int>(globalTriplet.y/brickElements_.y);
  int k = Teuchos::as<int>(globalTriplet.z/brickElements_.z);

  std::stringstream ss;
  ss << "eblock-" << i << "_" << j;
  if(dim_==3)
    ss << "_" << k;

  return ss.str();
}

typename CartesianConnManager::template Triplet<int>
CartesianConnManager::computeMyRankTriplet(int myRank,int dim,const Triplet<int> & procs)
{
  // for 2D the rank satisfies:  r = px * j + i
  // for 3D the rank satisfies:  r = px * py * k + px * j + i

  Triplet<int> index;

  if(dim==2) {
    TEUCHOS_ASSERT(procs.z==1); // sanity check

    index.z = 0; // by assumption
    index.y = myRank / procs.x;
    index.x = myRank - procs.x*index.y;
  }
  else {
    TEUCHOS_ASSERT(dim==3); // sanity check

    index.z =  myRank / (procs.x*procs.y);
    index.y = (myRank - procs.x*procs.y*index.z) / procs.x;
    index.x =  myRank - procs.x*procs.y*index.z - procs.x*index.y;
  }

  return index;
}

typename CartesianConnManager::template Triplet<ConnManager::GlobalOrdinal>
CartesianConnManager::
computeLocalBrickElementGlobalTriplet(int index,const Triplet<GlobalOrdinal> & myBrickElements,
                                           const Triplet<GlobalOrdinal> & myBrickOffset)
{
  Triplet<GlobalOrdinal> localIndex;

  localIndex.z =  index / (myBrickElements.x*myBrickElements.y);
  localIndex.y = (index - myBrickElements.x*myBrickElements.y*localIndex.z) / myBrickElements.x;
  localIndex.x =  index - myBrickElements.x*myBrickElements.y*localIndex.z - myBrickElements.x*localIndex.y;

  localIndex.x += myBrickOffset.x;
  localIndex.y += myBrickOffset.y;
  localIndex.z += myBrickOffset.z;

  return localIndex;
}

ConnManager::GlobalOrdinal
CartesianConnManager::computeGlobalBrickElementIndex(const Triplet<GlobalOrdinal> & bricklement,
                                                const Triplet<GlobalOrdinal> & brickShape)
{
  return brickShape.x * brickShape.y * bricklement.z + brickShape.x * bricklement.y + bricklement.x;
}

ConnManager::GlobalOrdinal
CartesianConnManager::computeGlobalBrickElementIndex(const Triplet<GlobalOrdinal> & brickElement) const
{
  return computeGlobalBrickElementIndex(brickElement,totalBrickElements_);
}

ConnManager::LocalOrdinal
CartesianConnManager::computeLocalBrickElementIndex(const Triplet<GlobalOrdinal> & brickElement,
                                               const Triplet<GlobalOrdinal> & myBrickElements,
                                               const Triplet<GlobalOrdinal> & myBrickOffset)
{
  // first make sure its in range
  GlobalOrdinal dx = brickElement.x-myBrickOffset.x;
  GlobalOrdinal dy = brickElement.y-myBrickOffset.y;
  GlobalOrdinal dz = brickElement.z-myBrickOffset.z;

  if(   dx>=myBrickElements.x || dx<0
     || dy>=myBrickElements.y || dy<0
     || dz>=myBrickElements.z || dz<0)
    return -1;

  return myBrickElements.x*myBrickElements.y*dz + myBrickElements.x*dy + dx;
}

ConnManager::LocalOrdinal
CartesianConnManager::computeLocalBrickElementIndex(const Triplet<GlobalOrdinal> & brickElement) const
{
  return computeLocalBrickElementIndex(brickElement,myBrickElements_,myBrickOffset_);
}

void
CartesianConnManager::buildLocalElements()
{
  // figure out largest possible even subdivision
  Triplet<GlobalOrdinal> evenBrickElements;
  evenBrickElements.x = totalBrickElements_.x / processors_.x;
  evenBrickElements.y = totalBrickElements_.y / processors_.y;
  evenBrickElements.z = totalBrickElements_.z / processors_.z;

  // compute the reaminder
  Triplet<GlobalOrdinal> remainderBrickElements;
  remainderBrickElements.x = totalBrickElements_.x - evenBrickElements.x * processors_.x;
  remainderBrickElements.y = totalBrickElements_.y - evenBrickElements.y * processors_.y;
  remainderBrickElements.z = totalBrickElements_.z - evenBrickElements.z * processors_.z;

  // compute the processor rank triplet
  Triplet<int> myRank = computeMyRankTriplet(myRank_,dim_,processors_);

  // how many elements in each direction owned by this proccessor
  Triplet<GlobalOrdinal> myBrickElements;
  myBrickElements.x = evenBrickElements.x + (remainderBrickElements.x > myRank.x ? 1 : 0);
  myBrickElements.y = evenBrickElements.y + (remainderBrickElements.y > myRank.y ? 1 : 0);
  myBrickElements.z = evenBrickElements.z + (remainderBrickElements.z > myRank.z ? 1 : 0);

  myBrickElements_ = myBrickElements;

  // figure out the offset for this processor
  Triplet<GlobalOrdinal> myBrickOffset;
  myBrickOffset.x = myRank.x*evenBrickElements.x + std::min(Teuchos::as<int>(remainderBrickElements.x),myRank.x);
  myBrickOffset.y = myRank.y*evenBrickElements.y + std::min(Teuchos::as<int>(remainderBrickElements.y),myRank.y);
  myBrickOffset.z = myRank.z*evenBrickElements.z + std::min(Teuchos::as<int>(remainderBrickElements.z),myRank.z);

  myBrickOffset_ = myBrickOffset;

  // loop over all your elements and assign them to an element block that owns them
  for(GlobalOrdinal i=0;i<numSubElemsPerBrickElement_*myBrickElements_.x*myBrickElements_.y*myBrickElements_.z;i++)
    localElements_[getBlockId(Teuchos::as<LocalOrdinal>(i/numSubElemsPerBrickElement_))].push_back(Teuchos::as<int>(i));
}

void
CartesianConnManager::updateConnectivity_2d(const panzer::FieldPattern & fp,
                                            int subcellDim,
                                            int localCellId,
                                            std::vector<GlobalOrdinal> & conn) const
{
  int localElementId = localCellId/numSubElemsPerBrickElement_;
  int subElementId = localCellId%numSubElemsPerBrickElement_;

  Triplet<GlobalOrdinal> index = computeLocalBrickElementGlobalTriplet(localElementId,myBrickElements_,myBrickOffset_);

  for(int c=0;c<fp.getSubcellCount(subcellDim);c++) {
    std::size_t num = fp.getSubcellIndices(subcellDim,c).size();
    TEUCHOS_ASSERT(num==1 || num==0);

    if(num==0) continue;

    // there is an index to add here
    if(subcellDim==0) {
      int n = subElemToBrickElementNodesMap_[subElementId][c];
      // node number     =  Number for the lower left
      GlobalOrdinal node =  (totalBrickElements_.x+1)*index.y + index.x
                           + (n==1 || n==2) * 1                    // shift for i+1
                           + (n==3 || n==2) * (totalBrickElements_.x+1); // shift for j+1
      conn.push_back(node);
    }
    else if(subcellDim==1) {
      int e = subElemToBrickElementEdgesMap_[subElementId][c];
      // node number     =  Number for the lower left
      GlobalOrdinal edge = totalNodes_ + (2*totalBrickElements_.x+1)*index.y + index.x
                                        + (        e==1) * 1
                                        + (        e==2) * (2*totalBrickElements_.x+1)
                                        + (e==1 || e==3) * totalBrickElements_.x;

      if((elemTopology_.getKey()==shards::Triangle<3>::key) && (e == 4)) {
        GlobalOrdinal totalElementEdges = (totalBrickElements_.x+1)*totalBrickElements_.y+totalBrickElements_.x*(totalBrickElements_.y+1);
        edge = totalNodes_ + totalElementEdges + computeGlobalBrickElementIndex(index);
      }
      conn.push_back(edge);
    }
    else if(subcellDim==2)
      conn.push_back(totalNodes_+totalEdges_+numSubElemsPerBrickElement_*computeGlobalBrickElementIndex(index)+subElementId);

    else {
      TEUCHOS_ASSERT(false);
    }
  }
}

void
CartesianConnManager::updateConnectivity_3d(const panzer::FieldPattern & fp,
                                            int subcellDim,
                                            int localCellId,
                                            std::vector<GlobalOrdinal> & conn) const
{
  int localElementId = localCellId/numSubElemsPerBrickElement_;
  int subElementId = localCellId%numSubElemsPerBrickElement_;

  Triplet<GlobalOrdinal> index = computeLocalBrickElementGlobalTriplet(localElementId,myBrickElements_,myBrickOffset_);

  //std::cout << "Totals: " << totalNodes_ << " " << totalEdges_ <<  " " << totalFaces_ << " | " ;
  for(int c=0;c<fp.getSubcellCount(subcellDim);c++) {
    std::size_t num = fp.getSubcellIndices(subcellDim,c).size();
    TEUCHOS_ASSERT(num==1 || num==0);

    if(num==0) continue;

    // there is an index to add here
    if(subcellDim==0) {
      int n = subElemToBrickElementNodesMap_[subElementId][c];
      // node number     =  Number for the lower left
      GlobalOrdinal node =  (totalBrickElements_.x+1)*(totalBrickElements_.y+1)*index.z + (totalBrickElements_.x+1)*index.y + index.x
                           + (n==1 || n==2 || n==5 || n==6) * 1                                          // shift for i+1
                           + (n==3 || n==2 || n==6 || n==7) * (totalBrickElements_.x+1)                       // shift for j+1
                           + (n==4 || n==5 || n==6 || n==7) * (totalBrickElements_.y+1)*(totalBrickElements_.x+1); // shift for k+1

      conn.push_back(node);
    }
    else if(subcellDim==1) {
      int e = subElemToBrickElementEdgesMap_[subElementId][c];
      GlobalOrdinal kshift = (totalBrickElements_.x+1)*totalBrickElements_.y +
                             totalBrickElements_.x*(totalBrickElements_.y+1) +
                             (totalBrickElements_.x+1)*(totalBrickElements_.y+1);

      GlobalOrdinal basePoint = index.x+index.y*(2*totalBrickElements_.x+1)+
                                index.z*kshift;

      GlobalOrdinal edge =  totalNodes_ + basePoint;

      // horizontal edges: bottom
      if(e== 0) edge += 0;
      if(e== 1) edge += totalBrickElements_.x+1;
      if(e== 2) edge += 2*totalBrickElements_.x+1;
      if(e== 3) edge += totalBrickElements_.x;

      // horizontal edges: top
      if(e== 4) edge += kshift + 0;
      if(e== 5) edge += kshift + totalBrickElements_.x+1;
      if(e== 6) edge += kshift + 2*totalBrickElements_.x+1;
      if(e== 7) edge += kshift + totalBrickElements_.x;

      // vertical edges
      if(e==8 || e==9 || e==10 || e==11)
        edge +=  (totalBrickElements_.x+1)*totalBrickElements_.y + totalBrickElements_.x*(totalBrickElements_.y+1)
                - index.y*totalBrickElements_.x;

      if(e== 8) edge += 0;
      if(e== 9) edge += 1;
      if(e==10) edge += totalBrickElements_.x+2;
      if(e==11) edge += totalBrickElements_.x+1;

      if((elemTopology_.getKey()==shards::Tetrahedron<4>::key) && (e >= 12)) {
        GlobalOrdinal totalElementEdges = (totalBrickElements_.x+1)*totalBrickElements_.y*(totalBrickElements_.z+1)+
            totalBrickElements_.x*(totalBrickElements_.y+1)*(totalBrickElements_.z+1)
            +(totalBrickElements_.x+1)*(totalBrickElements_.y+1)*totalBrickElements_.z;
        GlobalOrdinal faceKshift = totalBrickElements_.x*totalBrickElements_.y
                    +(totalBrickElements_.x+1)*totalBrickElements_.y
                    +totalBrickElements_.x*(totalBrickElements_.y+1);
        basePoint = totalElementEdges + index.x+index.y*totalBrickElements_.x+index.z*faceKshift;
        edge =  totalNodes_ + basePoint;

        // vertical faces
        if(e==12 || e==13 || e==14 || e==15)
          edge += totalBrickElements_.x*totalBrickElements_.y + (index.y+1)*(totalBrickElements_.x+1);

        if(e==12) edge += -totalBrickElements_.x-1;
        if(e==13) edge += 0;
        if(e==14) edge += totalBrickElements_.x;
        if(e==15) edge += -1;
        if(e==16) edge += 0;
        if(e==17) edge += faceKshift; // move it up a level

        if(e == 18)
          edge =  totalNodes_ + 2*totalElementEdges + computeGlobalBrickElementIndex(index);
      }
      conn.push_back(edge);
    }
    else if(subcellDim==2) {
      int f = subElemToBrickElementFacesMap_[subElementId][c];
      GlobalOrdinal kshift =  totalBrickElements_.x*totalBrickElements_.y
                             +(totalBrickElements_.x+1)*totalBrickElements_.y
                             +totalBrickElements_.x*(totalBrickElements_.y+1);
      GlobalOrdinal basePoint = index.x+index.y*totalBrickElements_.x+index.z*kshift;
      GlobalOrdinal face = totalNodes_+totalEdges_+basePoint;

      // vertical faces
      if(f==0 || f==1 || f==2 || f==3)
        face += totalBrickElements_.x*totalBrickElements_.y + (index.y+1)*(totalBrickElements_.x+1);

      if(f==0) face += -totalBrickElements_.x-1;
      if(f==1) face += 0;
      if(f==2) face += totalBrickElements_.x;
      if(f==3) face += -1;
      if(f==4) face += 0;
      if(f==5) face += kshift; // move it up a level

      if((elemTopology_.getKey()==shards::Tetrahedron<4>::key) && (f >= 6)) {

        GlobalOrdinal totalElementFaces = totalBrickElements_.x*totalBrickElements_.y*(totalBrickElements_.z+1)
            + totalBrickElements_.x*(totalBrickElements_.y+1)*totalBrickElements_.z
            + (totalBrickElements_.x+1)*totalBrickElements_.y*totalBrickElements_.z;

        face = totalNodes_+totalEdges_+totalElementFaces+basePoint;

        // vertical faces
        if(f==6 || f==7 || f==8 || f==9)
          face += totalBrickElements_.x*totalBrickElements_.y + (index.y+1)*(totalBrickElements_.x+1);

        if(f==6) face += -totalBrickElements_.x-1;
        if(f==7) face += 0;
        if(f==8) face += totalBrickElements_.x;
        if(f==9) face += -1;
        if(f==10) face += 0;
        if(f==11) face += kshift; // move it up a level

        if( f >= 12) {
          face = totalNodes_+totalEdges_+2*totalElementFaces;
          face += 6*computeGlobalBrickElementIndex(index)+(f-12);

        }
      }
      conn.push_back(face);
    }
    else if(subcellDim==3) {
      conn.push_back(totalNodes_+totalEdges_+totalFaces_+numSubElemsPerBrickElement_*computeGlobalBrickElementIndex(index)+subElementId);
    }
    else {
      TEUCHOS_ASSERT(false);
    }
  }
}

} // end unit test
} // end panzer
