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

#ifndef __CartesianConnManager_impl_hpp__
#define __CartesianConnManager_impl_hpp__

#include "CartesianConnManager.hpp"

#include "Panzer_FieldPattern.hpp"
#include "Shards_CellTopology.hpp"
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

namespace panzer {
namespace unit_test {

template <typename LocalOrdinal,typename GlobalOrdinal>
void CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
initialize(const Teuchos::MpiComm<int> & comm,GlobalOrdinal nx, GlobalOrdinal ny,
                                              int px, int py,
                                              int bx, int by)
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

  elements_.x = nx;
  elements_.y = ny;
  elements_.z = 1;

  processors_.x = px;
  processors_.y = py;
  processors_.z = 1;  // even though its 2D do this

  blocks_.x = bx;
  blocks_.y = by;
  blocks_.z = 1;  // even though its 2D do this

  totalElements_.x = nx*bx;
  totalElements_.y = ny*by;
  totalElements_.z = 1*1;

  totalNodes_ = (bx*nx+1)*(by*ny+1);
  totalEdges_ = (bx*nx+1)*(by*ny)+(bx*nx)*(by*ny+1);
  totalFaces_ = (bx*nx)*(by*ny);
  totalCells_ = totalFaces_;

  buildLocalElements();
}

template <typename LocalOrdinal,typename GlobalOrdinal>
void CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
initialize(const Teuchos::MpiComm<int> & comm,GlobalOrdinal nx, GlobalOrdinal ny, GlobalOrdinal nz,
                                              int px, int py,int pz,
                                              int bx, int by,int bz)
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

  elements_.x = nx;
  elements_.y = ny;
  elements_.z = nz;

  processors_.x = px;
  processors_.y = py;
  processors_.z = pz;

  blocks_.x = bx;
  blocks_.y = by;
  blocks_.z = bz;

  totalElements_.x = nx*bx;
  totalElements_.y = ny*by;
  totalElements_.z = nz*bz;

  totalNodes_ = (bx*nx+1)*(by*ny+1)*(bz*nz+1);
  totalEdges_ = (bx*nx+1)*(by*ny)*(bz*nz+1)+(bx*nx)*(by*ny+1)*(bz*nz+1)+(bx*nx+1)*(by*ny+1)*(bz*nz);
  totalFaces_ = (bx*nx)*(by*ny)*(bz*nz+1) + (bx*nx)*(by*ny+1)*(bz*nz) + (bx*nx+1)*(by*ny)*(bz*nz);
  totalCells_ = (bx*nx)*(by*ny)*(bz*nz);

  buildLocalElements();
}

template <typename LocalOrdinal,typename GlobalOrdinal>
void 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
buildConnectivity(const panzer::FieldPattern & fp)
{
  int numIds = fp.numberIds();

  connectivity_.clear();
  connectivity_.resize(myElements_.x*myElements_.y*myElements_.z);
  
  std::vector<std::string> elementBlockIds;
  getElementBlockIds(elementBlockIds);

  // loop over all element blocks
  for(std::size_t i=0;i<elementBlockIds.size();i++) {

    // loop over each element block 
    const std::vector<LocalOrdinal> & elmts = getElementBlock(elementBlockIds[i]);
    for(std::size_t e=0;e<elmts.size();e++) {
      LocalOrdinal element = elmts[e];
      std::vector<GlobalOrdinal> & conn = connectivity_[element];
      
      // loop over each dimension, add in coefficients
      //int ordered_dim[4] = { 0, 1, 3, 2 }; // this is the order of creation
                                           // because of the way shards orders
                                           // 3D elements
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


template <typename LocalOrdinal,typename GlobalOrdinal>
Teuchos::RCP<panzer::ConnManagerBase<LocalOrdinal> > 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
noConnectivityClone() const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<CartesianConnManager<LocalOrdinal,GlobalOrdinal> > clone
    = rcp(new CartesianConnManager<LocalOrdinal,GlobalOrdinal>);

  clone->numProc_ = numProc_;
  clone->myRank_ = myRank_;

  clone->dim_ = dim_;

  clone->elements_.x = elements_.x;
  clone->elements_.y = elements_.y;
  clone->elements_.z = elements_.z;

  clone->processors_.x = processors_.x;
  clone->processors_.y = processors_.y;
  clone->processors_.z = processors_.z;

  clone->blocks_.x = blocks_.x;
  clone->blocks_.y = blocks_.y;
  clone->blocks_.z = blocks_.z;

  clone->totalElements_.x = totalElements_.x;
  clone->totalElements_.y = totalElements_.y;
  clone->totalElements_.z = totalElements_.z;

  clone->totalNodes_ = totalNodes_;
  clone->totalEdges_ = totalEdges_;
  clone->totalFaces_ = totalFaces_;
  clone->totalCells_ = totalCells_;

  clone->buildLocalElements();
 
  return clone;
}

template <typename LocalOrdinal,typename GlobalOrdinal>
std::size_t 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
numElementBlocks() const
{
  return Teuchos::as<std::size_t>(blocks_.x * blocks_.y * blocks_.z);
}

template <typename LocalOrdinal,typename GlobalOrdinal>
void 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
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
template <typename LocalOrdinal,typename GlobalOrdinal>
void
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
getElementBlockTopologies(std::vector<shards::CellTopology> & elementBlockTopologies) const
{

  if ( dim_ == 2 ) {
    int nblocks = blocks_.x*blocks_.y;
    const CellTopologyData & myCellData = *shards::getCellTopologyData<shards::Quadrilateral<4> >();
    struct shards::CellTopology my_topo(&myCellData);
    elementBlockTopologies = std::vector<shards::CellTopology>(nblocks, my_topo);
  }
  if ( dim_ == 3 ) {
    int nblocks = blocks_.x*blocks_.y*blocks_.z;
    const CellTopologyData & myCellData = *shards::getCellTopologyData<shards::Hexahedron<8> >();
    struct shards::CellTopology my_topo(&myCellData);
    elementBlockTopologies = std::vector<shards::CellTopology>(nblocks, my_topo);
  }

}
template <typename LocalOrdinal,typename GlobalOrdinal>
const std::vector<LocalOrdinal> & 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
getElementBlock(const std::string & blockId) const
{
  // find the element block
  auto found = localElements_.find(blockId);
  if(found!=localElements_.end())
    return found->second;

  // if not found return an empty vector
  return emptyVector_;
}

template <typename LocalOrdinal,typename GlobalOrdinal>
std::string 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
getBlockId(LocalOrdinal localElmtId) const
{
  auto globalTriplet = computeLocalElementGlobalTriplet(localElmtId,myElements_,myOffset_);

  int i = Teuchos::as<int>(globalTriplet.x/elements_.x);
  int j = Teuchos::as<int>(globalTriplet.y/elements_.y);
  int k = Teuchos::as<int>(globalTriplet.z/elements_.z);

  std::stringstream ss;
  ss << "eblock-" << i << "_" << j;
  if(dim_==3) 
    ss << "_" << k;

  return ss.str();
}

template <typename LocalOrdinal,typename GlobalOrdinal>
typename CartesianConnManager<LocalOrdinal,GlobalOrdinal>::template Triplet<int> 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
computeMyRankTriplet(int myRank,int dim,const Triplet<int> & procs) 
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

template <typename LocalOrdinal,typename GlobalOrdinal>
typename CartesianConnManager<LocalOrdinal,GlobalOrdinal>::template Triplet<GlobalOrdinal> 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
computeLocalElementGlobalTriplet(int index,const Triplet<GlobalOrdinal> & myElements,
                                           const Triplet<GlobalOrdinal> & myOffset)
{
  Triplet<GlobalOrdinal> localIndex;

  localIndex.z =  index / (myElements.x*myElements.y);
  localIndex.y = (index - myElements.x*myElements.y*localIndex.z) / myElements.x;
  localIndex.x =  index - myElements.x*myElements.y*localIndex.z - myElements.x*localIndex.y;

  localIndex.x += myOffset.x;
  localIndex.y += myOffset.y;
  localIndex.z += myOffset.z;

  return localIndex;
}

template <typename LocalOrdinal,typename GlobalOrdinal>
GlobalOrdinal
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
computeGlobalElementIndex(const Triplet<GlobalOrdinal> & element,const Triplet<GlobalOrdinal> & shape)
{
  return shape.x * shape.y * element.z + shape.x * element.y + element.x;
}

template <typename LocalOrdinal,typename GlobalOrdinal>
GlobalOrdinal 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
computeGlobalElementIndex(const Triplet<GlobalOrdinal> & element) const
{
  return computeGlobalElementIndex(element,totalElements_);
}

template <typename LocalOrdinal,typename GlobalOrdinal>
LocalOrdinal 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
computeLocalElementIndex(const Triplet<GlobalOrdinal> & element,
                         const Triplet<GlobalOrdinal> & myElements,
                         const Triplet<GlobalOrdinal> & myOffset)
{
  // first make sure its in range
  GlobalOrdinal dx = element.x-myOffset.x;  
  GlobalOrdinal dy = element.y-myOffset.y;  
  GlobalOrdinal dz = element.z-myOffset.z;  

  if(   dx>=myElements.x || dx<0
     || dy>=myElements.y || dy<0
     || dz>=myElements.z || dz<0)
    return -1;

  return myElements.x*myElements.y*dz + myElements.x*dy + dx;
}

template <typename LocalOrdinal,typename GlobalOrdinal>
LocalOrdinal 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
computeLocalElementIndex(const Triplet<GlobalOrdinal> & element) const
{
  return computeLocalElementIndex(element,myElements_,myOffset_);
}

template <typename LocalOrdinal,typename GlobalOrdinal>
void 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
buildLocalElements()
{
  // figure out largest possible even subdivision
  Triplet<GlobalOrdinal> evenElements;
  evenElements.x = totalElements_.x / processors_.x;
  evenElements.y = totalElements_.y / processors_.y;
  evenElements.z = totalElements_.z / processors_.z;

  // compute the reaminder
  Triplet<GlobalOrdinal> remainderElements;
  remainderElements.x = totalElements_.x - evenElements.x * processors_.x;
  remainderElements.y = totalElements_.y - evenElements.y * processors_.y;
  remainderElements.z = totalElements_.z - evenElements.z * processors_.z;

  // compute the processor rank triplet
  Triplet<int> myRank = computeMyRankTriplet(myRank_,dim_,processors_);

  // how many elements in each direction owned by this proccessor
  Triplet<GlobalOrdinal> myElements;
  myElements.x = evenElements.x + (remainderElements.x > myRank.x ? 1 : 0);
  myElements.y = evenElements.y + (remainderElements.y > myRank.y ? 1 : 0);
  myElements.z = evenElements.z + (remainderElements.z > myRank.z ? 1 : 0);

  myElements_ = myElements;

  // figure out the offset for this processor
  Triplet<GlobalOrdinal> myOffset;
  myOffset.x = myRank.x*evenElements.x + std::min(Teuchos::as<int>(remainderElements.x),myRank.x);
  myOffset.y = myRank.y*evenElements.y + std::min(Teuchos::as<int>(remainderElements.y),myRank.y);
  myOffset.z = myRank.z*evenElements.z + std::min(Teuchos::as<int>(remainderElements.z),myRank.z);

  myOffset_ = myOffset;

  // loop over all your elements and assign them to an element block that owns them
  for(GlobalOrdinal i=0;i<myElements_.x*myElements_.y*myElements_.z;i++)
    localElements_[getBlockId(Teuchos::as<GlobalOrdinal>(i))].push_back(Teuchos::as<int>(i)); 
}

template <typename LocalOrdinal,typename GlobalOrdinal>
void 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
updateConnectivity_2d(const panzer::FieldPattern & fp,
                      int subcellDim,
                      int localElementId,
                      std::vector<GlobalOrdinal> & conn) const
{
  Triplet<GlobalOrdinal> index = computeLocalElementGlobalTriplet(localElementId,myElements_,myOffset_);
 
  for(int c=0;c<fp.getSubcellCount(subcellDim);c++) {
    std::size_t num = fp.getSubcellIndices(subcellDim,c).size();
    TEUCHOS_ASSERT(num==1 || num==0);
 
    if(num==0) continue;
 
    // there is an index to add here
    if(subcellDim==0) {
      // node number     =  Number for the lower left          
      GlobalOrdinal node =  (totalElements_.x+1)*index.y + index.x
                           + (c==1 || c==2) * 1                    // shift for i+1
                           + (c==3 || c==2) * (totalElements_.x+1); // shift for j+1
      conn.push_back(node);
    }
    else if(subcellDim==1) {
      // node number     =  Number for the lower left          
      GlobalOrdinal edge =  totalNodes_ + (2*totalElements_.x+1)*index.y + index.x
                                        + (        c==1) * 1                     
                                        + (        c==2) * (2*totalElements_.x+1)
                                        + (c==1 || c==3) * totalElements_.x;
      conn.push_back(edge);
    } 
    else if(subcellDim==2)
      conn.push_back(totalNodes_+totalEdges_+computeGlobalElementIndex(index));
    else {
      TEUCHOS_ASSERT(false);
    }
  }
}

template <typename LocalOrdinal,typename GlobalOrdinal>
void 
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::
updateConnectivity_3d(const panzer::FieldPattern & fp,
                      int subcellDim,
                      int localElementId,
                      std::vector<GlobalOrdinal> & conn) const
{
  Triplet<GlobalOrdinal> index = computeLocalElementGlobalTriplet(localElementId,myElements_,myOffset_);
 
  for(int c=0;c<fp.getSubcellCount(subcellDim);c++) {
    std::size_t num = fp.getSubcellIndices(subcellDim,c).size();
    TEUCHOS_ASSERT(num==1 || num==0);
 
    if(num==0) continue;
 
    // there is an index to add here
    if(subcellDim==0) {
      // node number     =  Number for the lower left          
      GlobalOrdinal node =  (totalElements_.x+1)*(totalElements_.y+1)*index.z + (totalElements_.x+1)*index.y + index.x
                           + (c==1 || c==2 || c==5 || c==6) * 1                                          // shift for i+1
                           + (c==3 || c==2 || c==6 || c==7) * (totalElements_.x+1)                       // shift for j+1
                           + (c==4 || c==5 || c==6 || c==7) * (totalElements_.y+1)*(totalElements_.x+1); // shift for k+1

      conn.push_back(node);
    }
    else if(subcellDim==1) {
      GlobalOrdinal basePoint = index.x+index.y*(2*totalElements_.x+1)+
                                index.z*(totalElements_.x+1)*totalElements_.y +
                                index.z*totalElements_.x*(totalElements_.y+1) +
                                index.z*(totalElements_.x+1)*(totalElements_.y+1);
      GlobalOrdinal kshift = (totalElements_.x+1)*totalElements_.y +
                             totalElements_.x*(totalElements_.y+1) +
                             (totalElements_.x+1)*(totalElements_.y+1);
      GlobalOrdinal edge =  totalNodes_ + basePoint;

      // horizontal edges: bottom
      if(c== 0) edge += 0;
      if(c== 1) edge += totalElements_.x+1;
      if(c== 2) edge += 2*totalElements_.x+1;
      if(c== 3) edge += totalElements_.x;

      // horizontal edges: top
      if(c== 4) edge += kshift + 0;
      if(c== 5) edge += kshift + totalElements_.x+1;
      if(c== 6) edge += kshift + 2*totalElements_.x+1;
      if(c== 7) edge += kshift + totalElements_.x;

      // vertical edges
      if(c==8 || c==9 || c==10 || c==11)
        edge +=  (totalElements_.x+1)*totalElements_.y + totalElements_.x*(totalElements_.y+1) 
                - index.y*totalElements_.x;

      if(c== 8) edge += 0;
      if(c== 9) edge += 1;
      if(c==10) edge += totalElements_.x+2;
      if(c==11) edge += totalElements_.x+1;

      conn.push_back(edge);
    } 
    else if(subcellDim==2) {
      GlobalOrdinal kshift =  totalElements_.x*totalElements_.y
                             +(totalElements_.x+1)*totalElements_.y
                             +totalElements_.x*(totalElements_.y+1);
      GlobalOrdinal basePoint = index.x+index.y*totalElements_.x+index.z*kshift;
      GlobalOrdinal face = totalNodes_+totalEdges_+basePoint;

      // vertical faces
      if(c==0 || c==1 || c==2 || c==3)
        face += totalElements_.x*totalElements_.y + (index.y+1)*(totalElements_.x+1);

      if(c==0) face += -totalElements_.x-1;
      if(c==1) face += 0;
      if(c==2) face += totalElements_.x;
      if(c==3) face += -1;
      if(c==4) face += 0;
      if(c==5) face += kshift; // move it up a level

      conn.push_back(face);
    }
    else if(subcellDim==3) {
      conn.push_back(totalNodes_+totalEdges_+totalFaces_+computeGlobalElementIndex(index));
    }
    else {
      TEUCHOS_ASSERT(false);
    }
  }
}

} // end unit test
} // end panzer

#endif
