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
  TEUCHOS_ASSERT(comm.getSize()==px*py);
  TEUCHOS_ASSERT(0<bx*by);

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
  TEUCHOS_ASSERT(comm.getSize()==px*py*pz);
  TEUCHOS_ASSERT(0<bx*by*bz);

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

  totalNodes_ = (bx*nx+1)*(by*ny+1)*(bz*nx+1);
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
  int dim = fp.getDimension();
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
      // ASSUME: Field pattern has logical ordering nodes, edges, faces, cells
      for(int d=0;d<=dim;d++)
        updateConnectivity(fp,d,element,conn);

      TEUCHOS_ASSERT(numIds==Teuchos::as<int>(conn.size()));
    }
  }
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

  // int i = Teuchos::as<int>(globalTriplet.x/blocks_.x);
  // int j = Teuchos::as<int>(globalTriplet.y/blocks_.y);
  // int k = Teuchos::as<int>(globalTriplet.z/blocks_.z);
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
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::Triplet<int> 
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
CartesianConnManager<LocalOrdinal,GlobalOrdinal>::Triplet<GlobalOrdinal> 
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
updateConnectivity(const panzer::FieldPattern & fp,
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
  }
}

} // end unit test
} // end panzer

#endif
