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
#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_MESHDATABASE_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_MESHDATABASE_HPP

#include <iostream>
#include <iterator>
#include <set>
#include <sstream>

#include "Teuchos_Comm.hpp"

#include "fem_assembly_typedefs.hpp"


namespace TpetraExamples
{


struct LLA
{
  LLA() {data[0]=0; data[1]=0;}

  LLA(global_ordinal_t i, global_ordinal_t j) {data[0]=i; data[1]=j;}

  global_ordinal_t operator[](global_ordinal_t i) const {return data[i];}
  global_ordinal_t & operator[](global_ordinal_t i) {return data[i];}

  global_ordinal_t data[2];
};



class MeshDatabase
{
public:
  MeshDatabase(Teuchos::RCP<const Teuchos::Comm<int> > comm,
               global_ordinal_t global_elements_x,
               global_ordinal_t global_elements_y,
               global_ordinal_t procs_x,
               global_ordinal_t procs_y);

  ~MeshDatabase(){}

  // Size accessors
  size_t getNumOwnedElements() const {return ownedElementGlobalIDs_.extent(0);}

  size_t getNumGhostElements() const {return ghostElementGlobalIDs_.extent(0);}

  size_t getNumOwnedNodes() const {return ownedNodeGlobalIDs_.extent(0);}

  size_t getNumGhostNodes() const {return ghostNodeGlobalIDs_.extent(0);}

  size_t getNumOwnedAndGhostNodes() const {return ownedAndGhostNodeGlobalIDs_.extent(0);}

  size_t getNumOwnedAndGhostElements() const {return ownedAndGhostElementGlobalIDs_.extent(0);}

  // Data accessors
  global_ordinal_view_t getOwnedElementGlobalIDs() {return ownedElementGlobalIDs_;}
  global_ordinal_view_t getGhostElementGlobalIDs() {return ghostElementGlobalIDs_;}
  global_ordinal_view_t getOwnedAndGhostElementGlobalIDs() {return ownedAndGhostElementGlobalIDs_;}

  global_ordinal_view_t getOwnedNodeGlobalIDs() {return ownedNodeGlobalIDs_;}
  global_ordinal_view_t getGhostNodeGlobalIDs() {return ghostNodeGlobalIDs_;}
  global_ordinal_view_t getOwnedAndGhostNodeGlobalIDs() {return ownedAndGhostNodeGlobalIDs_;}

  global_ordinal_2d_array_t getOwnedElementToNode() {return ownedElementToNode_;}
  global_ordinal_2d_array_t getGhostElementToNode() {return ghostElementToNode_;}

  // Debugging output
  void print(std::ostream & oss);

  inline bool nodeIsOwned(global_ordinal_t idx) {
    global_ordinal_t i,j;
    ij_from_idx(globalNodes_[0],idx,i,j);
    return nodeIsOwned(i,j);
  }

  inline bool elementIsOwned(global_ordinal_t idx) {
    global_ordinal_t i,j;
    ij_from_idx(globalElements_[0], idx, i, j);
    return elementIsOwned(i,j);
  }


private:

  inline bool nodeIsOwned(global_ordinal_t i, global_ordinal_t j) {
    return myNodeStart_[0] <= i &&  i < myNodeStop_[0] && myNodeStart_[1] <= j &&  j < myNodeStop_[1];
  }

  inline bool elementIsOwned(global_ordinal_t i, global_ordinal_t j) {
    return myElementStart_[0] <= i && i < myElementStop_[0] && myElementStart_[1] <= j && myElementStop_[1];
  }

  inline global_ordinal_t idx_from_ij(global_ordinal_t num_x, global_ordinal_t i, global_ordinal_t j) const {
    return j*num_x+i;
  }

  inline void ij_from_idx(global_ordinal_t num_x, global_ordinal_t idx, global_ordinal_t &i, global_ordinal_t&j) const {
    i = idx%num_x;
    j = (global_ordinal_t)((idx-i)/num_x);
  }

  void initializeOwnedAndGhostNodeGlobalIDs(void);

  void initializeOwnedAndGhostElementGlobalIDs(void);

  global_ordinal_view_t ownedElementGlobalIDs_;
  global_ordinal_view_t ghostElementGlobalIDs_;

  global_ordinal_view_t ownedNodeGlobalIDs_;
  global_ordinal_view_t ghostNodeGlobalIDs_;

  global_ordinal_view_t ownedAndGhostNodeGlobalIDs_;
  global_ordinal_view_t ownedAndGhostElementGlobalIDs_;

  global_ordinal_2d_array_t ownedElementToNode_;
  global_ordinal_2d_array_t ghostElementToNode_;

  // Global Mesh Info
  LLA globalElements_;
  LLA globalNodes_;
  LLA globalProcs_;

  // Local Mesh Info
  global_ordinal_t myElementStart_[2];
  global_ordinal_t myElementStop_[2];
  global_ordinal_t myNodeStart_[2];
  global_ordinal_t myNodeStop_[2];

  // Comm info
  Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  size_t MyRank_;
  LLA myProcIJ_;
};



MeshDatabase::MeshDatabase(Teuchos::RCP<const Teuchos::Comm<int> > comm,
                           global_ordinal_t global_elements_x,
                           global_ordinal_t global_elements_y,
                           global_ordinal_t procs_x,
                           global_ordinal_t procs_y):
      globalElements_(global_elements_x,global_elements_y),
      globalNodes_(global_elements_x+1,global_elements_y+1),
      globalProcs_(procs_x,procs_y),
      comm_(comm)
{

  // NOTE: Elements/nodes are numbered sequentially with x as the "fast" direction

  // Get processor decomp information
  MyRank_ = comm_->getRank();
  ij_from_idx(globalProcs_[0],MyRank_,myProcIJ_[0],myProcIJ_[1]);

  // Get local element & node start / stop
  global_ordinal_t num_my_elements=1, num_my_nodes=1;
  for(int k=0; k<2; k++) {
    global_ordinal_t eper = globalElements_[k] / globalProcs_[k];

    myElementStart_[k] = myProcIJ_[k] * eper;
    myElementStop_[k]  = (myProcIJ_[k] == globalProcs_[k]-1) ? globalElements_[k] : (myProcIJ_[k]+1) * eper;
    num_my_elements *= (myElementStop_[k]-myElementStart_[k]);

    myNodeStart_[k] = myProcIJ_[k] * eper;
    myNodeStop_[k]  = (myProcIJ_[k] == globalProcs_[k]-1) ? globalNodes_[k] : (myProcIJ_[k]+1) * eper;
    num_my_nodes *= (myNodeStop_[k]-myNodeStart_[k]);
  }

  // Generate the owned element ids
  Kokkos::resize(ownedElementGlobalIDs_,num_my_elements);
  int ect=0;
  for(global_ordinal_t j=myElementStart_[1]; j<myElementStop_[1]; j++) {
    for(global_ordinal_t i=myElementStart_[0]; i<myElementStop_[0]; i++) {
      global_ordinal_t idx=idx_from_ij(globalElements_[0],i,j);
      ownedElementGlobalIDs_(ect) = idx;
      ect++;
    }
  }

  // Generate the owned node ids
  Kokkos::resize(ownedNodeGlobalIDs_,num_my_nodes);
  int nct=0;
  for(global_ordinal_t j=myNodeStart_[1]; j<myNodeStop_[1]; j++) {
    for(global_ordinal_t i=myNodeStart_[0]; i<myNodeStop_[0]; i++) {
      global_ordinal_t idx=idx_from_ij(globalNodes_[0],i,j);
      ownedNodeGlobalIDs_(nct) = idx;
      nct++;
    }
  }

  // Generate the element-to-node map
  // NOTE: Hardwired to QUAD4's.  Nodes are ordered exodus-style (counter-clockwise) within an element
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  Kokkos::resize(ownedElementToNode_,num_my_elements,4);
#else
  Kokkos::resize(ownedElementToNode_,num_my_elements);
#endif
  int cct=0;
  for(global_ordinal_t j=myElementStart_[1]; j<myElementStop_[1]; j++) {
    for(global_ordinal_t i=myElementStart_[0]; i<myElementStop_[0]; i++) {
      // The (i,j) of the bottom left corner matches for elements & nodes
      global_ordinal_t nidx=idx_from_ij(globalNodes_[0],i,j);

      ownedElementToNode_(cct,0) = nidx;
      ownedElementToNode_(cct,1) = nidx+1;
      ownedElementToNode_(cct,2) = nidx+globalNodes_[0]+1;
      ownedElementToNode_(cct,3) = nidx+globalNodes_[0];
      cct++;
    }
  }

  // Generate the list of "ghost" elements & ghostElement2NodeMap
  // NOTE: This only generates a halo for elements where I own at least one node.  On the x/y hi sides,
  // the highers element does not own all the nodes on that proc.  Ergo, no halo in that direction
  std::vector<global_ordinal_t> my_ghost_elements;
  for(global_ordinal_t j=myElementStart_[1]-1; j<myElementStop_[1]; j++) {
    if(j<0 || j>=globalElements_[1]) continue; // Ignore stuff off the mesh
    for(global_ordinal_t i=myElementStart_[0]-1; i<myElementStop_[0]; i++) {
      if(i<0 || i>=globalElements_[0]) continue; // Ignore stuff off the mesh

      // Ignore proc interior
      if( j>myElementStart_[1]-1 && j<myElementStop_[1] && i>myElementStart_[0]-1 && i<myElementStop_[0])
        continue;

      global_ordinal_t idx=idx_from_ij(globalElements_[0],i,j);
      my_ghost_elements.push_back(idx);
    }
  }

  // NOTE: This are not recorded in Aztec/Ifpack/ML ordering.  Because most apps don't do that.
  Kokkos::resize(ghostElementGlobalIDs_,my_ghost_elements.size());
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  Kokkos::resize(ghostElementToNode_,my_ghost_elements.size(),4);
#else
  Kokkos::resize(ghostElementToNode_,my_ghost_elements.size());
#endif
  for(size_t k=0; k<my_ghost_elements.size(); k++) {
    global_ordinal_t i,j, eidx= my_ghost_elements[k];
    ghostElementGlobalIDs_(k) = eidx;
    ij_from_idx(globalElements_[0],eidx,i,j);

    // The (i,j) of the bottom left corner matches for elements & nodes
    global_ordinal_t nidx=idx_from_ij(globalNodes_[0],i,j);

    ghostElementToNode_(k,0) = nidx;
    ghostElementToNode_(k,1) = nidx+1;
    ghostElementToNode_(k,2) = nidx+globalNodes_[0]+1;
    ghostElementToNode_(k,3) = nidx+globalNodes_[0];
  }

  // Generate the list of "ghost" nodes (aka any node that exists on the ownedElement list that isn't owned
  std::set<global_ordinal_t> my_ghost_nodes;
  for(size_t k=0; k<ownedElementToNode_.extent(0); k++) {
    for(size_t l=0; l<ownedElementToNode_.extent(1); l++) {
      global_ordinal_t nidx=ownedElementToNode_(k,l);
      if(!nodeIsOwned(nidx)) {
        my_ghost_nodes.insert(nidx);
      }
    }
  }

  Kokkos::resize(ghostNodeGlobalIDs_,my_ghost_nodes.size());
  for(auto k=my_ghost_nodes.begin(); k!=my_ghost_nodes.end(); k++) {
    size_t kk = std::distance(my_ghost_nodes.begin(),k);
    ghostNodeGlobalIDs_(kk) = *k;
  }

  initializeOwnedAndGhostNodeGlobalIDs();
  initializeOwnedAndGhostElementGlobalIDs();
}


void MeshDatabase::initializeOwnedAndGhostNodeGlobalIDs(void)
{
  size_t total_size = getNumOwnedNodes() + getNumGhostNodes();

  Kokkos::resize(ownedAndGhostNodeGlobalIDs_, total_size);

  size_t insert_idx = 0;
  for(size_t idx=0; idx < getNumOwnedNodes(); idx++)
  {
    ownedAndGhostNodeGlobalIDs_(insert_idx++) = getOwnedNodeGlobalIDs()(idx);
  }
  for(size_t idx=0; idx < getNumGhostNodes(); idx++)
  {
    ownedAndGhostNodeGlobalIDs_(insert_idx++) = getGhostNodeGlobalIDs()(idx);
  }
}


void MeshDatabase::initializeOwnedAndGhostElementGlobalIDs(void)
{
  size_t total_size = getNumOwnedElements() + getNumGhostElements();
  Kokkos::resize(ownedAndGhostElementGlobalIDs_, total_size);

  size_t insert_idx = 0;
  for(size_t idx=0; idx<getNumOwnedElements(); idx++)
  {
    ownedAndGhostElementGlobalIDs_(insert_idx++) = getOwnedElementGlobalIDs()(idx);
  }
  for(size_t idx=0; idx<getNumGhostElements(); idx++)
  {
    ownedAndGhostElementGlobalIDs_(insert_idx++) = getGhostElementGlobalIDs()(idx);
  }
}



void MeshDatabase::print(std::ostream & oss)
{
  std::ostringstream ss;
  ss<<"["<<MyRank_<<","<<myProcIJ_[0]<<","<<myProcIJ_[1]<<"]";
  oss<<ss.str()<<" Global Elements = ["<<globalElements_[0]<<"x"<<globalElements_[1]<<"] Nodes ="<<globalNodes_[0]<<"x"<<globalNodes_[1]<<"]\n";
  oss<<ss.str()<<" Stop/Start Elements   = ["<<myElementStart_[0]<<","<<myElementStop_[0]<<")x["<<myElementStart_[1]<<","<<myElementStop_[1]<<")\n";
  oss<<ss.str()<<" Stop/Start Nodes      = ["<<myNodeStart_[0]<<","<<myNodeStop_[0]<<")x["<<myNodeStart_[1]<<","<<myNodeStop_[1]<<")\n";

  oss<<ss.str()<<" Owned Global Elements = ";
  for(size_t i=0; i<ownedElementGlobalIDs_.extent(0); i++) {
    oss<<ownedElementGlobalIDs_[i]<<" ";
  }

  oss<<"\n"<<ss.str()<<" Owned Global Nodes    = ";
  for(size_t i=0; i<ownedNodeGlobalIDs_.extent(0); i++) {
    oss<<ownedNodeGlobalIDs_[i]<<" ";
  }

  oss<<"\n"<<ss.str()<<" Owned Element2Node    = ";
  for(size_t i=0; i<ownedElementToNode_.extent(0); i++) {
    oss<<"(";
    for(size_t j=0; j<ownedElementToNode_.extent(1); j++) {
      oss<<ownedElementToNode_(i,j)<<" ";
    }
    oss<<") ";
  }

  oss<<"\n"<<ss.str()<<" Ghost Global Elements = ";
  for(size_t i=0; i<ghostElementGlobalIDs_.extent(0); i++) {
    oss<<ghostElementGlobalIDs_[i]<<" ";
  }
  oss<<"\n"<<ss.str()<<" Ghost Global Nodes    = ";
  for(size_t i=0; i<ghostNodeGlobalIDs_.extent(0); i++) {
    oss<<ghostNodeGlobalIDs_[i]<<" ";
  }

  oss<<"\n"<<ss.str()<<" Ghost Element2Node    = ";
  for(size_t i=0; i<ghostElementToNode_.extent(0); i++) {
    oss<<"(";
    for(size_t j=0; j<ghostElementToNode_.extent(1); j++) {
      oss<<ghostElementToNode_(i,j)<<" ";
    }
    oss<<") ";
  }

  oss << "\n"<<ss.str()<<" Owned And Ghost Nodes = ";
  for(size_t i=0; i<ownedAndGhostNodeGlobalIDs_.extent(0); i++) {
    oss << ownedAndGhostNodeGlobalIDs_[i]<<" ";
  }

  oss << "\n"<<ss.str()<<" Owned And Ghost Elements = ";
  for(size_t i=0; i<ownedAndGhostElementGlobalIDs_.extent(0); i++) {
    oss << ownedAndGhostElementGlobalIDs_[i]<<" ";
  }

  oss<<std::endl;
}


} // end of namespace TpetraExamples

#endif // TPETRAEXAMPLES_FEM_ASSEMBLY_MESH_DATABASE


