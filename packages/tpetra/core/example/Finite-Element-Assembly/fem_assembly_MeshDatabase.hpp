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


namespace TpetraExamples {

struct LLA
{
  LLA() {data[0]=0; data[1]=0;}

  LLA(global_ordinal_type i, global_ordinal_type j) {data[0]=i; data[1]=j;}

  global_ordinal_type operator[](global_ordinal_type i) const {return data[i];}
  global_ordinal_type & operator[](global_ordinal_type i) {return data[i];}

  global_ordinal_type data[2];
};

class MeshDatabase {
public:
  MeshDatabase(Teuchos::RCP<const Teuchos::Comm<int> > comm,
               global_ordinal_type global_elements_x,
               global_ordinal_type global_elements_y,
               global_ordinal_type procs_x,
               global_ordinal_type procs_y);

  ~MeshDatabase(){}

  // Size accessors
  size_t getNumOwnedElements() const {return ownedElementGlobalIDs_.extent(0);}

  size_t getNumGhostElements() const {return ghostElementGlobalIDs_.extent(0);}

  size_t getNumOwnedNodes() const {return ownedNodeGlobalIDs_.extent(0);}

  size_t getNumGhostNodes() const {return ghostNodeGlobalIDs_.extent(0);}

  size_t getNumOwnedAndGhostNodes() const {return ownedAndGhostNodeGlobalIDs_.extent(0);}

  size_t getNumOwnedAndGhostElements() const {return ownedAndGhostElementGlobalIDs_.extent(0);}

  // Data accessors
  global_ordinal_view_type getOwnedElementGlobalIDs() {return ownedElementGlobalIDs_;}
  global_ordinal_view_type getGhostElementGlobalIDs() {return ghostElementGlobalIDs_;}
  global_ordinal_view_type getOwnedAndGhostElementGlobalIDs() {return ownedAndGhostElementGlobalIDs_;}

  global_ordinal_view_type getOwnedNodeGlobalIDs() {return ownedNodeGlobalIDs_;}
  global_ordinal_view_type getGhostNodeGlobalIDs() {return ghostNodeGlobalIDs_;}
  global_ordinal_view_type getOwnedAndGhostNodeGlobalIDs() {return ownedAndGhostNodeGlobalIDs_;}

  global_ordinal_2d_array_type getOwnedElementToNode() {return ownedElementToNode_;}
  global_ordinal_2d_array_type getGhostElementToNode() {return ghostElementToNode_;}

  // Debugging output
  void print(std::ostream & oss);

  inline bool nodeIsOwned(global_ordinal_type idx) const{
    global_ordinal_type i,j;
    ij_from_idx(globalNodes_[0],idx,i,j);
    return nodeIsOwned(i,j);
  }

  inline bool elementIsOwned(global_ordinal_type idx) const{
    global_ordinal_type i,j;
    ij_from_idx(globalElements_[0], idx, i, j);
    return elementIsOwned(i,j);
  }


private:

  inline bool nodeIsOwned(global_ordinal_type i, global_ordinal_type j) const{
    return myNodeStart_[0] <= i &&  i < myNodeStop_[0] && myNodeStart_[1] <= j &&  j < myNodeStop_[1];
  }

  inline  bool elementIsOwned(global_ordinal_type i, global_ordinal_type j) const{
    return myElementStart_[0] <= i && i < myElementStop_[0] && myElementStart_[1] <= j && myElementStop_[1];
  }

  inline global_ordinal_type idx_from_ij(global_ordinal_type num_x, global_ordinal_type i, global_ordinal_type j) const {
    return j*num_x+i;
  }

  inline void ij_from_idx(global_ordinal_type num_x, global_ordinal_type idx, global_ordinal_type &i, global_ordinal_type&j) const {
    i = idx%num_x;
    j = (global_ordinal_type)((idx-i)/num_x);
  }

  void initializeOwnedAndGhostNodeGlobalIDs(void);

  void initializeOwnedAndGhostElementGlobalIDs(void);

  //wrapped dual views 
  global_ordinal_view_type ownedElementGlobalIDs_;
  global_ordinal_view_type ghostElementGlobalIDs_;

  global_ordinal_view_type ownedNodeGlobalIDs_;
  global_ordinal_view_type ghostNodeGlobalIDs_;

  global_ordinal_view_type ownedAndGhostNodeGlobalIDs_;
  global_ordinal_view_type ownedAndGhostElementGlobalIDs_;

  global_ordinal_2d_array_type ownedElementToNode_;
  global_ordinal_2d_array_type ghostElementToNode_;

  // Global Mesh Info
  LLA globalElements_;
  LLA globalNodes_;
  LLA globalProcs_;

  // Local Mesh Info
  global_ordinal_type myElementStart_[2];
  global_ordinal_type myElementStop_[2];
  global_ordinal_type myNodeStart_[2];
  global_ordinal_type myNodeStop_[2];

  // Comm info
  Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  size_t MyRank_;
  LLA myProcIJ_;
};



MeshDatabase::MeshDatabase(Teuchos::RCP<const Teuchos::Comm<int> > comm,
                           global_ordinal_type global_elements_x,
                           global_ordinal_type global_elements_y,
                           global_ordinal_type procs_x,
                           global_ordinal_type procs_y):
      globalElements_(global_elements_x,global_elements_y),
      globalNodes_(global_elements_x+1,global_elements_y+1),
      globalProcs_(procs_x,procs_y),
      comm_(comm)
{

  // NOTE: Elements/nodes are numbered sequentially with x as the "fast" direction
  // NOTE: assembly is all on host, so the overall scopeguard is sufficient here
  // Get processor decomp information
  MyRank_ = comm_->getRank();
  ij_from_idx(globalProcs_[0],MyRank_,myProcIJ_[0],myProcIJ_[1]);

  // Get local element & node start / stop
  global_ordinal_type num_my_elements=1, num_my_nodes=1;
  for(int k=0; k<2; k++) {
    global_ordinal_type eper = globalElements_[k] / globalProcs_[k];

    myElementStart_[k] = myProcIJ_[k] * eper;
    myElementStop_[k]  = (myProcIJ_[k] == globalProcs_[k]-1) ? globalElements_[k] : (myProcIJ_[k]+1) * eper;
    num_my_elements *= (myElementStop_[k]-myElementStart_[k]);

    myNodeStart_[k] = myProcIJ_[k] * eper;
    myNodeStop_[k]  = (myProcIJ_[k] == globalProcs_[k]-1) ? globalNodes_[k] : (myProcIJ_[k]+1) * eper;
    num_my_nodes *= (myNodeStop_[k]-myNodeStart_[k]);
  }

  // Generate the owned element ids
  auto ownedElementGlobalIDs = ownedElementGlobalIDs_.getHostView(Tpetra::Access::ReadWrite);
  Kokkos::resize(ownedElementGlobalIDs,num_my_elements);
  int ect=0;
  for(global_ordinal_type j=myElementStart_[1]; j<myElementStop_[1]; j++) {
    for(global_ordinal_type i=myElementStart_[0]; i<myElementStop_[0]; i++) {
      global_ordinal_type idx=idx_from_ij(globalElements_[0],i,j);
      ownedElementGlobalIDs(ect) = idx;
      ect++;
    }
  }

  // Generate the owned node ids
  auto _ownedNodeGlobalIDs = ownedNodeGlobalIDs_.getHostView(Tpetra::Access::ReadWrite);
  Kokkos::resize(_ownedNodeGlobalIDs,num_my_nodes);
  int nct=0;
  for(global_ordinal_type j=myNodeStart_[1]; j<myNodeStop_[1]; j++) {
    for(global_ordinal_type i=myNodeStart_[0]; i<myNodeStop_[0]; i++) {
      global_ordinal_type idx=idx_from_ij(globalNodes_[0],i,j);
      _ownedNodeGlobalIDs(nct) = idx;
      nct++;
    }
  }

  // Generate the element-to-node map
  // NOTE: Hardwired to QUAD4's.  Nodes are ordered exodus-style (counter-clockwise) within an element
  auto _ownedElementToNode = ownedElementToNode_.getHostView(Tpetra::Access::ReadWrite);
  Kokkos::resize(_ownedElementToNode,num_my_elements);
  int cct=0;
  for(global_ordinal_type j=myElementStart_[1]; j<myElementStop_[1]; j++) {
    for(global_ordinal_type i=myElementStart_[0]; i<myElementStop_[0]; i++) {
      // The (i,j) of the bottom left corner matches for elements & nodes
      global_ordinal_type nidx=idx_from_ij(globalNodes_[0],i,j);

      _ownedElementToNode(cct,0) = nidx;
      _ownedElementToNode(cct,1) = nidx+1;
      _ownedElementToNode(cct,2) = nidx+globalNodes_[0]+1;
      _ownedElementToNode(cct,3) = nidx+globalNodes_[0];
      cct++;
    }
  }

  // Generate the list of "ghost" elements & ghostElement2NodeMap
  // NOTE: This only generates a halo for elements where I own at least one node.  On the x/y hi sides,
  // the highers element does not own all the nodes on that proc.  Ergo, no halo in that direction
  std::vector<global_ordinal_type> my_ghost_elements;
  for(global_ordinal_type j=myElementStart_[1]-1; j<myElementStop_[1]; j++) {
    if(j<0 || j>=globalElements_[1]) continue; // Ignore stuff off the mesh
    for(global_ordinal_type i=myElementStart_[0]-1; i<myElementStop_[0]; i++) {
      if(i<0 || i>=globalElements_[0]) continue; // Ignore stuff off the mesh

      // Ignore proc interior
      if( j>myElementStart_[1]-1 && j<myElementStop_[1] && i>myElementStart_[0]-1 && i<myElementStop_[0])
        continue;

      global_ordinal_type idx=idx_from_ij(globalElements_[0],i,j);
      my_ghost_elements.push_back(idx);
    }
  }

  // NOTE: This are not recorded in Aztec/Ifpack/ML ordering.  Because most apps don't do that.
  auto _ghostElementGlobalIDs = ghostElementGlobalIDs_.getHostView(Tpetra::Access::ReadWrite);
  auto _ghostElementToNode = ghostElementToNode_.getHostView(Tpetra::Access::ReadWrite);
  Kokkos::resize(_ghostElementGlobalIDs,my_ghost_elements.size());
  Kokkos::resize(_ghostElementToNode,my_ghost_elements.size());
  for(size_t k=0; k<my_ghost_elements.size(); k++) {
    global_ordinal_type i,j, eidx= my_ghost_elements[k];
    _ghostElementGlobalIDs(k) = eidx;
    ij_from_idx(globalElements_[0],eidx,i,j);

    // The (i,j) of the bottom left corner matches for elements & nodes
    global_ordinal_type nidx=idx_from_ij(globalNodes_[0],i,j);

    _ghostElementToNode(k,0) = nidx;
    _ghostElementToNode(k,1) = nidx+1;
    _ghostElementToNode(k,2) = nidx+globalNodes_[0]+1;
    _ghostElementToNode(k,3) = nidx+globalNodes_[0];
  }

  // Generate the list of "ghost" nodes (aka any node that exists on the ownedElement list that isn't owned
  std::set<global_ordinal_type> my_ghost_nodes;
  auto ownedElementToNodeView = ownedElementToNode_.getHostView(Tpetra::Access::ReadOnly);  
  for(size_t k=0; k<ownedElementToNodeView.extent(0); k++) {
    for(size_t l=0; l<ownedElementToNodeView.extent(1); l++) {
      global_ordinal_type nidx=ownedElementToNodeView(k,l);
      if(!nodeIsOwned(nidx)) {
        my_ghost_nodes.insert(nidx);
      }
    }
  }

  auto _ghostNodeGlobalIDs = ghostNodeGlobalIDs_.getHostView(Tpetra::Access::ReadWrite);
  Kokkos::resize(_ghostNodeGlobalIDs,my_ghost_nodes.size());
  for(auto k=my_ghost_nodes.begin(); k!=my_ghost_nodes.end(); k++) {
    size_t kk = std::distance(my_ghost_nodes.begin(),k);
    _ghostNodeGlobalIDs(kk) = *k;
  }

  initializeOwnedAndGhostNodeGlobalIDs();
  initializeOwnedAndGhostElementGlobalIDs();
}


void MeshDatabase::initializeOwnedAndGhostNodeGlobalIDs(void)
{
  size_t total_size = getNumOwnedNodes() + getNumGhostNodes();

  auto _ownedAndGhostNodeGlobalIDs = ownedAndGhostNodeGlobalIDs_.getHostView(Tpetra::Access::ReadWrite);
  Kokkos::resize(_ownedAndGhostNodeGlobalIDs, total_size);

  {
    size_t insert_idx = 0;
    auto ownedNodeGlobalIDs = getOwnedNodeGlobalIDs().getHostView(Tpetra::Access::ReadOnly);
    auto ghostNodeGlobalIDs = getGhostNodeGlobalIDs().getHostView(Tpetra::Access::ReadOnly);
    for(size_t idx=0; idx < getNumOwnedNodes(); idx++)
    {
      _ownedAndGhostNodeGlobalIDs(insert_idx++) = ownedNodeGlobalIDs(idx);
    }
    for(size_t idx=0; idx < getNumGhostNodes(); idx++)
    {
      _ownedAndGhostNodeGlobalIDs(insert_idx++) = ghostNodeGlobalIDs(idx);
    }
  }
}


void MeshDatabase::initializeOwnedAndGhostElementGlobalIDs(void)
{
  size_t total_size = getNumOwnedElements() + getNumGhostElements();
  auto _ownedAndGhostElementGlobalIDs = ownedAndGhostElementGlobalIDs_.getHostView(Tpetra::Access::ReadWrite);
  Kokkos::resize(_ownedAndGhostElementGlobalIDs, total_size);

  {
    size_t insert_idx = 0;
    auto ownedElementGlobalIDs = getOwnedElementGlobalIDs().getHostView(Tpetra::Access::ReadOnly);
    auto ghostElementGlobalIDs = getGhostElementGlobalIDs().getHostView(Tpetra::Access::ReadOnly);
    for(size_t idx=0; idx<getNumOwnedElements(); idx++)
    {
      _ownedAndGhostElementGlobalIDs(insert_idx++) = ownedElementGlobalIDs(idx);
    }
    for(size_t idx=0; idx<getNumGhostElements(); idx++)
    {
      _ownedAndGhostElementGlobalIDs(insert_idx++) = ghostElementGlobalIDs(idx);
    }
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
  {
    auto IDs = ownedElementGlobalIDs_.getHostView(Tpetra::Access::ReadOnly);
    for(size_t i=0; i<IDs.extent(0); i++) {
      oss<<IDs[i]<<" ";
    }
  }

  oss<<"\n"<<ss.str()<<" Owned Global Nodes    = ";
  {
    auto IDs = ownedNodeGlobalIDs_.getHostView(Tpetra::Access::ReadOnly);
    for(size_t i=0; i<IDs.extent(0); i++) {
      oss<<IDs[i]<<" ";
    }
  }

  oss<<"\n"<<ss.str()<<" Owned Element2Node    = ";
  {
    auto IDs = ownedElementToNode_.getHostView(Tpetra::Access::ReadOnly);
    for(size_t i=0; i<IDs.extent(0); i++) {
      oss<<"(";
      for(size_t j=0; j<IDs.extent(1); j++) {
	oss<<IDs(i,j)<<" ";
      }
      oss<<") ";
    }
  }

  oss<<"\n"<<ss.str()<<" Ghost Global Elements = ";
  {
    auto IDs = ghostElementGlobalIDs_.getHostView(Tpetra::Access::ReadOnly);
    for(size_t i=0; i<IDs.extent(0); i++) {
      oss<<IDs[i]<<" ";
    }
  }

  oss<<"\n"<<ss.str()<<" Ghost Global Nodes    = ";
  {
    auto IDs = ghostNodeGlobalIDs_.getHostView(Tpetra::Access::ReadOnly);
    for(size_t i=0; i<IDs.extent(0); i++) {
      oss<<IDs[i]<<" ";
    }
  }

  oss<<"\n"<<ss.str()<<" Ghost Element2Node    = ";
  {
    auto IDs = ghostElementToNode_.getHostView(Tpetra::Access::ReadOnly);
    for(size_t i=0; i<IDs.extent(0); i++) {
      oss<<"(";
      for(size_t j=0; j<IDs.extent(1); j++) {
	oss<<IDs(i,j)<<" ";
      }
      oss<<") ";
    }
  }

  oss << "\n"<<ss.str()<<" Owned And Ghost Nodes = ";
  {
    auto IDs = ownedAndGhostNodeGlobalIDs_.getHostView(Tpetra::Access::ReadOnly);
    for(size_t i=0; i<IDs.extent(0); i++) {
      oss << IDs[i]<<" ";
    }
  }

  oss << "\n"<<ss.str()<<" Owned And Ghost Elements = ";
  {
    auto IDs = ownedAndGhostElementGlobalIDs_.getHostView(Tpetra::Access::ReadOnly);
    for(size_t i=0; i<IDs.extent(0); i++) {
      oss << IDs[i]<<" ";
    }
  }

  oss<<std::endl;
}


// A toy class designed to simulate the cost of ghosting element state
// information.  This class does no actual computation, but moves an
// amount of data similar to what would happen from ghosting element
// state information.

class GhostState {
public:
  GhostState (Teuchos::RCP<const import_type> importer,
              const int amountOfStatePerElement) :
    elementImporter_ (importer)
  {
    using Teuchos::rcp;
    using MV = multivector_type;
    if (! elementImporter_.is_null ()) {
      ownedState_ = rcp (new MV (importer->getSourceMap (),
                                 amountOfStatePerElement));
      ghostState_ = rcp (new MV (importer->getTargetMap (),
                                 amountOfStatePerElement));
    }
  }

  void doGhost () {
    if (! elementImporter_.is_null ()) {
      ghostState_->doImport(*ownedState_,*elementImporter_,Tpetra::REPLACE);
    }
  }

private:
  Teuchos::RCP<const import_type> elementImporter_;
  Teuchos::RCP<multivector_type> ownedState_;
  Teuchos::RCP<multivector_type> ghostState_;
};

} // namespace TpetraExamples

#endif // TPETRAEXAMPLES_FEM_ASSEMBLY_MESH_DATABASE


