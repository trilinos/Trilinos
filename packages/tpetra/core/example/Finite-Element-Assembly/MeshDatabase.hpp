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


#ifndef MESHDATABASE_HPP
#define MESHDATABASE_HPP
#include <iostream>
#include <sstream>

#include "typedefs.hpp"
#include "Teuchos_Comm.hpp"

struct LLA {
  LLA() {data[0]=0; data[1]=0;}

  LLA(GlobalOrdinal i, GlobalOrdinal j) {data[0]=i; data[1]=j;}

  GlobalOrdinal operator[](GlobalOrdinal i) const {return data[i];}
  GlobalOrdinal & operator[](GlobalOrdinal i){return data[i];}

  GlobalOrdinal data[2];
};


class MeshDatabase {
public:
  MeshDatabase(Teuchos::RCP<const Teuchos::Comm<int> > comm, int global_elements_x, int global_elements_y, int procs_x, int procs_y):
    globalElements_(global_elements_x,global_elements_y),
    globalNodes_(global_elements_x+1,global_elements_y+1),
    globalProcs_(procs_x,procs_y),
    comm_(comm) {

    // NOTE: Elements/nodes are numbered sequentially with x as the "fast" direction

    // Get processor decomp information
    MyRank_ = comm_->getRank();
    ij_from_idx(globalProcs_[0],MyRank_,myProcIJ_[0],myProcIJ_[1]);

    // Get local element & node start / stop
    GlobalOrdinal num_my_elements=1, num_my_nodes=1;
    for(int k=0; k<2; k++) {
      GlobalOrdinal eper = globalElements_[k] / globalProcs_[k];

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
    for(GlobalOrdinal j=myElementStart_[1]; j<myElementStop_[1]; j++) {
      for(GlobalOrdinal i=myElementStart_[0]; i<myElementStop_[0]; i++) {
	GlobalOrdinal idx=idx_from_ij(globalElements_[0],i,j);
	ownedElementGlobalIDs_[ect] = idx;
	ect++;	
      }
    }
   
    // Generate the owned node ids
    Kokkos::resize(ownedNodeGlobalIDs_,num_my_nodes);
    int nct=0;
    for(GlobalOrdinal j=myNodeStart_[1]; j<myNodeStop_[1]; j++) {
      for(GlobalOrdinal i=myNodeStart_[0]; i<myNodeStop_[0]; i++) {
	GlobalOrdinal idx=idx_from_ij(globalNodes_[0],i,j);
	ownedNodeGlobalIDs_[nct] = idx;
	nct++;	
      }
    }
   
    // Generate the element-to-node map
    // NOTE: Hardwired to QUAD4's.  Nodes are ordered exodus-style (counter-clockwise) within an element
    Kokkos::resize(ownedElementToNode_,num_my_elements,4);
    int cct=0;
    for(GlobalOrdinal j=myElementStart_[1]; j<myElementStop_[1]; j++) {
      for(GlobalOrdinal i=myElementStart_[0]; i<myElementStop_[0]; i++) {
	// The (i,j) of the bottom left corner matches for elements & nodes
	GlobalOrdinal idx=idx_from_ij(globalNodes_[0],i,j);

	ownedElementToNode_(cct,0) = idx;
	ownedElementToNode_(cct,1) = idx+1;
	ownedElementToNode_(cct,2) = idx+globalNodes_[0]+1;
	ownedElementToNode_(cct,3) = idx+globalNodes_[0];
	cct++;
      }
    }

  }

  ~MeshDatabase(){}

  // Size accessors
  size_t getNumOwnedElements() const {return ownedElementGlobalIDs_.dimension(0);}

  size_t getNumGhostElements() const {return ghostElementGlobalIDs_.dimension(0);}

  size_t getNumOwnedNodes() const {return ownedNodeGlobalIDs_.dimension(0);}

  size_t getNumGhostNodes() const {return ghostNodeGlobalIDs_.dimension(0);}
  

  // Data accessors
  global_ordinal_view_type getOwnedElementGlobalIDs() {return ownedElementGlobalIDs_;}
  global_ordinal_view_type getGhostElementGlobalIDs() {return ghostElementGlobalIDs_;}

  global_ordinal_view_type getOwnedNodeGlobalIDs() {return ownedNodeGlobalIDs_;}
  global_ordinal_view_type getGhostNodeGloballDs() {return ghostNodeGlobalIDs_;}

  local_ordinal_2d_array_type getOwnedElementToNode() {return ownedElementToNode_;}
  local_ordinal_2d_array_type getGhostElementToNode() {return ghostElementToNode_;}
  

  void print(std::ostream & oss) { 
    std::ostringstream ss;
    ss<<"["<<MyRank_<<","<<myProcIJ_[0]<<","<<myProcIJ_[1]<<"]";
    oss<<ss.str()<<" Global Elements = ["<<globalElements_[0]<<"x"<<globalElements_[1]<<"] Nodes ="<<globalNodes_[0]<<"x"<<globalNodes_[1]<<"]\n";
    oss<<ss.str()<<" Stop/Start Elements = ["<<myElementStart_[0]<<","<<myElementStop_[0]<<")x["<<myElementStart_[1]<<","<<myElementStop_[1]<<")\n";
    oss<<ss.str()<<" Stop/Start Nodes    = ["<<myNodeStart_[0]<<","<<myNodeStop_[0]<<")x["<<myNodeStart_[1]<<","<<myNodeStop_[1]<<")\n";

    oss<<ss.str()<<" Owned Global Elements = ";
    for(size_t i=0; i<ownedElementGlobalIDs_.dimension(0); i++)
      oss<<ownedElementGlobalIDs_[i]<<" ";

    oss<<"\n"<<ss.str()<<" Owned Global Nodes   = ";
    for(size_t i=0; i<ownedNodeGlobalIDs_.dimension(0); i++)
      oss<<ownedNodeGlobalIDs_[i]<<" ";

    oss<<"\n"<<ss.str()<<" Owned Element2Node   = ";
    for(size_t i=0; i<ownedElementToNode_.dimension(0); i++) {
      oss<<"(";
      for(size_t j=0; j<ownedElementToNode_.dimension(1); j++)
	oss<<ownedElementToNode_(i,j)<<" ";
      oss<<") ";
    }      

    oss<<std::endl;
  }

private:

  inline GlobalOrdinal idx_from_ij(int num_x, GlobalOrdinal i, GlobalOrdinal j) const {return j*num_x+i;}
  inline void ij_from_idx(int num_x, GlobalOrdinal idx, GlobalOrdinal &i, GlobalOrdinal&j) const { 
    i = idx%num_x;
    j = (GlobalOrdinal)((idx-i)/num_x);
  }

  global_ordinal_view_type ownedElementGlobalIDs_;
  global_ordinal_view_type ghostElementGlobalIDs_;
  global_ordinal_view_type ownedNodeGlobalIDs_;
  global_ordinal_view_type ghostNodeGlobalIDs_;

  local_ordinal_2d_array_type ownedElementToNode_;  
  local_ordinal_2d_array_type ghostElementToNode_;

  // Global Mesh Info
  LLA globalElements_;
  LLA globalNodes_;
  LLA globalProcs_;

  // Local Mesh Info
  GlobalOrdinal myElementStart_[2];
  GlobalOrdinal myElementStop_[2];
  GlobalOrdinal myNodeStart_[2];
  GlobalOrdinal myNodeStop_[2];

  // Comm info
  Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  int MyRank_;
  LLA myProcIJ_;


  };


// Generates a dummy finite element stiffness matrix for quads
//scalar_2d_array_type generateFiniteElementMatrix() {



//}

#endif
