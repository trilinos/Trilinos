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

#include "typedefs.hpp"
#include "Teuchos_Comm.hpp"

class MeshDatabase {
public:
  MeshDatabase(Teuchos::RCP<const Teuchos::Comm<int> > comm, int global_elements_x, int global_elements_y):comm_(comm){

    int rank     = comm->getRank();
    int numProcs = comm->getSize();

    

  }

  ~MeshDatabase(){}

  // Size accessors
  size_t getNumOwnedElements() const {return ownedElementGlobalIDs_.dimension(0);}

  size_t getNumGhostElements() const {return ghostElementGlobalIDs_.dimension(0);}

  size_t getNumOwnedNodes() const {return ownedNodeGlobalIDs_.dimension(0);}

  size_t getNumGhostNodes() const {return ghostNodeGlobalIDs_.dimension(0);}
  

  // Data accessors
  global_ordinal_view_type getOwnedElementGlobalIDs() {return ownedElementGlobalIDs_;}
  global_ordinal_view_type getGhostGlobalIDs() {return ghostElementGlobalIDs_;}

  global_ordinal_view_type getOwnedElementNodeIDs() {return ownedNodeGlobalIDs_;}
  global_ordinal_view_type getGhostNodeIDs() {return ghostNodeGlobalIDs_;}

  local_ordinal_2d_array_type getOwnedElementToNode() {return ownedElementToNode_;}
  local_ordinal_2d_array_type getGhostElementToNode() {return ghostElementToNode_;}
  

private:
  global_ordinal_view_type ownedElementGlobalIDs_;
  global_ordinal_view_type ghostElementGlobalIDs_;
  global_ordinal_view_type ownedNodeGlobalIDs_;
  global_ordinal_view_type ghostNodeGlobalIDs_;

  local_ordinal_2d_array_type ownedElementToNode_;  
  local_ordinal_2d_array_type ghostElementToNode_;

  Teuchos::RCP<const Teuchos::Comm<int> > comm_;
};


// Generates a dummy finite element stiffness matrix for quads
//scalar_2d_array_type generateFiniteElementMatrix() {



//}

#endif
