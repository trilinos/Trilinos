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

#ifndef __Panzer_STK_Quad8ToQuad4MeshFactory_hpp__
#define __Panzer_STK_Quad8ToQuad4MeshFactory_hpp__

#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>
#include <vector>
#include <string>

namespace panzer_stk {

class STK_Interface;

/** This class reads in a Quad8 mesh and converts it to a Quad4
    mesh. It will create the cells/nodes/edges, copy the sideset and
    nodeset data, and copy a list of user provided cell (not nodal)
    field data into the mesh.
  */
class Quad8ToQuad4MeshFactory : public STK_MeshFactory {
public:

  Quad8ToQuad4MeshFactory(const std::string& quad8MeshFileName,
                          stk::ParallelMachine mpi_comm = MPI_COMM_WORLD,
                          const bool print_debug = false);

  Quad8ToQuad4MeshFactory(const Teuchos::RCP<panzer_stk::STK_Interface>& quad8Mesh,
                          const bool print_debug = false);

  Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const;
  virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const;
  virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const;

  //! Derived from ParameterListAcceptor
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

  //! Derived from ParameterListAcceptor
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

protected:

  void buildMetaData(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
  void buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;

  void addSideSets(STK_Interface & mesh) const;
  void addNodeSets(STK_Interface & mesh) const;
  void addEdgeBlocks(STK_Interface & mesh) const;
  void copyCellFieldData(STK_Interface & mesh) const;

  Teuchos::RCP<panzer_stk::STK_Interface> quad8Mesh_;

  mutable unsigned int machRank_, machSize_;

  bool createEdgeBlocks_;

  /// If true, offset mesh GIDs to exercise 32-bit limits.
  bool offsetGIDs_;
  mutable stk::mesh::EntityId offset_;

  bool print_debug_;

  std::string edgeBlockName_;
};

}

#endif
