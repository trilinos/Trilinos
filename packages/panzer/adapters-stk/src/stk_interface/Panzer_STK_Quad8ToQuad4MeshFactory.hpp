// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
                          stk::ParallelMachine mpi_comm = MPI_COMM_WORLD, // CHECK: ALLOW MPI_COMM_WORLD
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
