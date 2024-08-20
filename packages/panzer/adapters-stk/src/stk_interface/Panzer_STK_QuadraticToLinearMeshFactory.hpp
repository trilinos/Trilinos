// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_QuadraticToLinearMeshFactory_hpp__
#define __Panzer_STK_QuadraticToLinearMeshFactory_hpp__

#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>
#include <vector>
#include <string>

namespace panzer_stk {

class STK_Interface;

/** This class reads in a second-order (quadratic) mesh and converts it 
 *  to a first-order (linear) mesh. It will create the cells/nodes/edges, copy the sideset and
 *  nodeset data, and copy a list of user provided cell (not nodal) field data into the mesh.
 */
class QuadraticToLinearMeshFactory : public STK_MeshFactory {
public:

  QuadraticToLinearMeshFactory(const std::string& quadMeshFileName,
                          stk::ParallelMachine mpi_comm = MPI_COMM_WORLD, // CHECK: ALLOW MPI_COMM_WORLD
                          const bool print_debug = false);

  QuadraticToLinearMeshFactory(const Teuchos::RCP<panzer_stk::STK_Interface>& quadMesh,
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

  /**
   * Infer the output topology from the given input mesh.
   * Also, ensure the input mesh topology matches the parameter list.
   * Currently, each element block must have the same topology.
   */
  void getOutputTopology();

  Teuchos::RCP<panzer_stk::STK_Interface> quadMesh_; //! Second order mesh

  mutable unsigned int machRank_, machSize_;

  bool createEdgeBlocks_;

  /// If true, offset mesh GIDs to exercise 32-bit limits.
  bool offsetGIDs_;
  mutable stk::mesh::EntityId offset_;

  bool print_debug_;

  std::string edgeBlockName_;

  const CellTopologyData * outputTopoData_; //! Output mesh topology data 

  unsigned int nDim_; //! Dimension of the mesh
  unsigned int nNodes_; //! Nodes in one element of the linear basis

  //! List of currently supported input topologies.
  std::vector<shards::CellTopology> supportedInputTopos_ = {
    shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<8>>()),
    shards::CellTopology(shards::getCellTopologyData<shards::Triangle<6>>()),
    shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<10>>()),
    shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<20>>())
    };

  //! Map from input topology to the output shards topology data. The list here is
  //! currently supported. Right now, this is one-to-one and may need to be expanded.
  std::map<const std::string,const CellTopologyData *> outputTopoMap_ = {
    {shards::getCellTopologyData<shards::Quadrilateral<8>>()->name,
     shards::getCellTopologyData<shards::Quadrilateral<4>>()},
    {shards::getCellTopologyData<shards::Triangle<6>>()->name,
     shards::getCellTopologyData<shards::Triangle<3>>()},
    {shards::getCellTopologyData<shards::Tetrahedron<10>>()->name,
     shards::getCellTopologyData<shards::Tetrahedron<4>>()},
    {shards::getCellTopologyData<shards::Hexahedron<20>>()->name,
     shards::getCellTopologyData<shards::Hexahedron<8>>()}
    };
};

}

#endif