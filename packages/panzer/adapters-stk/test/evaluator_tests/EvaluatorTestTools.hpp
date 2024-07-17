// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EVALUATOR_TEST_TOOLS
#define EVALUATOR_TEST_TOOLS

#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_DOFManager.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_STK_WorksetFactory.hpp"

#include "Panzer_OrientationsInterface.hpp"

namespace EvaluatorTestTools {

struct WorksetsAndOrts {

  Teuchos::RCP<std::vector<panzer::Workset> > worksets;
  Teuchos::RCP<const std::vector<Intrepid2::Orientation> > orientations;

};

/* @brief Create worksets and orientations for purposes of evaluator testing.
*
* @param[in] mesh STK Mesh
* @param[in] fmap Map of field names to basis descriptors
* @param[in] workset_size Maximal workset size. Use \c -1 to use one workset (all elements).
*
* @warning Will likely need to be modified for future evaluator tests.
*/

WorksetsAndOrts getWorksetsAndOrtsForFields(
  Teuchos::RCP<panzer_stk::STK_Interface> mesh, 
  std::map<std::string,panzer::BasisDescriptor> fmap,
  const int workset_size)
{

  WorksetsAndOrts wksOrts;

  // Dof manager is, at least, a global indexer
  // and is passed in/out, if needed

  // And a needs map (block to needs)
  std::map<std::string, panzer::WorksetNeeds> needs_map;

  // Assuming one element block for now...
  std::vector<std::string> eblocks;
  mesh->getElementBlockNames(eblocks);

  Teuchos::RCP<panzer::DOFManager> dof_manager;
  {

    // Build a connectivity manager for the given mesh
    const auto conn_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

    // Initialize the dof manager with the conn manager
    dof_manager = Teuchos::rcp(new panzer::DOFManager(conn_manager,MPI_COMM_WORLD));

    // Build basis and fields
    const auto & cell_topology = mesh->getCellTopology(eblocks[0]);
    // Set cell data
    // numCells will get overwritten with our particular path thru getWorksets below
    needs_map[eblocks[0]].cellData = panzer::CellData(1,cell_topology);

    for (auto & map : fmap)
    {
      auto name = map.first;
      auto bd = map.second;

      auto intrepid_basis = panzer::createIntrepid2Basis<PHX::Device::execution_space,double,double>(bd.getType(),bd.getOrder(),cell_topology);
      auto field_pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(intrepid_basis));

      // Add field to dof manager to lock in the field pattern
      dof_manager->addField(eblocks[0], name, field_pattern);

      needs_map[eblocks[0]].addBasis(bd);
    }
  }

  // Finalize the dof manager
  dof_manager->buildGlobalUnknowns();

  // Build workset factory
  auto factory = Teuchos::rcp(new panzer_stk::WorksetFactory(mesh));

  auto workset_container = Teuchos::rcp(new panzer::WorksetContainer(factory, needs_map));
  workset_container->setGlobalIndexer(dof_manager);
  wksOrts.worksets = workset_container->getWorksets(panzer::WorksetDescriptor(eblocks[0],
                                                    workset_size, false, true));

  wksOrts.orientations = workset_container->getOrientations();

  return wksOrts;

}

/* @brief Create a basic inline mesh for testing purposes 
* 
* @param[in] pl Mesh parameter list
*
* @note Limited to Square, Cube geometry and Quad, Tri, Tet, Hex 
*/

Teuchos::RCP<panzer_stk::STK_Interface> createInlineMesh(Teuchos::RCP<Teuchos::ParameterList> pl)
{

  auto dim = pl->get<int>("Mesh Dimension");
  auto type = pl->get<std::string>("Type");

  Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;

  if (type == "Quad" && dim == 2) {
    mesh_factory = Teuchos::rcp(new panzer_stk::SquareQuadMeshFactory);
  } else if (type == "Tri" && dim == 2) {
    mesh_factory = Teuchos::rcp(new panzer_stk::SquareTriMeshFactory);
  } else if (type == "Tet" && dim == 3) {
    mesh_factory = Teuchos::rcp(new panzer_stk::CubeTetMeshFactory);
  } else if (type == "Hex" && dim == 3) {
    mesh_factory = Teuchos::rcp(new panzer_stk::CubeHexMeshFactory);
  } else {
    // Throw an error
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
      "ERROR: Type and/or dimension of inline mesh not valid!");
  }

  mesh_factory->setParameterList(Teuchos::rcp(new Teuchos::ParameterList(pl->sublist("Mesh Factory Parameter List"))));
  auto mesh = mesh_factory->buildMesh(MPI_COMM_WORLD);

  return mesh;
}

}

#endif