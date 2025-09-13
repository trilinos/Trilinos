// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_ParameterList.hpp"

#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Selector.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"

TEUCHOS_UNIT_TEST(GenerateInitialExodusMesh, WriteFunction)
{
  using Teuchos::RCP;

  // To run in parallel we need to set an ioss property to split
  // genesis file across mpi ranks.
  // setenv("IOSS_PROPERTIES", "DECOMPOSITION_METHOD=rib", 1);

  // ********************************
  // Put analytic field on exodus file and write exodus file
  // ********************************
  {
    std::string file_name = "quad8.gen";
    panzer_stk::SquareQuadMeshFactory factory;
    auto params = Teuchos::parameterList("Mesh Parameters");
    params->set("X0",-1.0);
    params->set("Y0",-1.0);
    params->set("Xf",1.0);
    params->set("Yf",1.0);
    params->set("X Blocks",1);
    params->set("Y Blocks",1);
    params->set("X Elements",10);
    params->set("Y Elements",10);
    factory.setParameterList(params);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildUncommitedMesh(MPI_COMM_WORLD);
    
    std::vector<std::string> block_names;
    mesh->getElementBlockNames(block_names);
    for (const auto& block : block_names)
      mesh->addCellField("TEST_FIELD",block);

    const bool setupIO = true;
    const bool doPerceptRefinement = false;
    mesh->initialize(MPI_COMM_WORLD,setupIO,doPerceptRefinement);
    factory.completeMeshConstruction(*mesh,MPI_COMM_WORLD);

    const int num_time_steps = 10;
    const double dt = 0.1;
    double time = 0.0;
    const double source_multiplier = 5.0;
    bool first = true;
    
    for (int time_index = 0; time_index < num_time_steps+1; ++time_index) {
      for (const auto& block : block_names) {
        std::vector<stk::mesh::Entity> elements;
        mesh->getMyElements(block,elements);
        auto field = mesh->getCellField("TEST_FIELD",block);
        for (size_t elem_idx=0; elem_idx < elements.size(); ++elem_idx) {

          // Compute element centroid
          double x = 0.0;
          double y = 0.0;
          for (int node_idx=0; node_idx < 4; ++node_idx) {
            auto node = mesh->findConnectivityById(elements[elem_idx],mesh->getNodeRank(),node_idx);
            const double * const coords = mesh->getNodeCoordinates(node);
            x += coords[0];
            y += coords[1];
          }
          x /= 4.0;
          y /= 4.0;
          double tmp = x*x + y*y;
          double r = (tmp > 0) ? std::sqrt(tmp) : 0.0;

          double* value = static_cast<double*>(stk::mesh::field_data(*field,elements[elem_idx]));
          if ( (time > 0.3) && (r < 0.5) ) {
            *value = time * source_multiplier + 4.0 * r;
          }
          else {
            *value = 0.0;
          }

          // std::cout  << "block=" << block
          //            << ", element=" << mesh->getBulkData()->identifier(elements[elem_idx])
          //            << ", x=" << x
          //            << ", y=" << y
          //            << ", val=" << *value << std::endl;
        }
      }
      if (first) {
        mesh->setInitialStateTime(0.0);
        mesh->setupExodusFile("energy-load-cell-data-SOURCE.exo");
        first = false;
      }
      mesh->writeToExodus(time);
      time += dt;
    }

  }
}
