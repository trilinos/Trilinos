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
#include "Panzer_STK_ExodusReaderFactory.hpp"
#include "Panzer_STK_Quad8ToQuad4MeshFactory.hpp"

bool is_equal(const double& x, const double& x_gold) {
  constexpr double tol = 10.0*std::numeric_limits<double>::epsilon();
  return std::fabs(x-x_gold) < tol ? true : false;
}
bool in_range(const double& x, const double& left, const double& right) {
  TEUCHOS_ASSERT(left<right);
  constexpr double tol = 10.0*std::numeric_limits<double>::epsilon();
  if ( (x >= left-tol) && (x <= right+tol) )
    return true;

  return false;
}

TEUCHOS_UNIT_TEST(tQuad8ToQuad4Factory, periodic_input)
{
  using Teuchos::RCP;

  // To run in parallel we need to set an ioss property to split
  // genesis file across mpi ranks.
  setenv("IOSS_PROPERTIES", "DECOMPOSITION_METHOD=rib", 1);

  // ********************************
  // Put analytic field on Quad8 genesis file and write exodus file
  // ********************************
  {
    std::string file_name = "meshes/quad8.gen";
    panzer_stk::STK_ExodusReaderFactory factory(file_name);
    RCP<panzer_stk::STK_Interface> mesh = factory.buildMesh(MPI_COMM_WORLD);

    std::vector<std::string> block_names;
    mesh->getElementBlockNames(block_names);
    for (const auto& block : block_names)
      mesh->addCellField("TEST_FIELD",block);

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

        double* value = static_cast<double*>(stk::mesh::field_data(*field,elements[elem_idx]));
        *value = x + 4.0 * y * y;

        // std::cout  << "block=" << block
        //            << ", element=" << mesh->getBulkData()->identifier(elements[elem_idx])
        //            << ", x=" << x
        //            << ", y=" << y
        //            << ", val=" << *value << std::endl;
      }
    }
    mesh->setupExodusFile("quad8.exo");
    mesh->writeToExodus(0.0);
  }

  // ********************************
  // Convert the Quad8 to Quad4
  // ********************************
  std::string file_name = "quad8.exo";
  RCP<panzer_stk::STK_Interface> mesh;
  {
    // The quad8.exo file output above is in split form, not a single
    // file. need to disable the ioss flag above to read the split
    // file correctly.
    unsetenv("IOSS_PROPERTIES");

    const bool print_debug = false;
    panzer_stk::Quad8ToQuad4MeshFactory factory(file_name,MPI_COMM_WORLD,print_debug);
    mesh = factory.buildMesh(MPI_COMM_WORLD);

    // Write the mesh
    mesh->setupExodusFile("quad4.exo");
    mesh->writeToExodus(0.0);
  }

  // ********************************
  // Test the inline mesh field values against analytic function
  // ********************************
  {
    std::vector<std::string> block_names;
    mesh->getElementBlockNames(block_names);

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

        double* value = static_cast<double*>(stk::mesh::field_data(*field,elements[elem_idx]));
        double gold_value = x + 4.0 * y * y;

        double tol = 10.0*std::numeric_limits<double>::epsilon();
        TEST_FLOATING_EQUALITY(*value,gold_value,tol);

        // std::cout  << "block=" << block
        //            << ", element=" << mesh->getBulkData()->identifier(elements[elem_idx])
        //            << ", x=" << x
        //            << ", y=" << y
        //            << ", val=" << *value << std::endl;
      }
    }
  }

  // ********************************
  // Test that the sidesets are correct
  // ********************************

  // outer_box
  {
    std::vector<stk::mesh::Entity> sides;
    mesh->getMySides("outer_box",sides);

    int local_size = sides.size();
    int global_size = 0;
    MPI_Allreduce(&local_size,&global_size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    TEST_EQUALITY(global_size,100);

    for (const auto& side : sides) {
      for (int node=0; node < 2; ++node) {
        auto node_entity = mesh->findConnectivityById(side,mesh->getNodeRank(),node);
        const double * const coords = mesh->getNodeCoordinates(node_entity);
        const double x = coords[0];
        const double y = coords[1];
        const double x_bound = 0.145;
        const double y_bound = 0.1;
        // std::cout << "x=" << x << ", y=" << y << std::endl;
        TEST_ASSERT(
                    (is_equal(x, x_bound) && in_range(y,-y_bound,y_bound)) ||
                    (is_equal(x,-x_bound) && in_range(y,-y_bound,y_bound)) ||
                    (is_equal(y, y_bound) && in_range(x,-x_bound,x_bound)) ||
                    (is_equal(y,-y_bound) && in_range(x,-x_bound,x_bound))
                     );
      }
    }
  }

  // Inner box
  {
    std::vector<stk::mesh::Entity> sides;
    mesh->getMySides("inner_box",sides);

    int local_size = sides.size();
    int global_size = 0;
    MPI_Allreduce(&local_size,&global_size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    TEST_EQUALITY(global_size,84);

    for (const auto& side : sides) {
      for (int node=0; node < 2; ++node) {
        auto node_entity = mesh->findConnectivityById(side,mesh->getNodeRank(),node);
        const double * const coords = mesh->getNodeCoordinates(node_entity);
        const double x = coords[0];
        const double y = coords[1];
        const double x_bound = 0.125;
        const double y_bound = 0.08;
        // std::cout << "x=" << x << ", y=" << y << std::endl;
        TEST_ASSERT(
                    (is_equal(x, x_bound) && in_range(y,-y_bound,y_bound)) ||
                    (is_equal(x,-x_bound) && in_range(y,-y_bound,y_bound)) ||
                    (is_equal(y, y_bound) && in_range(x,-x_bound,x_bound)) ||
                    (is_equal(y,-y_bound) && in_range(x,-x_bound,x_bound))
                     );
      }
    }
  }

}
