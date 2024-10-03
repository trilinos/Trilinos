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
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Tpetra_Core.hpp"

#include "Teuchos_GlobalMPISession.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"

#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_STKConnManager.hpp"
#include "PanzerSTK_UnitTest_BuildMesh.hpp"

#include <string>

namespace panzer_stk {

namespace {

Teuchos::RCP<const panzer::LocalMeshInfo>
buildParallelLocalMeshInfo(const std::vector<int> & N,   // Cells per dimension
                           const std::vector<int> & B,   // Blocks per dimension
                           const std::vector<int> & P,   // Processors per dimension
                           const std::vector<double> &L, // Domain length per block
                           const std::vector<int> periodic_dims) // Periodic dimensions
{
  Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildParallelMesh(N,B,P,L,periodic_dims);
  return generateLocalMeshInfo(*mesh);
}

void
runConnectivityTest(const std::vector<int> & N,   // Cells per dimension
                    const std::vector<int> & B,   // Blocks per dimension
                    const std::vector<int> & P,   // Processors per dimension
                    const std::vector<double> &L, // Domain length per block
                    const std::vector<int> periodic_dims, // Periodic dimensions
                    Teuchos::FancyOStream &out,
                    bool &success)
{

  // The point of this test is to make sure that the LocalMeshInfo object connectivity
  // agrees with the ConnManager.
  // Note the ConnManager does not contain all of the connectivity found in the LocalMeshInfo (associated with periodic boundary conditions)

  Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildParallelMesh(N, B, P, L, periodic_dims);
  auto mesh_info = generateLocalMeshInfo(*mesh);

  auto conn = Teuchos::rcp(new panzer_stk::STKConnManager(mesh));

  std::vector<shards::CellTopology> topologies;
  conn->getElementBlockTopologies(topologies);
  const auto topology = topologies[0];

  Teuchos::RCP<panzer::FieldPattern> cell_pattern;
  if(topology.getDimension() == 1){
    cell_pattern = Teuchos::rcp(new panzer::EdgeFieldPattern(topology));
  } else if(topology.getDimension() == 2){
    cell_pattern = Teuchos::rcp(new panzer::FaceFieldPattern(topology));
  } else if(topology.getDimension() == 3){
    cell_pattern = Teuchos::rcp(new panzer::ElemFieldPattern(topology));
  }
  conn->buildConnectivity(*cell_pattern);

  std::vector<std::string> element_blocks;
  conn->getElementBlockIds(element_blocks);

  // TODO: "cell_sets" is an advanced feature that is not yet setup

//  // Owned cells
//  for(const auto & element_block : element_blocks){
//    const auto block_id = mesh_info->cell_sets->getSetIndex(element_block);
//    for(const auto & local_cell_id : conn->getElementBlock(element_block)){
//      const auto global_cell_id = conn->getConnectivity(local_cell_id)[0];
//
//      // Make sure the local mesh info thinks the same thing
//      TEST_EQUALITY(mesh_info->global_cells(local_cell_id), global_cell_id);
//    }
//  }
//
//  // Whatever ghost cells are in the ConnManager
//  for(const auto & element_block : element_blocks){
//    const auto block_id = mesh_info->cell_sets->getSetIndex(element_block);
//    for(const auto & local_cell_id : conn->getNeighborElementBlock(element_block)){
//      const auto global_cell_id = conn->getConnectivity(local_cell_id)[0];
//
//      // Make sure the local mesh info thinks the same thing
//      TEST_EQUALITY(mesh_info->global_cells(local_cell_id), global_cell_id);
//    }
//  }
}

}

TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 1D_mesh)
{

  int myRank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  auto mesh_info = buildParallelLocalMeshInfo({3},{2},{2},{2.},{0});

  // Make sure there are two blocks (eblock-0, and eblock-1)
  TEST_EQUALITY(mesh_info->element_blocks.size(), 2);
  TEST_EQUALITY(mesh_info->sidesets.size(), 2);

  TEST_ASSERT(mesh_info->has_connectivity);
//  TEST_ASSERT(! mesh_info->cell_sets.is_null());

  if(myRank == 0){

    {
      const auto & block = mesh_info->element_blocks.at("eblock-0");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Element Block eblock-0" << std::endl;

      TEST_ASSERT(block.num_virtual_cells == 0);
//      TEST_ASSERT(not block.cell_sets.is_null());
//      TEST_EQUALITY(block.cell_sets->getNumElements(), block.num_owned_cells+block.num_ghstd_cells+block.num_virtual_cells);
//      for(int i=0; i<block.num_owned_cells; ++i){
//        TEST_EQUALITY(block.cell_sets->getElementSet(i), std::string("eblock-0"));
//      }
//      int num_block_0 = 0, num_block_1 = 0;
//      for(int i=block.num_owned_cells; i<block.num_owned_cells+block.num_ghstd_cells; ++i){
//        if(block.cell_sets->getElementSet(i) == std::string("eblock-0"))
//          ++num_block_0;
//        else if(block.cell_sets->getElementSet(i) == std::string("eblock-1"))
//          ++num_block_1;
//      }
//      TEST_EQUALITY(num_block_0, 1);
//      TEST_EQUALITY(num_block_1, 1);

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 2);
      TEST_EQUALITY(block.num_virtual_cells, 0);
      TEST_EQUALITY((int) global_cells_h(0), 0);
      TEST_EQUALITY((int) global_cells_h(1), 1);
      TEST_EQUALITY((int) global_cells_h(2), 2);
      TEST_EQUALITY((int) global_cells_h(3), 5);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->element_blocks.at("eblock-1");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Element Block eblock-1" << std::endl;

      TEST_ASSERT(block.num_virtual_cells == 0);
//      TEST_ASSERT(not block.cell_sets.is_null());
//      TEST_EQUALITY(block.cell_sets->getNumElements(), block.num_owned_cells+block.num_ghstd_cells+block.num_virtual_cells);
//      for(int i=0; i<block.num_owned_cells; ++i){
//        TEST_EQUALITY(block.cell_sets->getElementSet(i), std::string("eblock-1"));
//      }
//      int num_block_0 = 0, num_block_1 = 0;
//      for(int i=block.num_owned_cells; i<block.num_owned_cells+block.num_ghstd_cells; ++i){
//        if(block.cell_sets->getElementSet(i) == std::string("eblock-0"))
//          ++num_block_0;
//        else if(block.cell_sets->getElementSet(i) == std::string("eblock-1"))
//          ++num_block_1;
//      }
//      TEST_EQUALITY(num_block_0, 1);
//      TEST_EQUALITY(num_block_1, 1);

      // Block should be empty
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 2);
      TEST_EQUALITY(block.num_virtual_cells, 0);
      TEST_EQUALITY((int) global_cells_h(0), 3);
      TEST_EQUALITY((int) global_cells_h(1), 4);
      TEST_EQUALITY((int) global_cells_h(2), 5);
      TEST_EQUALITY((int) global_cells_h(3), 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0").at("left");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Sideset eblock-0 left" << std::endl;

      TEST_ASSERT(block.num_virtual_cells == 0);
//      TEST_ASSERT(not block.cell_sets.is_null());
//      TEST_EQUALITY(block.cell_sets->getNumElements(), block.num_owned_cells+block.num_ghstd_cells+block.num_virtual_cells);
//      TEST_EQUALITY(block.cell_sets->getElementSet(0), std::string("eblock-0"));
//      TEST_EQUALITY(block.cell_sets->getElementSet(1), std::string("eblock-1"));

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 1);
      TEST_EQUALITY(block.num_virtual_cells, 0);
      TEST_EQUALITY((int) global_cells_h(0), 0);
      TEST_EQUALITY((int) global_cells_h(1), 5);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.find("eblock-1") == mesh_info->sidesets.end());

  } else {
    {
      const auto & block = mesh_info->element_blocks.at("eblock-0");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Element Block eblock-0" << std::endl;

      TEST_ASSERT(block.num_virtual_cells == 0);
//      TEST_ASSERT(not block.cell_sets.is_null());
//      TEST_EQUALITY(block.cell_sets->getNumElements(), block.num_owned_cells+block.num_ghstd_cells+block.num_virtual_cells);
//      for(int i=0; i<block.num_owned_cells; ++i){
//        TEST_EQUALITY(block.cell_sets->getElementSet(i), std::string("eblock-0"));
//      }
//      int num_block_0 = 0, num_block_1 = 0;
//      for(int i=block.num_owned_cells; i<block.num_owned_cells+block.num_ghstd_cells; ++i){
//        if(block.cell_sets->getElementSet(i) == std::string("eblock-0"))
//          ++num_block_0;
//        else if(block.cell_sets->getElementSet(i) == std::string("eblock-1"))
//          ++num_block_1;
//      }
//      TEST_EQUALITY(num_block_0, 1);
//      TEST_EQUALITY(num_block_1, 1);

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 2);
      TEST_EQUALITY(block.num_virtual_cells, 0);
      TEST_EQUALITY((int) global_cells_h(0), 2);
      TEST_EQUALITY((int) global_cells_h(1), 3);
      TEST_EQUALITY((int) global_cells_h(2), 1);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->element_blocks.at("eblock-1");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Element Block eblock-1" << std::endl;

      TEST_ASSERT(block.num_virtual_cells == 0);
//      TEST_ASSERT(not block.cell_sets.is_null());
//      TEST_EQUALITY(block.cell_sets->getNumElements(), block.num_owned_cells+block.num_ghstd_cells+block.num_virtual_cells);
//      for(int i=0; i<block.num_owned_cells; ++i){
//        TEST_EQUALITY(block.cell_sets->getElementSet(i), std::string("eblock-1"));
//      }
//      int num_block_0 = 0, num_block_1 = 0;
//      for(int i=block.num_owned_cells; i<block.num_owned_cells+block.num_ghstd_cells; ++i){
//        if(block.cell_sets->getElementSet(i) == std::string("eblock-0"))
//          ++num_block_0;
//        else if(block.cell_sets->getElementSet(i) == std::string("eblock-1"))
//          ++num_block_1;
//      }
//      TEST_EQUALITY(num_block_0, 1);
//      TEST_EQUALITY(num_block_1, 1);

      // Block should be empty
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 2);
      TEST_EQUALITY(block.num_virtual_cells, 0);
      TEST_EQUALITY((int) global_cells_h(0), 5);
      TEST_EQUALITY((int) global_cells_h(1), 0);
      TEST_EQUALITY((int) global_cells_h(2), 4);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.find("eblock-0") == mesh_info->sidesets.end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-1").at("right");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Sideset eblock-1 right" << std::endl;

      TEST_ASSERT(block.num_virtual_cells == 0);
//      TEST_ASSERT(not block.cell_sets.is_null());
//      TEST_EQUALITY(block.cell_sets->getNumElements(), block.num_owned_cells+block.num_ghstd_cells+block.num_virtual_cells);
//      TEST_EQUALITY(block.cell_sets->getElementSet(0), std::string("eblock-1"));
//      TEST_EQUALITY(block.cell_sets->getElementSet(1), std::string("eblock-0"));

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 1);
      TEST_EQUALITY(block.num_virtual_cells, 0);
      TEST_EQUALITY((int) global_cells_h(0), 5);
      TEST_EQUALITY((int) global_cells_h(1), 0);
      TEST_ASSERT(block.has_connectivity);
    }
  }

}

// 1D

TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 1D_connectivity_bc)
{
  runConnectivityTest({3},{2},{2},{1.},{},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 1D_connectivity_bc0)
{
  runConnectivityTest({3},{2},{2},{1.},{0},out,success);
}

// 2D

TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 2D_connectivity_bc)
{
  runConnectivityTest({3,3},{2,1},{2,1},{1.,1.},{},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 2D_connectivity_bc0)
{
  runConnectivityTest({3,3},{2,1},{2,1},{1.,1.},{0},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 2D_connectivity_bc1)
{
  runConnectivityTest({3,3},{2,1},{2,1},{1.,1.},{1},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 2D_connectivity_bc01)
{
  runConnectivityTest({3,3},{2,1},{2,1},{1.,1.},{0,1},out,success);
}

// 3D - not all combinations are tested

TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 3D_connectivity_bc)
{
  runConnectivityTest({3,3,2},{2,1,1},{2,1,1},{1.,1.,1.},{},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 3D_connectivity_bc0)
{
  runConnectivityTest({3,3,2},{2,1,1},{2,1,1},{1.,1.,1.},{0},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 3D_connectivity_bc1)
{
  runConnectivityTest({3,3,2},{2,1,1},{2,1,1},{1.,1.,1.},{1},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 3D_connectivity_bc2)
{
  runConnectivityTest({3,3,2},{2,1,1},{2,1,1},{1.,1.,1.},{2},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 3D_connectivity_bc01)
{
  runConnectivityTest({3,3,2},{2,1,1},{2,1,1},{1.,1.,1.},{0,1},out,success);
}
TEUCHOS_UNIT_TEST(parallelPeriodicLocalMeshUtilities, 3D_connectivity_bc012)
{
  runConnectivityTest({3,3,2},{2,1,1},{2,1,1},{1.,1.,1.},{0,1,2},out,success);
}

}
