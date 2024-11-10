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
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Tpetra_Core.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"

#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_STK_ExodusReaderFactory.hpp"

#include "PanzerSTK_UnitTest_BuildMesh.hpp"

#include <string>

namespace panzer_stk {
namespace {

Teuchos::RCP<const panzer::LocalMeshInfo>
buildParallelLocalMeshInfo(const std::vector<int> & N,   // Cells per dimension
                           const std::vector<int> & B,   // Blocks per dimension
                           const std::vector<int> & P,   // Processors per dimension
                           const std::vector<double> &L) // Domain length per block
{
  std::vector<int> p;
//  for(unsigned int i=0; i<N.size(); ++i)
//    p.push_back(i);
  Teuchos::RCP<panzer_stk::STK_Interface> mesh = buildParallelMesh(N,B,P,L,p);
  return generateLocalMeshInfo(*mesh);
}

}

TEUCHOS_UNIT_TEST(parallelLocalMeshUtilities, no_mesh)
{
  // Make sure if fails when you pass in a null mesh
  panzer_stk::STK_Interface mesh;
  TEST_THROW((generateLocalMeshInfo(mesh)),std::logic_error);
}

TEUCHOS_UNIT_TEST(parallelLocalMeshUtilities, 1D_mesh)
{

  int myRank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  auto mesh_info = buildParallelLocalMeshInfo({3},{2},{2},{2.});

  // Make sure there are two blocks (eblock-0, and eblock-1)
  TEST_EQUALITY(mesh_info->element_blocks.size(), 2);
  TEST_EQUALITY(mesh_info->sidesets.size(), 2);

  TEST_ASSERT(mesh_info->has_connectivity);

  if(myRank == 0){

    {
      const auto & block = mesh_info->element_blocks.at("eblock-0");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Element Block eblock-0" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 1);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_EQUALITY((int) global_cells_h(0), 0);
      TEST_EQUALITY((int) global_cells_h(1), 1);
      TEST_EQUALITY((int) global_cells_h(2), 2);
      TEST_EQUALITY((int) global_cells_h(3), 6);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->element_blocks.at("eblock-1");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Element Block eblock-1" << std::endl;

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

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_EQUALITY((int) global_cells_h(0), 0);
      TEST_EQUALITY((int) global_cells_h(1), 6);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.find("eblock-1") == mesh_info->sidesets.end());

  } else {
    {
      const auto & block = mesh_info->element_blocks.at("eblock-0");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Element Block eblock-0" << std::endl;

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

      // Block should be empty
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 1);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_EQUALITY((int) global_cells_h(0), 5);
      TEST_EQUALITY((int) global_cells_h(1), 4);
      TEST_EQUALITY((int) global_cells_h(2), 7);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.find("eblock-0") == mesh_info->sidesets.end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-1").at("right");
      auto global_cells_h = Kokkos::create_mirror_view(block.global_cells);
      Kokkos::deep_copy(global_cells_h, block.global_cells);

      out << "Sideset eblock-1 right" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_EQUALITY((int) global_cells_h(0), 5);
      TEST_EQUALITY((int) global_cells_h(1), 7);
      TEST_ASSERT(block.has_connectivity);
    }
  }

}

TEUCHOS_UNIT_TEST(parallelLocalMeshUtilities, 2D_mesh)
{

  int myRank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  auto mesh_info = buildParallelLocalMeshInfo({2,2},{2,1},{2,1},{2.,2.});

  // Make sure there are two blocks (eblock-0, and eblock-1)
  TEST_EQUALITY(mesh_info->element_blocks.size(), 2);
  TEST_EQUALITY(mesh_info->sidesets.size(), 2);
  TEST_ASSERT(mesh_info->has_connectivity);

  if(myRank == 0){

    {
      const auto & block = mesh_info->element_blocks.at("eblock-0_0");

      out << "Element Block eblock-0_0" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 2);
      TEST_EQUALITY(block.num_virtual_cells, 4);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->element_blocks.at("eblock-1_0");

      out << "Element Block eblock-1_0" << std::endl;

      // Block should be empty
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 4);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0").at("left");

      out << "Sideset eblock-0_0 left" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0").at("top");

      out << "Sideset eblock-0_0 top" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0").at("bottom");

      out << "Sideset eblock-0_0 bottom" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-0_0").find("right") == mesh_info->sidesets.at("eblock-0_0").end());

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-1_0").find("left") == mesh_info->sidesets.at("eblock-1_0").end());

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-1_0").find("right") == mesh_info->sidesets.at("eblock-1_0").end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0").at("top");

      out << "Sideset eblock-1_0 top" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0").at("bottom");

      out << "Sideset eblock-1_0 bottom" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_ASSERT(block.has_connectivity);
    }

  } else {

    {
      const auto & block = mesh_info->element_blocks.at("eblock-0_0");

      out << "Element Block eblock-0_0" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 4);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->element_blocks.at("eblock-1_0");

      out << "Element Block eblock-1_0" << std::endl;

      // Block should be empty
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 2);
      TEST_EQUALITY(block.num_virtual_cells, 4);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-0_0").find("left") == mesh_info->sidesets.at("eblock-0_0").end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0").at("top");

      out << "Sideset eblock-0_0 top" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0").at("bottom");

      out << "Sideset eblock-0_0 bottom" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-0_0").find("right") == mesh_info->sidesets.at("eblock-0_0").end());

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-1_0").find("left") == mesh_info->sidesets.at("eblock-1_0").end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0").at("right");

      out << "Sideset eblock-1_0 right" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0").at("top");

      out << "Sideset eblock-1_0 top" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0").at("bottom");

      out << "Sideset eblock-1_0 bottom" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 1);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 1);
      TEST_ASSERT(block.has_connectivity);
    }
  }

}


TEUCHOS_UNIT_TEST(parallelLocalMeshUtilities, 3D_mesh)
{

  int myRank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  auto mesh_info = buildParallelLocalMeshInfo({2,2,2},{2,1,1},{2,1,1},{2.,2.,2.});

  // Make sure there are two blocks (eblock-0, and eblock-1)
  TEST_EQUALITY(mesh_info->element_blocks.size(), 2);
  TEST_EQUALITY(mesh_info->sidesets.size(), 2);
  TEST_ASSERT(mesh_info->has_connectivity);

  if(myRank == 0){

    {
      const auto & block = mesh_info->element_blocks.at("eblock-0_0_0");

      out << "Element Block eblock-0_0_0" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 4);
      TEST_EQUALITY(block.num_ghstd_cells, 4);
      TEST_EQUALITY(block.num_virtual_cells, 12);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->element_blocks.at("eblock-1_0_0");

      out << "Element Block eblock-1_0_0" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 4);
      TEST_EQUALITY(block.num_ghstd_cells, 8);
      TEST_EQUALITY(block.num_virtual_cells, 8);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("left");

      out << "Sideset eblock-0_0_0 left" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 4);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 4);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("top");

      out << "Sideset eblock-0_0_0 top" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("bottom");

      out << "Sideset eblock-0_0_0 bottom" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("front");

      out << "Sideset eblock-0_0_0 front" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("back");

      out << "Sideset eblock-0_0_0 back" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-0_0_0").find("right") == mesh_info->sidesets.at("eblock-0_0_0").end());

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-1_0_0").find("left") == mesh_info->sidesets.at("eblock-1_0_0").end());

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-1_0_0").find("right") == mesh_info->sidesets.at("eblock-1_0_0").end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("top");

      out << "Sideset eblock-1_0_0 top" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("bottom");

      out << "Sideset eblock-1_0_0 bottom" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("front");

      out << "Sideset eblock-1_0_0 front" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("back");

      out << "Sideset eblock-1_0_0 back" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

  } else {

    {
      const auto & block = mesh_info->element_blocks.at("eblock-0_0_0");

      out << "Element Block eblock-0_0_0" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 4);
      TEST_EQUALITY(block.num_ghstd_cells, 8);
      TEST_EQUALITY(block.num_virtual_cells, 8);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->element_blocks.at("eblock-1_0_0");

      out << "Element Block eblock-1_0_0" << std::endl;

      // Block should be empty
      TEST_EQUALITY(block.num_owned_cells, 4);
      TEST_EQUALITY(block.num_ghstd_cells, 4);
      TEST_EQUALITY(block.num_virtual_cells, 12);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-0_0_0").find("left") == mesh_info->sidesets.at("eblock-0_0_0").end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("top");

      out << "Sideset eblock-0_0_0 top" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("bottom");

      out << "Sideset eblock-0_0_0 bottom" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-0_0_0").find("right") == mesh_info->sidesets.at("eblock-0_0_0").end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("front");

      out << "Sideset eblock-0_0_0 front" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-0_0_0").at("back");

      out << "Sideset eblock-0_0_0 back" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

//    TEUCHOS_ASSERT(mesh_info->sidesets.at("eblock-1_0_0").find("left") == mesh_info->sidesets.at("eblock-1_0_0").end());

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("right");

      out << "Sideset eblock-1_0_0 right" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 4);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 4);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("top");

      out << "Sideset eblock-1_0_0 top" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("bottom");

      out << "Sideset eblock-1_0_0 bottom" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("front");

      out << "Sideset eblock-1_0_0 front" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }

    {
      const auto & block = mesh_info->sidesets.at("eblock-1_0_0").at("back");

      out << "Sideset eblock-1_0_0 back" << std::endl;

      // Block should have some basic stuff working
      TEST_EQUALITY(block.num_owned_cells, 2);
      TEST_EQUALITY(block.num_ghstd_cells, 0);
      TEST_EQUALITY(block.num_virtual_cells, 2);
      TEST_ASSERT(block.has_connectivity);
    }
  }

}

}
