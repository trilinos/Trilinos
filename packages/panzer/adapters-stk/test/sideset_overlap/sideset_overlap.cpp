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

#include "Panzer_STK_CheckSidesetOverlap.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"

namespace panzer_test {

  TEUCHOS_UNIT_TEST(periodic_bcs, sidesetOverlap_2D)
  {
    // setup mesh
    Teuchos::RCP<panzer_stk::STK_Interface> mesh;
    {
      Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
      pl->set("X Blocks",2);
      pl->set("Y Blocks",1);
      pl->set("X Elements",2);
      pl->set("Y Elements",2);
      pl->set("X Procs",2);
      pl->set("Y Procs",2);
      panzer_stk::SquareQuadMeshFactory mesh_factory;
      mesh_factory.setParameterList(pl);
      mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // Test is hard coded to four processes to make sure we cover
    // parallel split of sidesets
    TEST_ASSERT(mesh->getComm()->getSize() == 4);

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("left","right",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","bottom",*mesh) );

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("top","bottom",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","bottom",*mesh) );
  }

  TEUCHOS_UNIT_TEST(periodic_bcs, sidesetOverlap_3D)
  {
    // setup mesh
    Teuchos::RCP<panzer_stk::STK_Interface> mesh;
    {
      Teuchos::RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList);
       pl->set("X Blocks",2);
       pl->set("Y Blocks",2);
       pl->set("Z Blocks",2);
       pl->set("X Elements",2);
       pl->set("Y Elements",2);
       pl->set("Z Elements",2);
       pl->set("X Procs",1);
       pl->set("Y Procs",2);
       pl->set("Z Procs",2);
       panzer_stk::CubeHexMeshFactory mesh_factory;
       mesh_factory.setParameterList(pl);
       mesh = mesh_factory.buildMesh(MPI_COMM_WORLD);
    }

    // Test is hard coded to four processes to make sure we cover
    // parallel split of sidesets
    TEST_ASSERT(mesh->getComm()->getSize() == 4);

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("left","right",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","bottom",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","front",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("left","back",*mesh) );

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("top","bottom",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","bottom",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","front",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("right","back",*mesh) );

    TEST_ASSERT( !panzer_stk::checkSidesetOverlap("back","front",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("back","left",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("back","right",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("back","top",*mesh) );
    TEST_ASSERT( panzer_stk::checkSidesetOverlap("back","bottom",*mesh) );
  }
}
