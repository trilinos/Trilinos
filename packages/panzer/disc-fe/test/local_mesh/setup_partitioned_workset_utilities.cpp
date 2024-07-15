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

#include "PanzerCore_config.hpp"

#include "Panzer_Workset.hpp"

#include "Panzer_WorksetNeeds.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_IntegrationDescriptor.hpp"

#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_SetupPartitionedWorksetUtilities.hpp"

#include "Panzer_UnitTest_LocalMeshUtilities.hpp"

#include <vector>
#include <string>
#include <iostream>

namespace panzer {

TEUCHOS_UNIT_TEST(setupPartitionedWorksetUtilities, basic)
{

  // Test empty mesh
  {
    panzer::LocalMeshInfo empty_mesh;

    panzer::WorksetDescriptor description("block",panzer::WorksetSizeType::ALL_ELEMENTS,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(empty_mesh, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    TEST_EQUALITY(worksets->size(),0);

  }

  // Generate a local mesh info object
  Teuchos::RCP<panzer::LocalMeshInfo> mesh_info = generateLocalMeshInfo();

  // Test bad descriptors
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::CLASSIC_MODE,true);

    // Should throw an error
    TEST_THROW(buildPartitionedWorksets(*mesh_info, description),std::logic_error);

  }
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::NO_ELEMENTS,true);

    // It shouldn't return any worksets
    TEST_THROW(buildPartitionedWorksets(*mesh_info, description),std::logic_error);
  }
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::ALL_ELEMENTS,false);

    // This should trow an error since the description doesn't allow partitioning
    TEST_THROW(buildPartitionedWorksets(*mesh_info, description),std::logic_error);
  }

  // Test the worksets

  // Non-existant block
  {
    panzer::WorksetDescriptor description("block32",panzer::WorksetSizeType::ALL_ELEMENTS,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Nothing should have been found
    TEST_EQUALITY(worksets->size(), 0);

  }

  // Non-existant sideset
  {
    panzer::WorksetDescriptor description("block32","sideset0",panzer::WorksetSizeType::ALL_ELEMENTS,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Nothing should have been found
    TEST_EQUALITY(worksets->size(), 0);

  }
  {
    panzer::WorksetDescriptor description("block0","sideset32",panzer::WorksetSizeType::ALL_ELEMENTS,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Nothing should have been found
    TEST_EQUALITY(worksets->size(), 0);

  }
  {
    panzer::WorksetDescriptor description("block0","sideset1",panzer::WorksetSizeType::ALL_ELEMENTS,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Nothing should have been found
    TEST_EQUALITY(worksets->size(), 0);
  }

  // Existing block partition
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::ALL_ELEMENTS,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(worksets->size(), 1);

    const auto & workset = (*worksets)[0];
    TEST_EQUALITY(workset.numOwnedCells(),2);
    TEST_EQUALITY(workset.numGhostCells(),1);
    TEST_EQUALITY(workset.numVirtualCells(),1);
  }

  // Existing sideset partition
  {
    panzer::WorksetDescriptor description("block0","sideset0",panzer::WorksetSizeType::ALL_ELEMENTS,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(worksets->size(), 1);

    const auto & workset = (*worksets)[0];
    TEST_EQUALITY(workset.numOwnedCells(),1);
    TEST_EQUALITY(workset.numGhostCells(),0);
    TEST_EQUALITY(workset.numVirtualCells(),1);
  }
  {
    panzer::WorksetDescriptor description("block1","sideset2",panzer::WorksetSizeType::ALL_ELEMENTS,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(worksets->size(), 1);

    const auto & workset = (*worksets)[0];
    TEST_EQUALITY(workset.numOwnedCells(),1);
    TEST_EQUALITY(workset.numGhostCells(),1);
    TEST_EQUALITY(workset.numVirtualCells(),0);
  }

  // Existing block double partition
  {
    // We want two partitions, each with a size of 1 cell
    panzer::WorksetDescriptor description("block1",1,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Two partitions
    TEST_EQUALITY(worksets->size(), 2);

    {
      const auto & workset = (*worksets)[0];
      TEST_EQUALITY(workset.numOwnedCells(),1);

      // Either both extra cells are ghosts, or one is virtual and one is ghost
      TEST_EQUALITY(workset.numGhostCells()+workset.numVirtualCells(),2);
    }

    {
      const auto & workset = (*worksets)[1];
      TEST_EQUALITY(workset.numOwnedCells(),1);

      // Either both extra cells are ghosts, or one is virtual and one is ghost
      TEST_EQUALITY(workset.numGhostCells()+workset.numVirtualCells(),2);
    }
  }

  // Existing sideset double partition, but only one partition is available
  {
    // We want two partitions, each with a size of 1 cell
    panzer::WorksetDescriptor description("block1","sideset1",1,true);

    // It shouldn't return any worksets
    auto worksets = buildPartitionedWorksets(*mesh_info, description);

    TEST_INEQUALITY(worksets,Teuchos::null);

    // Only one partition should have been built
    TEST_EQUALITY(worksets->size(), 1);

    {
      const auto & workset = (*worksets)[0];
      TEST_EQUALITY(workset.numOwnedCells(),1);
      TEST_EQUALITY(workset.numGhostCells(),0);
      TEST_EQUALITY(workset.numVirtualCells(),1);
    }
  }


}

} // end namespace panzer
