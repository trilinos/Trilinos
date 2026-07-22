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

#include "Panzer_WorksetDescriptor.hpp"

#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_LocalPartitioningUtilities.hpp"

#include "Panzer_UnitTest_LocalMeshUtilities.hpp"

#include <vector>
#include <string>
#include <iostream>

namespace panzer {

TEUCHOS_UNIT_TEST(localMeshPartitioningUtilities, basic)
{

  // Test empty mesh
  {
    panzer::LocalMeshInfo empty_mesh;

    panzer::WorksetDescriptor description("block",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    // It should not return any partitions
    generateLocalMeshPartitions(empty_mesh, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);
  }

  // Generate a local mesh info object
  Teuchos::RCP<panzer::LocalMeshInfo> mesh_info = generateLocalMeshInfo();

  // Test bad descriptors
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::CLASSIC_MODE,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    // It should throw an error
    TEST_THROW(generateLocalMeshPartitions(*mesh_info, description, partitions), std::logic_error);
  }
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::NO_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    // It should throw an error
    TEST_THROW(generateLocalMeshPartitions(*mesh_info, description, partitions), std::logic_error);
  }
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::ALL_ELEMENTS,false);
    std::vector<panzer::LocalMeshPartition> partitions;

    // It should throw an error
    TEST_THROW(generateLocalMeshPartitions(*mesh_info, description, partitions), std::logic_error);
  }

  // Test the partitions

  // Non-existant block partition
  {
    panzer::WorksetDescriptor description("block32",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }

  // Non-existant sideset partition
  {
    panzer::WorksetDescriptor description("block32","sideset0",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }
  {
    panzer::WorksetDescriptor description("block0","sideset32",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }
  {
    panzer::WorksetDescriptor description("block0","sideset1",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }

  // Existing block partition
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(partitions.size(), 1);

    const auto & partition = partitions[0];
    TEST_EQUALITY(partition.num_owned_cells,2);
    TEST_EQUALITY(partition.num_ghstd_cells,1);
    TEST_EQUALITY(partition.num_virtual_cells,1);
  }

  // Existing sideset partition
  {
    panzer::WorksetDescriptor description("block0","sideset0",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(partitions.size(), 1);

    const auto & partition = partitions[0];
    TEST_EQUALITY(partition.num_owned_cells,1);
    TEST_EQUALITY(partition.num_ghstd_cells,0);
    TEST_EQUALITY(partition.num_virtual_cells,1);
  }
  {
    panzer::WorksetDescriptor description("block1","sideset2",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(partitions.size(), 1);

    const auto & partition = partitions[0];
    TEST_EQUALITY(partition.num_owned_cells,1);
    TEST_EQUALITY(partition.num_ghstd_cells,1);
    TEST_EQUALITY(partition.num_virtual_cells,0);
  }

  // Existing block double partition
  {
    // We want two partitions, each with a size of 1 cell
    panzer::WorksetDescriptor description("block1",1,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Two partitions
    TEST_EQUALITY(partitions.size(), 2);

    {
      const auto & partition = partitions[0];
      TEST_EQUALITY(partition.num_owned_cells,1);

      // Either both extra cells are ghosts, or one is virtual and one is ghost
      TEST_EQUALITY(partition.num_ghstd_cells+partition.num_virtual_cells,2);
    }

    {
      const auto & partition = partitions[1];
      TEST_EQUALITY(partition.num_owned_cells,1);

      // Either both extra cells are ghosts, or one is virtual and one is ghost
      TEST_EQUALITY(partition.num_ghstd_cells+partition.num_virtual_cells,2);
    }
  }

  // Existing sideset double partition, but only one partition is available
  {
    // We want two partitions, each with a size of 1 cell
    panzer::WorksetDescriptor description("block1","sideset1",1,true);
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Only one partition should have been built
    TEST_EQUALITY(partitions.size(), 1);

    {
      const auto & partition = partitions[0];
      TEST_EQUALITY(partition.num_owned_cells,1);
      TEST_EQUALITY(partition.num_ghstd_cells,0);
      TEST_EQUALITY(partition.num_virtual_cells,1);
    }
  }

}

} // end namespace panzer
