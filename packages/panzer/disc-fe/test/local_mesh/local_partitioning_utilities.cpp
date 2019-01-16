// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
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


    std::vector<panzer::LocalMeshPartition > partitions;

    // It should not return any partitions
    generateLocalMeshPartitions(empty_mesh, panzer::blockGhostedDescriptor("block"), partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);
  }

  // Generate a local mesh info object
  Teuchos::RCP<panzer::LocalMeshInfo> mesh_info = generateLocalMeshInfo();

  {
    std::vector<panzer::LocalMeshPartition > partitions;

    // It should throw an error
    TEST_THROW(generateLocalMeshPartitions(*mesh_info, panzer::blockGhostedDescriptor("block0",panzer::WorksetSizeType::NO_ELEMENTS), partitions), std::logic_error);
  }

  // Test the partitions

  // Non-existant block partition
  {
    std::vector<panzer::LocalMeshPartition> partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::blockGhostedDescriptor("block32"), partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }

  // Non-existant sideset partition
  {
    std::vector<panzer::LocalMeshPartition > partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::sidesetGhostedDescriptor("block32","sideset0"), partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }
  {
    std::vector<panzer::LocalMeshPartition > partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::sidesetGhostedDescriptor("block0","sideset32"), partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }
  {
    std::vector<panzer::LocalMeshPartition > partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::sidesetGhostedDescriptor("block0","sideset1"), partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }

  // Existing block partition
  {
    std::vector<panzer::LocalMeshPartition > partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::blockGhostedDescriptor("block1"), partitions);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(partitions.size(), 1);

    const auto & partition = partitions[0];
    TEST_EQUALITY(partition.num_owned_cells,2);
    TEST_EQUALITY(partition.num_ghost_cells,1);
    TEST_EQUALITY(partition.num_virtual_cells,1);
  }
  {
    std::vector<panzer::LocalMeshPartition > partitions;

    // The partitions won't be contiguous, but we will still get partitions
    generateLocalMeshPartitions(*mesh_info, panzer::blockDescriptor("block0"), partitions);

    // We should have only one partition
    TEST_EQUALITY(partitions.size(), 1);

    // NOTE: since we haven't turned on partitioning, there are no ghost or virtual cells
    const auto & partition = partitions[0];
    TEST_EQUALITY(partition.num_owned_cells,2);
    TEST_EQUALITY(partition.num_ghost_cells,0);
    TEST_EQUALITY(partition.num_virtual_cells,0);
  }

  // Existing sideset partition
  {
    std::vector<panzer::LocalMeshPartition > partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::sidesetGhostedDescriptor("block0","sideset0"), partitions);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(partitions.size(), 1);

    const auto & partition = partitions[0];
    TEST_EQUALITY(partition.num_owned_cells,1);
    TEST_EQUALITY(partition.num_ghost_cells,0);
    TEST_EQUALITY(partition.num_virtual_cells,1);
  }
  {
    std::vector<panzer::LocalMeshPartition > partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::sidesetGhostedDescriptor("block1","sideset2"), partitions);

    // Only one partition since we requested a full partition
    TEST_EQUALITY(partitions.size(), 1);

    const auto & partition = partitions[0];
    TEST_EQUALITY(partition.num_owned_cells,1);
    TEST_EQUALITY(partition.num_ghost_cells,1);
    TEST_EQUALITY(partition.num_virtual_cells,0);
  }

  // Existing block double partition
  {
    // We want two partitions, each with a size of 1 cell
    std::vector<panzer::LocalMeshPartition > partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::blockGhostedDescriptor("block1",1), partitions);

    // Two partitions
    TEST_EQUALITY(partitions.size(), 2);

    {
      const auto & partition = partitions[0];
      TEST_EQUALITY(partition.num_owned_cells,1);

      // Either both extra cells are ghosts, or one is virtual and one is ghost
      TEST_EQUALITY(partition.num_ghost_cells+partition.num_virtual_cells,2);
    }

    {
      const auto & partition = partitions[1];
      TEST_EQUALITY(partition.num_owned_cells,1);

      // Either both extra cells are ghosts, or one is virtual and one is ghost
      TEST_EQUALITY(partition.num_ghost_cells+partition.num_virtual_cells,2);
    }
  }

  // Existing sideset double partition, but only one partition is available
  {
    // We want two partitions, each with a size of 1 cell
    std::vector<panzer::LocalMeshPartition > partitions;

    generateLocalMeshPartitions(*mesh_info, panzer::sidesetGhostedDescriptor("block1","sideset1",1), partitions);

    // Only one partition should have been built
    TEST_EQUALITY(partitions.size(), 1);

    {
      const auto & partition = partitions[0];
      TEST_EQUALITY(partition.num_owned_cells,1);
      TEST_EQUALITY(partition.num_ghost_cells,0);
      TEST_EQUALITY(partition.num_virtual_cells,1);
    }
  }

}

} // end namespace panzer
