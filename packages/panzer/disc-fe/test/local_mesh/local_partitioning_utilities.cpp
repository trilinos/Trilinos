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

using LO=int;
using GO=panzer::Ordinal64;

TEUCHOS_UNIT_TEST(localMeshPartitioningUtilities, basic)
{

  // Test empty mesh
  {
    panzer::LocalMeshInfo<LO,GO> empty_mesh;

    panzer::WorksetDescriptor description("block",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

    // It should not return any partitions
    generateLocalMeshPartitions(empty_mesh, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);
  }

  // Generate a local mesh info object
  Teuchos::RCP<panzer::LocalMeshInfo<LO,GO>> mesh_info = generateLocalMeshInfo<LO,GO>();

  // Test bad descriptors
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::CLASSIC_MODE,true);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

    // It should throw an error
    TEST_THROW(generateLocalMeshPartitions(*mesh_info, description, partitions), std::logic_error);
  }
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::NO_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

    // It should throw an error
    TEST_THROW(generateLocalMeshPartitions(*mesh_info, description, partitions), std::logic_error);
  }
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::ALL_ELEMENTS,false);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

    // It should throw an error
    TEST_THROW(generateLocalMeshPartitions(*mesh_info, description, partitions), std::logic_error);
  }

  // Test the partitions

  // Non-existant block partition
  {
    panzer::WorksetDescriptor description("block32",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }

  // Non-existant sideset partition
  {
    panzer::WorksetDescriptor description("block32","sideset0",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }
  {
    panzer::WorksetDescriptor description("block0","sideset32",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }
  {
    panzer::WorksetDescriptor description("block0","sideset1",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

    generateLocalMeshPartitions(*mesh_info, description, partitions);

    // Nothing should have been found
    TEST_EQUALITY(partitions.size(), 0);

  }

  // Existing block partition
  {
    panzer::WorksetDescriptor description("block0",panzer::WorksetSizeType::ALL_ELEMENTS,true);
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

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
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

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
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

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
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

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
    std::vector<panzer::LocalMeshPartition<LO,GO> > partitions;

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
