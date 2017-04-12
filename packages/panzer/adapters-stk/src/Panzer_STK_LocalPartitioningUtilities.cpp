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

#include "Panzer_STK_LocalPartitioningUtilities.hpp"

#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"

#include "Panzer_STK_Interface.hpp"
#include "Panzer_Workset_Builder.hpp"

#include "Panzer_WorksetDescriptor.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_Assert.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"

#include "Panzer_STK_LocalMeshUtilities.hpp"

#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <vector>

namespace panzer_stk
{

namespace partitioning_utils
{

template<typename LO, typename GO>
void
splitMeshInfo(const panzer::LocalMeshInfoBase<LO,GO> & mesh_info,
              const int splitting_size,
              std::vector<panzer::LocalMeshPartition<LO,GO> > & partitions)
{

  // This is not a partitioning scheme.
  // This just chunks the mesh_info into equally sized chunks and ignores connectivity

  const LO num_owned_cells = mesh_info.num_owned_cells;

  std::vector<LO> partition_cells;

  if(splitting_size < 0){
    // Just one chunk
    partition_cells.reserve(mesh_info.num_owned_cells);
    for(LO i=0;i<mesh_info.num_owned_cells;++i){
      partition_cells.push_back(i);
    }

    partitions.push_back(panzer::LocalMeshPartition<LO,GO>());
    panzer::LocalMeshPartition<LO,GO> & partition_info = partitions.back();

    tools::setupSubLocalMeshInfo(mesh_info,partition_cells,partition_info);
  } else {

    partition_cells.reserve(splitting_size);

    const LO num_partitions = mesh_info.num_owned_cells / splitting_size + ((mesh_info.num_owned_cells % splitting_size == 0) ? 0 : 1);
    LO partition_start_index = 0;

    for(LO partition=0;partition<num_partitions;++partition){
      const LO partition_end_index = std::min(partition_start_index + splitting_size, num_owned_cells);
      partition_cells.resize(partition_end_index - partition_start_index,-1);

      for(int i=partition_start_index;i<partition_end_index;++i){
        partition_cells[i] = i;
      }

      partitions.push_back(panzer::LocalMeshPartition<LO,GO>());
      panzer::LocalMeshPartition<LO,GO> & partition_info = partitions.back();

      tools::setupSubLocalMeshInfo(mesh_info,partition_cells,partition_info);

      partition_start_index = partition_end_index;
    }

  }

}

template<typename LO, typename GO>
void
partitionMeshInfo(const panzer::LocalMeshInfoBase<LO,GO> & mesh_info,
                 const size_t requested_partition_size,
                 std::vector<panzer::LocalMeshPartition<LO,GO> > & partitions)
{
  // Not yet sure how to do this
  TEUCHOS_ASSERT(false);
}

}

template <typename LO, typename GO>
void
generateLocalMeshPartitions(const panzer_stk::STK_Interface & mesh,
                            const panzer::WorksetDescriptor & description,
                            std::vector<panzer::LocalMeshPartition<LO,GO> > & partitions)
{

  // We have to make sure that the partitioning is possible
  TEUCHOS_ASSERT(description.getWorksetSize() != panzer::WorksetDescriptor::EMPTY);
  TEUCHOS_ASSERT(description.getWorksetSize() != 0);

  const std::string & element_block_name = description.getElementBlock();

  // Generate the mesh info describing the local mesh on this process
  // TODO: This needs to be a singleton, or belong to a class that is allocated once
  // Generating the local mesh info is an expensive process
  // It is largely topology stuff, and we can have a separate update call if the geometry of the mesh changes
  panzer::LocalMeshInfo<LO,GO> mesh_info;

  // Very expensive
  panzer_stk::generateLocalMeshInfo<LO,GO>(mesh, mesh_info);

  // We have two processes for in case this is a sideset or element block
  if(description.useSideset()){

    const std::string & sideset_name = description.getSideset();

    TEUCHOS_ASSERT(mesh_info.sidesets.find(element_block_name) != mesh_info.sidesets.end());

    const auto & sideset_map = mesh_info.sidesets.at(element_block_name);

    TEUCHOS_ASSERT(sideset_map.find(sideset_name) != sideset_map.end());

    const panzer::LocalMeshSidesetInfo<LO,GO> & sideset_info = sideset_map.at(sideset_name);

    // Partitioning is not important for sidesets
    panzer_stk::partitioning_utils::splitMeshInfo<LO,GO>(sideset_info, description.getWorksetSize(), partitions);

    for(auto & partition : partitions){
      partition.sideset_name = sideset_name;
      partition.element_block_name = element_block_name;
      partition.cell_topology = mesh.getCellTopology(element_block_name);
    }

  } else {

    // Make sure the element block of interest exists
    TEUCHOS_ASSERT(mesh_info.element_blocks.find(element_block_name) != mesh_info.element_blocks.end());

    // Grab the element block we're interested in
    const panzer::LocalMeshBlockInfo<LO,GO> & block_info = mesh_info.element_blocks[element_block_name];


    if(description.getWorksetSize() == panzer::WorksetDescriptor::FULL){
      // We only have one partition describing the entire local mesh
      panzer_stk::partitioning_utils::splitMeshInfo(block_info, -1, partitions);
    } else {
      // We need to partition local mesh
      panzer_stk::partitioning_utils::partitionMeshInfo<LO,GO>(block_info, description.getWorksetSize(), partitions);
    }

    for(auto & partition : partitions){
      partition.element_block_name = element_block_name;
      partition.cell_topology = mesh.getCellTopology(element_block_name);
    }
  }

}

}

template
void
panzer_stk::generateLocalMeshPartitions<int,int>(const panzer_stk::STK_Interface & mesh,
                                                 const panzer::WorksetDescriptor & description,
                                                 std::vector<panzer::LocalMeshPartition<int,int> > & partitions);

#ifndef PANZER_ORDINAL64_IS_INT
template
void
panzer_stk::generateLocalMeshPartitions<int,panzer::Ordinal64>(const panzer_stk::STK_Interface & mesh,
                                                               const panzer::WorksetDescriptor & description,
                                                               std::vector<panzer::LocalMeshPartition<int,panzer::Ordinal64> > & partitions);
#endif
