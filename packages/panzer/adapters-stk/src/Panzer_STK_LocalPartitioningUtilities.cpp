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
convertMeshInfoToSingleChunk(const panzer::LocalMeshInfo<LO,GO> & mesh_info,
                             panzer::LocalMeshChunk<LO,GO> & chunk)
{

  const GO num_owned_cells = mesh_info.owned_cells.dimension(0);
  const GO num_ghost_cells = mesh_info.ghstd_cells.dimension(0);
  const GO num_virtual_cells = mesh_info.virtual_cells.dimension(0);
  const GO num_total_cells = num_owned_cells+num_ghost_cells+num_virtual_cells;
  const GO num_vertices_per_cell = mesh_info.owned_vertices.dimension(1);
  const GO num_dims_per_vertex = mesh_info.owned_vertices.dimension(2);

  chunk.cell_topology = mesh_info.cell_topology;

  chunk.element_block_name = mesh_info.element_block_name;

  chunk.local_mesh_cell_indexes.resize(num_total_cells,-1);
  for(GO i=0;i<num_total_cells;++i){chunk.local_mesh_cell_indexes[i]=LO(i);}

  chunk.owned_cell_global_indexes.resize(num_owned_cells,-1);
  for(GO i=0;i<num_owned_cells;++i){chunk.owned_cell_global_indexes[i]=mesh_info.owned_cells(i);}

  chunk.ghost_cell_global_indexes.resize(num_ghost_cells,-1);
  for(GO i=0;i<num_ghost_cells;++i){chunk.ghost_cell_global_indexes[i]=mesh_info.ghstd_cells(i);}

  chunk.virtual_cell_local_indexes.resize(num_virtual_cells,-1);
  for(GO i=0;i<num_virtual_cells;++i){chunk.virtual_cell_local_indexes[i]=mesh_info.virtual_cells(i);}

  chunk.cell_vertices = Kokkos::View<double***,PHX::Device>("cell_vertices",num_total_cells,num_vertices_per_cell,num_dims_per_vertex);
  Kokkos::deep_copy(chunk.cell_vertices, 0.);
  for(GO i=0;i<num_owned_cells;++i){
    for(GO j=0;j<num_vertices_per_cell;++j){
      for(GO k=0;k<num_dims_per_vertex;++k){
        chunk.cell_vertices(i,j,k) = mesh_info.owned_vertices(i,j,k);
      }
    }
  }
  for(GO i=0;i<num_ghost_cells;++i){
    for(GO j=0;j<num_vertices_per_cell;++j){
      for(GO k=0;k<num_dims_per_vertex;++k){
        chunk.cell_vertices(i+num_owned_cells,j,k) = mesh_info.ghstd_vertices(i,j,k);
      }
    }
  }

  chunk.face_to_cells = mesh_info.face_to_cells;
  chunk.cell_to_faces = mesh_info.cell_to_face;
  chunk.face_to_local_faces = mesh_info.face_to_lidx;
  chunk.num_cells = num_total_cells;
  chunk.num_faces = chunk.face_to_cells.dimension_0();
}

template<typename LO, typename GO>
void
paritionMeshInfoIntoChunks(const panzer::LocalMeshInfo<LO,GO> & mesh_info, const size_t & requested_partition_size, std::vector<panzer::LocalMeshChunk<LO,GO> > & chunks)
{
  // Not yet sure how to do this
  TEUCHOS_ASSERT(false);
}

}

template <typename LO, typename GO>
std::vector<panzer::LocalMeshChunk<LO,GO> >
generateLocalMeshChunks(const panzer_stk::STK_Interface & mesh,
                        const panzer::WorksetDescriptor & description)
{

  // We have to make sure that the partitioning is possible
  TEUCHOS_ASSERT(description.getWorksetSize() != panzer::WorksetDescriptor::EMPTY);
  TEUCHOS_ASSERT(description.getWorksetSize() != 0);

  // Generate the mesh info describing the local mesh on this process
  panzer::LocalMeshInfo<LO,GO> mesh_info;
  if(description.useSideset()){
    // We don't have this interface yet
    TEUCHOS_ASSERT(! description.connectsElementBlocks());

    // This is a sideset
    mesh_info = panzer_stk::generateLocalSidesetInfo<LO,GO>(mesh, description.getElementBlock(), description.getSideset());
  } else {
    // This is a cellset
    mesh_info = panzer_stk::generateLocalMeshInfo<LO,GO>(mesh, description.getElementBlock());
  }

  // Generate a set of chunks for this process
  std::vector<panzer::LocalMeshChunk<LO,GO> > chunks;
  if(description.getWorksetSize() == panzer::WorksetDescriptor::FULL){
    // We only have one chunk describing the entire local mesh
    chunks.push_back(panzer::LocalMeshChunk<LO,GO>());
    partitioning_utils::convertMeshInfoToSingleChunk(mesh_info, chunks.back());
  } else {
    // We need to partition local mesh into smaller chunks
    partitioning_utils::paritionMeshInfoIntoChunks(mesh_info, description.getWorksetSize(), chunks);
  }

  return chunks;
}

}

template
std::vector<panzer::LocalMeshChunk<int,int> >
panzer_stk::generateLocalMeshChunks<int,int>(const panzer_stk::STK_Interface & mesh,
                        const panzer::WorksetDescriptor & description);

#ifndef PANZER_ORDINAL64_IS_INT
template
std::vector<panzer::LocalMeshChunk<int,panzer::Ordinal64> >
panzer_stk::generateLocalMeshChunks<int,panzer::Ordinal64>(const panzer_stk::STK_Interface & mesh,
                        const panzer::WorksetDescriptor & description);
#endif
