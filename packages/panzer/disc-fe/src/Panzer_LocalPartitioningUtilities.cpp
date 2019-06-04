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

#include "Panzer_LocalPartitioningUtilities.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_Assert.hpp"

#include "Panzer_Workset_Builder.hpp"
#include "Panzer_WorksetDescriptor.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include <unordered_set>
#include <unordered_map>

namespace panzer
{

namespace partitioning_utilities
{

template<typename LO, typename GO>
void
setupSubLocalMeshInfo(const panzer::LocalMeshInfoBase<LO,GO> & parent_info,
                      const std::vector<LO> & owned_parent_cells,
                      panzer::LocalMeshInfoBase<LO,GO> & sub_info)
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::partitioning_utilities::setupSubLocalMeshInfo",setupSLMI);
  // The goal of this function is to fill a LocalMeshInfoBase (sub_info) with
  // a subset of cells from a given parent LocalMeshInfoBase (parent_info)

  // Note: owned_parent_cells are the owned cells for sub_info in the parent_info's indexing scheme
  // We need to generate sub_info's ghosts and figure out the virtual cells

  // Note: We will only handle a single ghost layer

  // Note: We assume owned_parent_cells are owned cells of the parent
  // i.e. owned_parent_indexes cannot refer to ghost or virtual cells in parent_info

  // Note: This function works with inter-face connectivity. NOT node connectivity.

  const int num_owned_cells = owned_parent_cells.size();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(num_owned_cells == 0, "panzer::partitioning_utilities::setupSubLocalMeshInfo : Input parent subcells must exist (owned_parent_cells)");

  const int num_parent_owned_cells = parent_info.num_owned_cells;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(num_parent_owned_cells <= 0, "panzer::partitioning_utilities::setupSubLocalMeshInfo : Input parent info must contain owned cells");

  const int num_parent_ghstd_cells = parent_info.num_ghstd_cells;

  const int num_parent_total_cells = parent_info.num_owned_cells + parent_info.num_ghstd_cells + parent_info.num_virtual_cells;

  // Just as a precaution, make sure the parent_info is setup properly
  TEUCHOS_ASSERT(static_cast<int>(parent_info.cell_to_faces.extent(0)) == num_parent_total_cells);
  const int num_faces_per_cell = parent_info.cell_to_faces.extent(1);

  // The first thing to do is construct a vector containing the parent cell indexes of all
  // owned, ghstd, and virtual cells
  std::vector<LO> ghstd_parent_cells;
  std::vector<LO> virtual_parent_cells;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Construct parent cell vector",ParentCell);
    // We grab all of the owned cells and put their global indexes into sub_info
    // We also put all of the owned cell indexes in the parent's indexing scheme into a set to use for lookups
    std::unordered_set<LO> owned_parent_cells_set(owned_parent_cells.begin(), owned_parent_cells.end());

    // We need to create a list of ghstd and virtual cells
    // We do this by running through sub_cell_indexes
    // and looking at the neighbors to find neighbors that are not owned

    // Virtual cells are defined as cells with indexes outside of the range of owned_cells and ghstd_cells
    const int virtual_parent_cell_offset = num_parent_owned_cells + num_parent_ghstd_cells;

    std::unordered_set<LO> ghstd_parent_cells_set;
    std::unordered_set<LO> virtual_parent_cells_set;
    for(int i=0;i<num_owned_cells;++i){
      const LO parent_cell_index = owned_parent_cells[i];
      for(int local_face_index=0;local_face_index<num_faces_per_cell;++local_face_index){
        const LO parent_face = parent_info.cell_to_faces(parent_cell_index, local_face_index);

        // Sidesets can have owned cells that border the edge of the domain (i.e. parent_face == -1)
        // If we are at the edge of the domain, we can ignore this face.
        if(parent_face < 0)
          continue;

        // Find the side index for neighbor cell with respect to the face
        const LO neighbor_parent_side = (parent_info.face_to_cells(parent_face,0) == parent_cell_index) ? 1 : 0;

        // Get the neighbor cell index in the parent's indexing scheme
        const LO neighbor_parent_cell = parent_info.face_to_cells(parent_face, neighbor_parent_side);

        // If the face exists, then the neighbor should exist
        TEUCHOS_ASSERT(neighbor_parent_cell >= 0);

        // We can easily check if this is a virtual cell
        if(neighbor_parent_cell >= virtual_parent_cell_offset){
          virtual_parent_cells_set.insert(neighbor_parent_cell);
        } else if(neighbor_parent_cell >= num_parent_owned_cells){
          // This is a quick check for a ghost cell
          // This branch only exists to cut down on the number of times the next branch (much slower) is called
          ghstd_parent_cells_set.insert(neighbor_parent_cell);
        } else {
          // There is still potential for this to be a ghost cell with respect to 'our' cells
          // The only way to check this is with a super slow lookup call
          if(owned_parent_cells_set.find(neighbor_parent_cell) == owned_parent_cells_set.end()){
            // The neighbor cell is not owned by 'us', therefore it is a ghost
            ghstd_parent_cells_set.insert(neighbor_parent_cell);
          }
        }
      }
    }

    // We now have a list of the owned, ghstd, and virtual cells in the parent's indexing scheme.
    // We will take the 'unordered_set's ordering for the the sub-indexing scheme

    ghstd_parent_cells.assign(ghstd_parent_cells_set.begin(), ghstd_parent_cells_set.end());
    virtual_parent_cells.assign(virtual_parent_cells_set.begin(), virtual_parent_cells_set.end());

  }

  const int num_ghstd_cells = ghstd_parent_cells.size();
  const int num_virtual_cells = virtual_parent_cells.size();
  const int num_real_cells = num_owned_cells + num_ghstd_cells;
  const int num_total_cells = num_real_cells + num_virtual_cells;

  std::vector<std::pair<LO, LO> > all_parent_cells(num_total_cells);
  for (std::size_t i=0; i< owned_parent_cells.size(); ++i)
    all_parent_cells[i] = std::pair<LO, LO>(owned_parent_cells[i], i);

  for (std::size_t i=0; i< ghstd_parent_cells.size(); ++i) {
    LO insert = owned_parent_cells.size()+i;
    all_parent_cells[insert] = std::pair<LO, LO>(ghstd_parent_cells[i], insert);
  }

  for (std::size_t i=0; i< virtual_parent_cells.size(); ++i) {
    LO insert = owned_parent_cells.size()+ ghstd_parent_cells.size()+i;
    all_parent_cells[insert] = std::pair<LO, LO>(virtual_parent_cells[i], insert);
  }

  sub_info.num_owned_cells = owned_parent_cells.size();
  sub_info.num_ghstd_cells = ghstd_parent_cells.size();
  sub_info.num_virtual_cells = virtual_parent_cells.size();

  // We now have the indexing order for our sub_info

  // Just as a precaution, make sure the parent_info is setup properly
  TEUCHOS_ASSERT(static_cast<int>(parent_info.cell_vertices.extent(0)) == num_parent_total_cells);
  TEUCHOS_ASSERT(static_cast<int>(parent_info.local_cells.extent(0)) == num_parent_total_cells);
  TEUCHOS_ASSERT(static_cast<int>(parent_info.global_cells.extent(0)) == num_parent_total_cells);

  const int num_vertices_per_cell = parent_info.cell_vertices.extent(1);
  const int num_dims = parent_info.cell_vertices.extent(2);

  // Fill owned, ghstd, and virtual cells: global indexes, local indexes and vertices
  sub_info.global_cells = Kokkos::View<GO*>("global_cells", num_total_cells);
  sub_info.local_cells = Kokkos::View<LO*>("local_cells", num_total_cells);
  sub_info.cell_vertices = Kokkos::View<double***, PHX::Device>("cell_vertices", num_total_cells, num_vertices_per_cell, num_dims);
  for(int cell=0;cell<num_total_cells;++cell){
    const LO parent_cell = all_parent_cells[cell].first;
    sub_info.global_cells(cell) = parent_info.global_cells(parent_cell);
    sub_info.local_cells(cell) = parent_info.local_cells(parent_cell);
    for(int vertex=0;vertex<num_vertices_per_cell;++vertex){
      for(int dim=0;dim<num_dims;++dim){
        sub_info.cell_vertices(cell,vertex,dim) = parent_info.cell_vertices(parent_cell,vertex,dim);
      }
    }
  }

  // Now for the difficult part

  // We need to create a new face indexing scheme from the old face indexing scheme

  // Create an auxiliary list with all cells - note this preserves indexing

  struct face_t{
    face_t(LO c0, LO c1, LO sc0, LO sc1)
    {
      cell_0=c0;
      cell_1=c1;
      subcell_index_0=sc0;
      subcell_index_1=sc1;
    }
    LO cell_0;
    LO cell_1;
    LO subcell_index_0;
    LO subcell_index_1;
  };


  // First create the faces
  std::vector<face_t> faces;
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("Create faces",CreateFaces);
    // faces_set: cell_0, subcell_index_0, cell_1, subcell_index_1
    std::unordered_map<LO,std::unordered_map<LO, std::pair<LO,LO> > > faces_set;

    std::sort(all_parent_cells.begin(), all_parent_cells.end());

    for(int owned_cell=0;owned_cell<num_owned_cells;++owned_cell){
      const LO owned_parent_cell = owned_parent_cells[owned_cell];
      for(int local_face=0;local_face<num_faces_per_cell;++local_face){
        const LO parent_face = parent_info.cell_to_faces(owned_parent_cell,local_face);

        // Skip faces at the edge of the domain
        if(parent_face<0)
          continue;

        // Get the cell on the other side of the face
        const LO neighbor_side = (parent_info.face_to_cells(parent_face,0) == owned_parent_cell) ? 1 : 0;

        const LO neighbor_parent_cell = parent_info.face_to_cells(parent_face, neighbor_side);
        const LO neighbor_subcell_index = parent_info.face_to_lidx(parent_face, neighbor_side);

        // Convert parent cell index into sub cell index
        std::pair<LO, LO> search_point(neighbor_parent_cell, 0);
        auto itr = std::lower_bound(all_parent_cells.begin(), all_parent_cells.end(), search_point);

        TEUCHOS_TEST_FOR_EXCEPT_MSG(itr == all_parent_cells.end(), "panzer_stk::setupSubLocalMeshInfo : Neighbor cell was not found in owned, ghosted, or virtual cells");

        const LO neighbor_cell = itr->second;

        LO cell_0, cell_1, subcell_index_0, subcell_index_1;
        if(owned_cell < neighbor_cell){
          cell_0 = owned_cell;
          subcell_index_0 = local_face;
          cell_1 = neighbor_cell;
          subcell_index_1 = neighbor_subcell_index;
        } else {
          cell_1 = owned_cell;
          subcell_index_1 = local_face;
          cell_0 = neighbor_cell;
          subcell_index_0 = neighbor_subcell_index;
        }

        // Add this interface to the set of faces - smaller cell index is 'left' (or '0') side of face
        faces_set[cell_0][subcell_index_0] = std::pair<LO,LO>(cell_1, subcell_index_1);
      }
    }

    for(const auto & cell_pair : faces_set){
      const LO cell_0 = cell_pair.first;
      for(const auto & subcell_pair : cell_pair.second){
        const LO subcell_index_0 = subcell_pair.first;
        const LO cell_1 = subcell_pair.second.first;
        const LO subcell_index_1 = subcell_pair.second.second;
        faces.push_back(face_t(cell_0,cell_1,subcell_index_0,subcell_index_1));
      }
    }
  }

  const int num_faces = faces.size();

  sub_info.face_to_cells = Kokkos::View<LO*[2]>("face_to_cells", num_faces);
  sub_info.face_to_lidx = Kokkos::View<LO*[2]>("face_to_lidx", num_faces);
  sub_info.cell_to_faces = Kokkos::View<LO**>("cell_to_faces", num_total_cells, num_faces_per_cell);

  // Default the system with invalid cell index
  Kokkos::deep_copy(sub_info.cell_to_faces, -1);

  for(int face_index=0;face_index<num_faces;++face_index){
    const face_t & face = faces[face_index];

    sub_info.face_to_cells(face_index,0) = face.cell_0;
    sub_info.face_to_cells(face_index,1) = face.cell_1;

    sub_info.cell_to_faces(face.cell_0,face.subcell_index_0) = face_index;
    sub_info.cell_to_faces(face.cell_1,face.subcell_index_1) = face_index;

    sub_info.face_to_lidx(face_index,0) = face.subcell_index_0;
    sub_info.face_to_lidx(face_index,1) = face.subcell_index_1;

  }

  // Complete.

}





template<typename LO, typename GO>
void
splitMeshInfo(const panzer::LocalMeshInfoBase<LO,GO> & mesh_info,
              const int splitting_size,
              std::vector<panzer::LocalMeshPartition<LO,GO> > & partitions)
{

  // Make sure the splitting size makes sense
  TEUCHOS_ASSERT(splitting_size != 0);

  // This is not a partitioning scheme.
  // This just breaks the mesh_info into equally sized chunks and ignores connectivity
  // This means that the cells in the partition probably won't be nearby each other - leads to an excess of ghost cells

  const LO num_owned_cells = mesh_info.num_owned_cells;

  if(splitting_size < 0){

    // Just one chunk
    std::vector<LO> partition_cells;
    partition_cells.reserve(mesh_info.num_owned_cells);
    for(LO i=0;i<mesh_info.num_owned_cells;++i){
      partition_cells.push_back(i);
    }

    // It seems this can happen
    if(partition_cells.size() > 0){
      partitions.push_back(panzer::LocalMeshPartition<LO,GO>());
      setupSubLocalMeshInfo(mesh_info,partition_cells,partitions.back());
    }

  } else {

    std::vector<LO> partition_cells;
    partition_cells.reserve(splitting_size);

    // There should be at least one partition
    const LO num_partitions = mesh_info.num_owned_cells / splitting_size + ((mesh_info.num_owned_cells % splitting_size == 0) ? 0 : 1);

    LO partition_start_index = 0;
    for(LO partition=0;partition<num_partitions;++partition){

      // Make sure end of partition maxes out at num_owned_cells
      const LO partition_end_index = std::min(partition_start_index + splitting_size, num_owned_cells);

      partition_cells.resize(partition_end_index - partition_start_index,-1);
      for(int i=partition_start_index;i<partition_end_index;++i){
        partition_cells[i] = i;
      }

      // It seems this can happen
      if(partition_cells.size() > 0){
        partitions.push_back(panzer::LocalMeshPartition<LO,GO>());
        setupSubLocalMeshInfo(mesh_info,partition_cells,partitions.back());
      }

      partition_start_index = partition_end_index;
    }

  }

}

template<typename LO, typename GO>
void
partitionMeshInfo(const panzer::LocalMeshInfoBase<LO,GO>& /* mesh_info */,
                 const size_t /* requested_partition_size */,
                 std::vector<panzer::LocalMeshPartition<LO,GO>>& /* partitions */)
{
  // Not yet sure how to do this
  TEUCHOS_ASSERT(false);
}

}

template <typename LO, typename GO>
void
generateLocalMeshPartitions(const panzer::LocalMeshInfo<LO,GO> & mesh_info,
                            const panzer::WorksetDescriptor & description,
                            std::vector<panzer::LocalMeshPartition<LO,GO> > & partitions)
{

  // We have to make sure that the partitioning is possible
  TEUCHOS_ASSERT(description.getWorksetSize() != panzer::WorksetSizeType::CLASSIC_MODE);
  TEUCHOS_ASSERT(description.getWorksetSize() != 0);

  // This could just return, but it would be difficult to debug why no partitions were returned
  TEUCHOS_ASSERT(description.requiresPartitioning());

  const std::string & element_block_name = description.getElementBlock();

  // We have two processes for in case this is a sideset or element block
  if(description.useSideset()){

    // If the element block doesn't exist, there are no partitions to create
    if(mesh_info.sidesets.find(element_block_name) == mesh_info.sidesets.end())
      return;
    const auto & sideset_map = mesh_info.sidesets.at(element_block_name);

    const std::string & sideset_name = description.getSideset();

    // If the sideset doesn't exist, there are no partitions to create
    if(sideset_map.find(sideset_name) == sideset_map.end())
      return;

    const panzer::LocalMeshSidesetInfo<LO,GO> & sideset_info = sideset_map.at(sideset_name);

    // Partitioning is not important for sidesets
    panzer::partitioning_utilities::splitMeshInfo<LO,GO>(sideset_info, description.getWorksetSize(), partitions);

    for(auto & partition : partitions){
      partition.sideset_name = sideset_name;
      partition.element_block_name = element_block_name;
      partition.cell_topology = sideset_info.cell_topology;
    }

  } else {

    // If the element block doesn't exist, there are no partitions to create
    if(mesh_info.element_blocks.find(element_block_name) == mesh_info.element_blocks.end())
      return;

    // Grab the element block we're interested in
    const panzer::LocalMeshBlockInfo<LO,GO> & block_info = mesh_info.element_blocks.at(element_block_name);

    if(description.getWorksetSize() == panzer::WorksetSizeType::ALL_ELEMENTS){
      // We only have one partition describing the entire local mesh
      panzer::partitioning_utilities::splitMeshInfo(block_info, -1, partitions);
    } else {
      // We need to partition local mesh

      // TODO: This needs to be used, but first needs to be written
      //panzer::partitioning_utilities::partitionMeshInfo<LO,GO>(block_info, description.getWorksetSize(), partitions);

      // FIXME: Until the above function is fixed, we will use this hack - this will lead to horrible partitions
      panzer::partitioning_utilities::splitMeshInfo(block_info, description.getWorksetSize(), partitions);

    }

    for(auto & partition : partitions){
      partition.element_block_name = element_block_name;
      partition.cell_topology = block_info.cell_topology;
    }
  }

}

}


template
void
panzer::partitioning_utilities::setupSubLocalMeshInfo<int,int>(const panzer::LocalMeshInfoBase<int,int> & parent_info,
                                                               const std::vector<int> & owned_parent_cells,
                                                               panzer::LocalMeshInfoBase<int,int> & sub_info);

#ifndef PANZER_ORDINAL64_IS_INT
template
void
panzer::partitioning_utilities::setupSubLocalMeshInfo<int,panzer::Ordinal64>(const panzer::LocalMeshInfoBase<int,panzer::Ordinal64> & parent_info,
                                                                             const std::vector<int> & owned_parent_cells,
                                                                             panzer::LocalMeshInfoBase<int,panzer::Ordinal64> & sub_info);
#endif

template
void
panzer::generateLocalMeshPartitions<int,int>(const panzer::LocalMeshInfo<int,int> & mesh_info,
                                             const panzer::WorksetDescriptor & description,
                                             std::vector<panzer::LocalMeshPartition<int,int> > & partitions);

#ifndef PANZER_ORDINAL64_IS_INT
template
void
panzer::generateLocalMeshPartitions<int,panzer::Ordinal64>(const panzer::LocalMeshInfo<int,panzer::Ordinal64> & mesh_info,
                                                           const panzer::WorksetDescriptor & description,
                                                           std::vector<panzer::LocalMeshPartition<int,panzer::Ordinal64> > & partitions);
#endif
