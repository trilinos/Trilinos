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

#include "Panzer_WorksetUtilities.hpp"

#include "Panzer_LocalPartitioningUtilities.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Kokkos_ViewFactory.hpp"

#include <set>

namespace panzer
{

namespace
{
Teuchos::RCP<std::vector<panzer::Workset> >
buildGroupedSubcellWorksets(const panzer::LocalMeshInfo & mesh_info,
                            const std::string & element_block_a,
                            const std::string & element_block_b,
                            const std::string & sideset,
                            const int workset_size,
                            const Teuchos::RCP<const OrientationsInterface> & orientations)
{

  // Cannot build side assemblies without connectivity
  TEUCHOS_TEST_FOR_EXCEPT_MSG(not mesh_info.has_connectivity,
                              "panzer::BuildSideAssembledWorksets : Attempted to build a side assembled workset from a LocalMeshInfo object without connectivity enabled");


  // ===================================================
  // Get sideset a info

  // If the element block is not local, then we ignore it
  const auto element_block_itr = mesh_info.sidesets.find(element_block_a);
  if(element_block_itr == mesh_info.sidesets.end())
    return Teuchos::rcp(new std::vector<panzer::Workset>());

  const auto & sideset_map = element_block_itr->second;
  const auto sideset_itr = sideset_map.find(sideset);

  // If the sideset is not local for this element block, then we ignore it
  if(sideset_itr == sideset_map.end())
    return Teuchos::rcp(new std::vector<panzer::Workset>());

  const auto & sideset_info = sideset_itr->second;

  // Used for looking up element block ownership for cells
  const auto & cell_sets = *sideset_info.cell_sets;

  // Side assembly requires information about cell block ownership
  TEUCHOS_TEST_FOR_EXCEPT_MSG(cell_sets.getNumElements() != sideset_info.num_owned_cells + sideset_info.num_ghost_cells + sideset_info.num_virtual_cells,
                              "panzer::BuildSideAssembledWorksets : Attempted to build a side assembled workset and either (1) don't have any cells in the LocalMeshInfo or (2) forgot to add cell block ownership to the LocalMeshInfo.");

  // ===================================================
  // We now have a sideset info object, from which we need to find all faces joining block a to block b

  const size_t num_total_faces = sideset_info.face_to_cells.extent(0);

  // Make sure element block 'b' is valid - we will return no worksets if no block is found
//  TEUCHOS_TEST_FOR_EXCEPT_MSG(sideset_info.block_to_id_map.find(element_block_b) == sideset_info.block_to_id_map.end(),
//                              "panzer::BuildSideAssembledWorksets : Could not find element block '"<<element_block_b<<"' for sideset '"<<sideset<<"'.");
  if(not cell_sets.hasSet(element_block_b))
    return Teuchos::rcp(new std::vector<panzer::Workset>());

  const unsigned int element_block_b_id = cell_sets.getSetIndex(element_block_b);

  std::vector<std::size_t> cell_ids_a, cell_ids_b;
  std::vector<std::size_t> local_cell_ids_a, local_cell_ids_b;
  std::vector<std::size_t> local_side_ids_a, local_side_ids_b;
  for(size_t face=0; face<num_total_faces; ++face){

    // Get the external cell - 1 is always the external face for a sideset info object
    const size_t cell_b = sideset_info.face_to_cells(face,1);

    // If cell b is not in block b, then skip this face
    if(element_block_b_id != cell_sets.getElementSetIndex(cell_b))
      continue;

    const size_t cell_a = sideset_info.face_to_cells(face,0);

    cell_ids_a.push_back(cell_a);
    cell_ids_b.push_back(cell_b);

    local_cell_ids_a.push_back((size_t) sideset_info.local_cells(cell_a));
    local_cell_ids_b.push_back((size_t) sideset_info.local_cells(cell_b));

    // TODO: face_to_lidx uses -1 to represent a missing local index...
    local_side_ids_a.push_back((size_t) sideset_info.face_to_lidx(face,0));
    local_side_ids_b.push_back((size_t) sideset_info.face_to_lidx(face,1));

  }

  const size_t num_cells = local_side_ids_a.size();

  // If we didn't find any cells, return nothing
  if(num_cells == 0)
    return Teuchos::rcp(new std::vector<panzer::Workset>());

  // Now that we know the number of cells, we can copy over the vertices
  const auto & vertices = sideset_info.cell_vertices;
  const size_t num_vertices = vertices.extent(1);
  const size_t num_dims = vertices.extent(2);

  // Now we know how many cells exist, so we can allocate the cell vertices
  Kokkos::DynRankView<double,PHX::Device> vertex_coordinates_a, vertex_coordinates_b;
  vertex_coordinates_a = Kokkos::createDynRankView(vertex_coordinates_a, "vertex coordinates a", num_cells, num_vertices, num_dims);
  vertex_coordinates_b = Kokkos::createDynRankView(vertex_coordinates_b, "vertex coordinates b", num_cells, num_vertices, num_dims);

  // Fill the vertices
  for(size_t cell=0;cell<num_cells;++cell){

    const size_t cell_a = cell_ids_a[cell];
    for(int v=0; v<num_vertices; ++v)
      for(int dim=0; dim<num_dims; ++dim)
        vertex_coordinates_a(cell,v,dim) = vertices(cell_a,v,dim);

    const size_t cell_b = cell_ids_b[cell];
    for(int v=0; v<num_vertices; ++v)
      for(int dim=0; dim<num_dims; ++dim)
        vertex_coordinates_b(cell,v,dim) = vertices(cell_b,v,dim);

  }

  return panzer::buildGroupedSubcellWorksets(element_block_a,
                                             sideset,
                                             sideset_info.cell_topology,
                                             local_cell_ids_a,
                                             local_side_ids_a,
                                             vertex_coordinates_a,
                                             element_block_b,
                                             sideset,
                                             sideset_info.cell_topology,
                                             local_cell_ids_b,
                                             local_side_ids_b,
                                             vertex_coordinates_b,
                                             workset_size,
                                             orientations);

}


Teuchos::RCP<std::vector<panzer::Workset> >
buildGroupedSubcellWorksets(const panzer::LocalMeshInfo & mesh_info,
                            const std::string & element_block,
                            const std::string & sideset,
                            const bool force_side_assembly,
                            const int workset_size,
                            const Teuchos::RCP<const OrientationsInterface> & orientations)
{

  const auto element_block_itr = mesh_info.sidesets.find(element_block);

  // If the element block is not local, then we ignore it
  if(element_block_itr == mesh_info.sidesets.end())
    return Teuchos::rcp(new std::vector<panzer::Workset>());

  const auto & sideset_map = element_block_itr->second;
  const auto sideset_itr = sideset_map.find(sideset);

  // If the sideset is not local for this element block, then we ignore it
  if(sideset_itr == sideset_map.end())
    return Teuchos::rcp(new std::vector<panzer::Workset>());

  const auto & sideset_info = sideset_itr->second;
  const auto & sideset_coords = sideset_info.cell_vertices;

  // If there are no owned cells in the sideset_info, then ignore it
  if(sideset_info.num_owned_cells == 0)
    return Teuchos::rcp(new std::vector<panzer::Workset>());

  const int num_vertices = sideset_coords.extent_int(1);
  const int num_dims = sideset_coords.extent_int(2);

  const size_t num_faces = sideset_info.face_to_cells.extent(0);
  std::vector<size_t> cell_ids, local_side_ids;
  cell_ids.reserve(num_faces);
  local_side_ids.reserve(num_faces);
  for(size_t face=0; face<num_faces; ++face){

    // We only want the 0 cell
    size_t cell = sideset_info.face_to_cells(face,0);
    int lidx = sideset_info.face_to_lidx(face,0);

    // Sanity check - local face index 0 should always be an owned cell and should never be -1
    TEUCHOS_ASSERT(cell < sideset_info.num_owned_cells);

    // The local face index cannot be negative for an existing face
    TEUCHOS_ASSERT(lidx >= 0);

    cell_ids.push_back(cell);
    local_side_ids.push_back(lidx);

  }

  const size_t num_cells = cell_ids.size();
  Kokkos::DynRankView<double,PHX::Device> vertex_coordinates;
  vertex_coordinates = Kokkos::createDynRankView(vertex_coordinates, "vertex coordinates", num_cells, num_vertices, num_dims);

  // Now we copy over vertex coordinates and convert cell_ids to local_cell_ids
  for(size_t i=0; i<num_cells; ++i){

    const size_t cell = cell_ids[i];

    // Store the local cell id
    cell_ids[i] = sideset_info.local_cells(cell);

    for(int vertex=0; vertex<num_vertices; ++vertex)
      for(int dim=0; dim<num_dims; ++dim)
        vertex_coordinates(i, vertex, dim) = sideset_coords(cell,vertex,dim);
  }

  // We have enough information now to create the worksets
  return panzer::buildGroupedSubcellWorksets(element_block,
                                             sideset,
                                             sideset_info.cell_topology,
                                             cell_ids,
                                             local_side_ids,
                                             vertex_coordinates,
                                             force_side_assembly,
                                             workset_size,
                                             orientations);
}

}

Teuchos::RCP<std::vector<panzer::Workset> >
buildWorksets(const panzer::LocalMeshInfo & mesh_info,
              const panzer::WorksetDescriptor & description,
              const Teuchos::RCP<const OrientationsInterface> & orientations)
{

  // There are four options for worksets through this interface

  // 1) Volume partitions for an element block
  // 2) Volume partitions for an element block + sideset
  // 3) Surface partitions for an element block + sideset (side assembly)
  // 4) Volume partitions for an element block + element_block + sideset (boundary between two element blocks) (side assembly)

  // Option 1 and 2 are supported by the partition interface, but 3 and 4 require the older interface
  if(description.useSideset()){
    if(description.connectsElementBlocks()){

      // This would require a major overhaul of LocalMeshInfo
      TEUCHOS_ASSERT(description.getSideset(0) == description.getSideset(1));

      // Option 4
      return buildGroupedSubcellWorksets(mesh_info,
                                        description.getElementBlock(0), description.getElementBlock(1),
                                        description.getSideset(0),
                                        description.getWorksetSize(),
                                        orientations);

    }

    if(description.groupBySubcell()){

      // Option 3
      return buildGroupedSubcellWorksets(mesh_info,
                                         description.getElementBlock(), description.getSideset(),
                                         description.sideAssembly(),
                                         description.getWorksetSize(),
                                         orientations);

    }
  }


  Teuchos::RCP<std::vector<panzer::Workset> > worksets = Teuchos::rcp(new std::vector<panzer::Workset>());

  // Each partition represents a chunk of the local mesh
  std::vector<panzer::LocalMeshPartition > partitions;
  panzer::generateLocalMeshPartitions(mesh_info, description, partitions);

  WorksetOptions options;
  options.side_assembly_ = description.sideAssembly();
  options.align_side_points_ = options.side_assembly_ && description.connectsElementBlocks();
  options.orientations_ = orientations;

  // Each partition generates a single workset
  for(const auto & partition : partitions){
    worksets->push_back(panzer::Workset());
    worksets->back().setup(partition, options);
  }

  return worksets;

}


Teuchos::RCP<std::vector<Workset> >
buildGroupedSubcellWorksets(const std::string & element_block,
                            const std::string & sideset,
                            const Teuchos::RCP<const shards::CellTopology> & cell_topology,
                            const std::vector<std::size_t>& local_cell_ids,
                            const std::vector<std::size_t>& local_side_ids,
                            const Kokkos::DynRankView<double,PHX::Device> & vertex_coordinates,
                            const bool force_side_assembly,
                            const int workset_size,
                            const Teuchos::RCP<const OrientationsInterface> & orientations)
{
  const std::size_t num_total_cells = local_cell_ids.size();

  TEUCHOS_ASSERT(local_side_ids.size() == num_total_cells);
  TEUCHOS_ASSERT(static_cast<std::size_t>(vertex_coordinates.extent(0)) == num_total_cells);

  // Build a single interface object for all edges (it's mostly boilerplate)
  LocalMeshPartition info;
  info.element_block_name = element_block;
  info.sideset_name = sideset;
  info.cell_topology = cell_topology;
  info.subcell_dimension = cell_topology->getDimension()-1;

  WorksetOptions options;
  options.side_assembly_ = force_side_assembly;
  options.align_side_points_ = false;
  options.orientations_ = orientations;

  auto worksets_ptr = Teuchos::rcp(new std::vector<panzer::Workset>);
  auto & worksets = *worksets_ptr;

  // If there are not cells, then don't bother generating worksets
  if(num_total_cells==0)
    return Teuchos::rcp(new std::vector<Workset>);

  // Create a lookup map for each set of elements that share a subcell index
  std::map<unsigned int, std::vector<size_t>> subcell_map;
  for(size_t i=0; i<num_total_cells; ++i)
    subcell_map[local_side_ids[i]].push_back(i);

  unsigned int num_vertices_per_cell = vertex_coordinates.extent(1);
  unsigned int num_dims = vertex_coordinates.extent(2);

  for(const auto & pr : subcell_map){

    info.subcell_index = pr.first;

    const size_t num_subcell_cells = pr.second.size();
    const std::size_t work_size = (workset_size <= 0) ? num_subcell_cells : workset_size;

    // fill the worksets
    size_t cell_count = 0;
    while(cell_count < num_subcell_cells){

      size_t this_work_size = work_size;
      if(cell_count + this_work_size > num_subcell_cells)
        this_work_size = num_subcell_cells - cell_count;

      Kokkos::View<int*,PHX::Device> local_cells("local_cells", this_work_size);
      Kokkos::View<double***,PHX::Device> cell_vertices("cell_vertices",this_work_size,num_vertices_per_cell,num_dims);

      for(size_t i=0; i<this_work_size; ++i){
        const size_t cell = cell_count + i;
        local_cells(i) = local_cell_ids[cell];
        for(size_t j=0; j<num_vertices_per_cell; ++j)
          for(unsigned int dim=0; dim<num_dims; ++dim)
            cell_vertices(i,j,dim) = vertex_coordinates(cell,j,dim);
      }

      // Setup the info object for this edge
      info.local_cells = local_cells;
      info.cell_vertices = cell_vertices;
      info.num_owned_cells = this_work_size;

      // Create and setup the new workset
      panzer::Workset workset;
      workset.setup(info, options);

      // Add to list of worksets
      worksets.push_back(workset);

      // Increment cell count
      cell_count += this_work_size;

    }
  }

  return worksets_ptr;
}

Teuchos::RCP<std::vector<Workset> >
buildGroupedSubcellWorksets(const std::string & element_block_a,
                            const std::string & sideset_a,
                            const Teuchos::RCP<const shards::CellTopology> & cell_topology_a,
                            const std::vector<std::size_t>& local_cell_ids_a,
                            const std::vector<std::size_t>& local_side_ids_a,
                            const Kokkos::DynRankView<double,PHX::Device> & vertex_coordinates_a,
                            const std::string & element_block_b,
                            const std::string & sideset_b,
                            const Teuchos::RCP<const shards::CellTopology> & cell_topology_b,
                            const std::vector<std::size_t>& local_cell_ids_b,
                            const std::vector<std::size_t>& local_side_ids_b,
                            const Kokkos::DynRankView<double,PHX::Device> & vertex_coordinates_b,
                            const int workset_size,
                            const Teuchos::RCP<const OrientationsInterface> & orientations)
{

  const std::size_t total_num_cells = local_cell_ids_a.size();
  const std::size_t total_num_cells_b = local_cell_ids_b.size();

  TEUCHOS_ASSERT(total_num_cells==total_num_cells_b);
  TEUCHOS_ASSERT(local_side_ids_a.size() == local_cell_ids_a.size());
  TEUCHOS_ASSERT(local_side_ids_a.size() == static_cast<std::size_t>(vertex_coordinates_a.extent(0)));
  TEUCHOS_ASSERT(local_side_ids_b.size() == local_cell_ids_b.size());
  TEUCHOS_ASSERT(local_side_ids_b.size() == static_cast<std::size_t>(vertex_coordinates_b.extent(0)));

  std::vector<LocalMeshPartition> zones;
  zones.resize(2);

  zones[0].subcell_dimension = zones[1].subcell_dimension = cell_topology_a->getDimension()-1;

  zones[0].element_block_name = element_block_a;
  zones[1].element_block_name = element_block_b;

  zones[0].sideset_name = sideset_a;
  zones[1].sideset_name = sideset_b;

  zones[0].cell_topology = cell_topology_a;
  zones[1].cell_topology = cell_topology_b;

  WorksetOptions options;
  options.side_assembly_ = true;
  options.align_side_points_ = true;
  options.orientations_ = orientations;

  auto worksets_ptr = Teuchos::rcp(new std::vector<panzer::Workset>);
  auto & worksets = *worksets_ptr;

  // If there are no cells, then we don't need a workset
  if(total_num_cells==0)
    return Teuchos::rcp(new std::vector<Workset>());

  // This collects all the elements that share the same sub cell pairs, this makes it easier to
  // build the required worksets
  // key is the pair of local face indices, value is a vector of cell indices that satisfy this pair
  std::map<std::pair<unsigned,unsigned>,std::vector<std::size_t> > element_list;
  for (std::size_t cell=0; cell < local_cell_ids_a.size(); ++cell)
    element_list[std::make_pair(local_side_ids_a[cell],local_side_ids_b[cell])].push_back(cell);

  unsigned int num_vertices_per_cell_a = vertex_coordinates_a.extent(1);
  unsigned int num_vertices_per_cell_b = vertex_coordinates_b.extent(1);

  unsigned int num_dims = vertex_coordinates_a.extent(2);

  // fill the worksets
  for(auto edge=element_list.begin(); edge!=element_list.end();++edge) {
    // loop over each workset
    const std::vector<std::size_t> & cell_indices = edge->second;

    size_t num_edge_cells = cell_indices.size();

    const std::size_t work_size = (workset_size <= 0) ? num_edge_cells : workset_size;

    size_t cell_count = 0;
    while(cell_count < num_edge_cells){

      size_t this_work_size = work_size;
      if(cell_count + this_work_size > num_edge_cells)
        this_work_size = num_edge_cells - cell_count;

      Kokkos::View<int*,PHX::Device> local_cells_a("local_cells_a", this_work_size);
      Kokkos::View<int*,PHX::Device> local_cells_b("local_cells_b", this_work_size);

      Kokkos::View<double***,PHX::Device> cell_vertices_a("cell_vertices_a",this_work_size,num_vertices_per_cell_a,num_dims);
      Kokkos::View<double***,PHX::Device> cell_vertices_b("cell_vertices_b",this_work_size,num_vertices_per_cell_b,num_dims);

      for(size_t i=0; i<this_work_size; ++i){
        const size_t cell = cell_indices[cell_count + i];
        local_cells_a(i) = local_cell_ids_a[cell];
        local_cells_b(i) = local_cell_ids_b[cell];
        for(size_t j=0; j<num_vertices_per_cell_a; ++j)
          for(unsigned int dim=0; dim<num_dims; ++dim)
            cell_vertices_a(i,j,dim) = vertex_coordinates_a(cell,j,dim);
        for(size_t j=0; j<num_vertices_per_cell_b; ++j)
          for(unsigned int dim=0; dim<num_dims; ++dim)
            cell_vertices_b(i,j,dim) = vertex_coordinates_b(cell,j,dim);

      }

      // Setup the info object for this edge
      zones[0].num_owned_cells = this_work_size;
      zones[1].num_owned_cells = this_work_size;

      zones[0].subcell_index = edge->first.first;
      zones[1].subcell_index = edge->first.second;

      zones[0].local_cells = local_cells_a;
      zones[1].local_cells = local_cells_b;

      zones[0].cell_vertices = cell_vertices_a;
      zones[1].cell_vertices = cell_vertices_b;

      // Create and setup the new workset
      panzer::Workset workset;
      workset.setup(zones, options);

      // Add to list of worksets
      worksets.push_back(workset);

      // Increment cell count
      cell_count += this_work_size;

    }
  }

  return worksets_ptr;
}

}
