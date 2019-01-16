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

#include "Panzer_NodeType.hpp"
#include "Panzer_STK_LocalMeshUtilities.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"

//#include "Panzer_HashUtils.hpp"
#include "Panzer_LocalMeshInfo.hpp"
#include "Panzer_LocalPartitioningUtilities.hpp"
#include "Panzer_FaceToElement.hpp"

#include "Panzer_ConnManager.hpp"
#include "Panzer_FieldPattern.hpp"
#include "Panzer_NodalFieldPattern.hpp"
#include "Panzer_EdgeFieldPattern.hpp"
#include "Panzer_FaceFieldPattern.hpp"
#include "Panzer_ElemFieldPattern.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_OrdinalTraits.hpp"

//#include "Tpetra_Import.hpp"
//#include "Tpetra_CrsMatrix.hpp"
//#include "Tpetra_RowMatrixTransposer.hpp"

#include <string>
#include <map>
#include <vector>
#include <unordered_set>

namespace panzer_stk
{

// No external access
namespace
{

void
setupLocalMeshSidesetInfo(const panzer_stk::STK_Interface & mesh,
                          panzer::ConnManager&  conn ,
                          const panzer::LocalMeshInfo & mesh_info,
                          const std::string & element_block,
                          const std::string & sideset,
                          std::map<std::string, std::map<std::string, panzer::LocalMeshSidesetInfo > > & sideset_map)
{

  // This function identifies all cells in mesh_info that belong to element_block
  // and creates a block_info from it.

  const size_t face_subcell_dimension = mesh.getCellTopology(element_block)->getDimension()-1;

  Kokkos::DynRankView<double,PHX::Device> data_allocation;

  // This is a list of all entities that make up the 'side'
  std::vector<stk::mesh::Entity> side_entities;
  {

    // Note we don't use mesh.getMySides because is has the additional restriction of only allowing 'owned' sides,
    // which doesn't work for us because all elements/faces/edges/nodes are specifically owned by processes
    stk::mesh::Part * side_part = mesh.getSideset(sideset);
    stk::mesh::Part * elmt_part = mesh.getElementBlockPart(element_block);
    TEUCHOS_ASSERT(side_part != 0);
    TEUCHOS_ASSERT(elmt_part != 0);

    stk::mesh::Selector side = *side_part;
    stk::mesh::Selector block = *elmt_part;

    // grab sides that exist on the sideset and the element block
    // Note that this will give us sides between the ghost cells and the sideset
    stk::mesh::get_selected_entities(side & block,mesh.getBulkData()->buckets(mesh.getSideRank()),side_entities);

  }

  // We now have a list of sideset entities, lets unwrap them and create the sideset_info!

  // Start by getting the local cell id for each face on the sideset - indexed by subcell index (local face id)
  std::unordered_map<size_t,std::vector<size_t> > owned_parent_cell_index_map;
  {

    // Create a mapping from local cell ids to cell indexes
    std::unordered_map<size_t, size_t> local_id_map;
    for(unsigned int i=0; i<mesh_info.num_owned_cells + mesh_info.num_ghost_cells; ++i)
      local_id_map[mesh_info.local_cells(i)] = i;

    // Get owned cells associated with faces
    std::vector<stk::mesh::Entity> owned_elements;
    std::vector<size_t> subcell_indexes;
    panzer_stk::workset_utils::getSubcellElements(mesh, element_block, side_entities, subcell_indexes, owned_elements);

    // Create a map of owned_parent cells and their associated subcell indexes that make up the sideset
    for(size_t i=0; i<subcell_indexes.size(); ++i)
      owned_parent_cell_index_map[mesh.elementLocalId(owned_elements[i])].push_back(subcell_indexes[i]);
  }

  const size_t num_owned_cells = owned_parent_cell_index_map.size();

  // If there are no cells, then don't add a sideset info object
  if(num_owned_cells == 0)
    return;

  // Create a new info object for the sideset
  auto & sideset_info = sideset_map[element_block][sideset];

  sideset_info.element_block_name = element_block;
  sideset_info.sideset_name = sideset;
  sideset_info.cell_topology = mesh.getCellTopology(element_block);

  sideset_info.num_owned_cells = num_owned_cells;

  struct Face{
    Face(panzer::LocalOrdinal c0, panzer::LocalOrdinal c1, panzer::LocalOrdinal sc0, panzer::LocalOrdinal sc1)
    {
      cell_0=c0;
      cell_1=c1;
      subcell_index_0=sc0;
      subcell_index_1=sc1;
    }
    panzer::LocalOrdinal cell_0;
    panzer::LocalOrdinal cell_1;
    panzer::LocalOrdinal subcell_index_0;
    panzer::LocalOrdinal subcell_index_1;
  };

  // Figure out how many cells on the other side of the sideset are ghost or virtual
  std::set<size_t> ghstd_parent_cells_set, virtual_parent_cells_set;
  std::vector<Face> faces;
  {
    const size_t parent_virtual_cell_offset = mesh_info.num_owned_cells + mesh_info.num_ghost_cells;
    for(const auto & pr : owned_parent_cell_index_map){
      const auto & parent_cell = pr.first;
      const auto & subcell_indexes = pr.second;

      for(const auto & subcell_index : subcell_indexes){

        // Get the face associated with this subcell index
        const auto face = mesh_info.cell_to_faces(parent_cell, subcell_index);

        const auto other_side = (mesh_info.face_to_cells(face,0) == parent_cell) ? 1 : 0;

        TEUCHOS_ASSERT(subcell_index == mesh_info.face_to_lidx(face, 1-other_side));

        // The cell on the other side of this subcell index is needed to define the face
        const auto other_side_cell = mesh_info.face_to_cells(face, other_side);
        const auto other_side_subcell_index = mesh_info.face_to_lidx(face, other_side);

        faces.push_back(Face(parent_cell, other_side_cell, subcell_index, other_side_subcell_index));

        if(other_side_cell >= parent_virtual_cell_offset){
          virtual_parent_cells_set.insert(other_side_cell);
        } else {
          ghstd_parent_cells_set.insert(other_side_cell);
        }
      }
    }
  }

  // We construct this vector to define the indexing for the cells in this sideset
  std::vector<size_t> all_parent_cells;
  for(const auto & pr : owned_parent_cell_index_map)
    all_parent_cells.push_back(pr.first);
  all_parent_cells.insert(all_parent_cells.end(),ghstd_parent_cells_set.begin(),ghstd_parent_cells_set.end());
  all_parent_cells.insert(all_parent_cells.end(),virtual_parent_cells_set.begin(),virtual_parent_cells_set.end());

  const size_t num_total_cells = all_parent_cells.size();
  const size_t num_faces = faces.size();
  const unsigned int num_vertices_per_cell = mesh_info.cell_vertices.extent(1);
  const unsigned int num_dims_per_vertex = mesh_info.cell_vertices.extent(2);
  const unsigned int num_faces_per_cell = mesh_info.cell_to_faces.extent(1);

  sideset_info.num_ghost_cells = ghstd_parent_cells_set.size();
  sideset_info.num_virtual_cells = virtual_parent_cells_set.size();

  sideset_info.local_cells = Kokkos::View<panzer::LocalOrdinal*,PHX::Device>("local_cell_ids",num_total_cells);
  sideset_info.global_cells = Kokkos::View<panzer::GlobalOrdinal*,PHX::Device>("global_cell_ids",num_total_cells);

  sideset_info.cell_vertices = Kokkos::View<double***,PHX::Device>("cell_vertices",num_total_cells,num_vertices_per_cell,num_dims_per_vertex);


  // This interface is half-implemented so we need to only use it if it exists
  if(mesh_info.cell_sets->getNumElements() > 0){
    std::vector<int> sub_cells(num_total_cells);
    for(size_t i=0; i<num_total_cells; ++i)
      sub_cells[i] = all_parent_cells[i];

    sideset_info.cell_sets = mesh_info.cell_sets->buildSubsets(sub_cells);
  }

  // Copy basic stuff from parent to sideset info objects
  for(size_t i=0; i<num_total_cells; ++i){
    const size_t parent_cell = all_parent_cells[i];
    sideset_info.local_cells(i) = mesh_info.local_cells(parent_cell);
    sideset_info.global_cells(i) = mesh_info.global_cells(parent_cell);
    for(unsigned int j=0; j<num_vertices_per_cell; ++j)
      for(unsigned int k=0; k<num_dims_per_vertex; ++k)
        sideset_info.cell_vertices(i,j,k) = mesh_info.cell_vertices(parent_cell,j,k);
  }

  // Now we have to set the connectivity for the faces.
  sideset_info.has_connectivity = mesh_info.has_connectivity;
  sideset_info.face_to_cells = Kokkos::View<panzer::LocalOrdinal*[2]>("face_to_cells",num_faces);
  sideset_info.face_to_lidx = Kokkos::View<panzer::LocalOrdinal*[2]>("face_to_lidx",num_faces) ;
  sideset_info.cell_to_faces = Kokkos::View<panzer::LocalOrdinal**>("cell_to_faces",num_total_cells,num_faces_per_cell);

  Kokkos::deep_copy(sideset_info.face_to_cells,-1);
  Kokkos::deep_copy(sideset_info.face_to_lidx,-1);
  Kokkos::deep_copy(sideset_info.cell_to_faces,-1);

  // Default the system with invalid cell index - this will be most of the entries
  Kokkos::deep_copy(sideset_info.cell_to_faces, -1);

  for(int face_index=0;face_index<num_faces;++face_index){
    const auto & face = faces[face_index];
    const auto & cell_0 = face.cell_0;
    const auto & cell_1 = face.cell_1;

    // TODO: Super expensive... need to find a better way
    const auto sideset_cell_0 = std::distance(all_parent_cells.begin(), std::find(all_parent_cells.begin(), all_parent_cells.end(), cell_0));
    const auto sideset_cell_1 = std::distance(all_parent_cells.begin(), std::find(all_parent_cells.begin()+num_owned_cells, all_parent_cells.end(), cell_1));

    sideset_info.face_to_cells(face_index,0) = sideset_cell_0;
    sideset_info.face_to_cells(face_index,1) = sideset_cell_1;

    sideset_info.face_to_lidx(face_index,0) = face.subcell_index_0;
    sideset_info.face_to_lidx(face_index,1) = face.subcell_index_1;

    sideset_info.cell_to_faces(sideset_cell_0,face.subcell_index_0) = face_index;
    sideset_info.cell_to_faces(sideset_cell_1,face.subcell_index_1) = face_index;

  }

}

} // namespace

Teuchos::RCP<panzer::LocalMeshInfo>
generateLocalMeshInfo(const panzer_stk::STK_Interface & mesh)
{

  TEUCHOS_FUNC_TIME_MONITOR_DIFF("panzer_stk::generateLocalMeshInfo",GenerateLocalMeshInfo);

  using Teuchos::RCP;
  using Teuchos::rcp;

  auto mesh_info_rcp = Teuchos::rcp(new panzer::LocalMeshInfo());
  auto & mesh_info = *mesh_info_rcp;

  // 'include_connectivity' is a bit of a tricky thing right now.
  // eventually we will need to only build the internal connectivity if we require it,
  // but most of this algorithm requires connectivity to work in the first place
  // so we build it anyway, and only pass it down to the internal partitions if needed

  //typedef Tpetra::CrsMatrix<int,LO,GO> crs_type;
  typedef Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> map_type;
  typedef Tpetra::Import<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> import_type;

  // Make sure the STK interface is valid
  TEUCHOS_ASSERT(mesh.isInitialized());

  // This is required by some of the STK stuff
  TEUCHOS_ASSERT(typeid(panzer::LocalOrdinal) == typeid(int));

  // We're allowed to do this since the connection manager only exists in this scope... even though it is also an RCP...
  RCP<panzer::ConnManager > conn_rcp = rcp(new panzer_stk::STKConnManager(Teuchos::rcpFromRef(mesh)));
  panzer::ConnManager & conn = *conn_rcp;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = mesh.getComm();

  // To create the local mesh info, we need the geometry (cell vertices) and the topology (ConnManager)
  PHX::View<const double***> owned_cell_vertices;
  {
    // Lets fill some arrays - note local_cell_ids doesn't contain virtual cells yet
    std::vector<size_t> owned_cell_local_ids;
    {    // Get list of all elements that are owned
      std::vector<stk::mesh::Entity> owned_elements;
      stk::mesh::Selector owned_selector = mesh.getMetaData()->locally_owned_part();
      stk::mesh::get_selected_entities(owned_selector,mesh.getBulkData()->buckets(mesh.getElementRank()),owned_elements);

      // Add local cell ids for owned elements
      for(const auto & element : owned_elements)
        owned_cell_local_ids.push_back(mesh.elementLocalId(element));
    }

    std::vector<shards::CellTopology> element_block_topologies;
    conn.getElementBlockTopologies(element_block_topologies);

    TEUCHOS_ASSERT(element_block_topologies.size() > 0);

    const shards::CellTopology & topology = element_block_topologies[0];

    // FIXME: We assume that all element blocks have the same topology.
    for(const auto & other_topology : element_block_topologies){
      TEUCHOS_ASSERT(other_topology.getKey() == topology.getKey());
    }

    const size_t num_vertices_per_cell = topology.getNodeCount();
    const size_t num_dims_per_vertex = topology.getDimension();

    // Get the vertices for the owned cells on this process
    auto owned_cell_vertices_temp = PHX::View<double***>("owned_cell_vertices",owned_cell_local_ids.size(),num_vertices_per_cell,num_dims_per_vertex);
    mesh.getElementVerticesNoResize(owned_cell_local_ids,owned_cell_vertices_temp);
    owned_cell_vertices = owned_cell_vertices_temp;
  }


  ///////////////////////////////////////////////////////////
  // STK is only required for 3 things:
  // 1) Comm
  // 2) ConnManager
  // 3) Owned cell vertices

  // This builds most of the mesh info
  mesh_info.initialize(comm, conn, owned_cell_vertices);

  // STK is still required for some of it

  ////////////////////////////////////////////////////////////////////////////////////
  // At this point all the data structures have been built.
  // Let's store it in the mesh_info object
  {

    // Setup element blocks and sidesets
    std::vector<std::string> element_blocks;
    mesh.getElementBlockNames(element_blocks);

    // Setup element blocks and sidesets
    std::vector<std::string> sidesets;
    mesh.getSidesetNames(sidesets);

    // Setup sidesets
    for(const auto & element_block : element_blocks){
      for(const auto & sideset : sidesets){
        PANZER_FUNC_TIME_MONITOR_DIFF("Setup LocalMeshSidesetInfo",SetupLocalMeshSidesetInfo);
        setupLocalMeshSidesetInfo(mesh, conn, mesh_info, element_block, sideset, mesh_info.sidesets);
      }
    }
  }

  return mesh_info_rcp;

}

}
