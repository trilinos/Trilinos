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

template <typename LO, typename GO>
std::vector<panzer::LocalMeshChunk<LO,GO> >
generateLocalMeshChunks(const panzer_stk::STK_Interface & mesh,
                        const panzer::WorksetDescriptor & description)
{

  std::vector<panzer::LocalMeshChunk<LO,GO> > chunks;

  // Get element block, local chunk of mesh, and bulk data
//  const stk::mesh::Part & element_block = *(mesh.getElementBlockPart(description.getElementBlock()));
//  const stk::mesh::Part & local_mesh = *(mesh.getOwnedPart());
//  const stk::mesh::BulkData & bulk_data = *(mesh.getBulkData());

  if(description.useSideset()){

    // This isn't even close to working...
    TEUCHOS_ASSERT(false);


//    const GO face_subcell_dimension = mesh.getCellTopology(description.getElementBlock())->getDimension()-1;
//
//    const int requested_chunk_size = description.getWorksetSize();
//
//
//    Kokkos::DynRankView<double,PHX::Device> data_allocation;
//


//
//    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//    const bool is_domain_boundary = true;
//
//    //===================================================
//    // Sideset
//
//    // This is a sideset, therefore we must... do something...
//    const std::string & sideset_name = description.getSideset();
//
//    std::vector<stk::mesh::Entity> side_entities;
//
//    // Grab the side entities associated with this sideset on the mesh
//    // Note: Throws exception if element block or sideset doesn't exist
//    try{
//
//      mesh.getMySides(sideset_name,description.getElementBlock(),side_entities);
//
//    } catch(STK_Interface::SidesetException & e) {
//       std::stringstream ss;
//       std::vector<std::string> sideset_names;
//       mesh.getSidesetNames(sideset_names);
//
//       // build an error message
//       ss << e.what() << "\nChoose one of:\n";
//       for(const auto & optional_sideset_name : sideset_names){
//          ss << "\"" << optional_sideset_name << "\"\n";
//       }
//
//       TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
//
//    } catch(STK_Interface::ElementBlockException & e) {
//       std::stringstream ss;
//       std::vector<std::string> element_block_names;
//       mesh.getElementBlockNames(element_block_names);
//
//       // build an error message
//       ss << e.what() << "\nChoose one of:\n";
//       for(const auto & optional_element_block_name : element_block_names){
//          ss << "\"" << optional_element_block_name << "\"\n";
//       }
//
//       TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
//
//    } catch(std::logic_error & e) {
//       std::stringstream ss;
//       ss << e.what() << "\nUnrecognized logic error.\n";
//
//       TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(true,std::logic_error,ss.str());
//
//    }
//
//    // We now have a list of sideset entities, lets unwrap them and create some chunks!
//
//    // Lets go all out and grab everything associated with that side
//
//    std::map<std::pair<int,int>, std::vector<GO> > local_cell_indexes_by_subcell;
//
//    // List of all elements - not really needed
//    std::vector<stk::mesh::Entity> elements;
//
//    // List of local subcell indexes indexed by element:
//    // For example: a Tet (element) would have
//    //  - 4 triangular faces (subcell_index 0-3, subcell_dimension=2)
//    //  - 6 edges (subcell_index 0-5, subcell_dimension=1)
//    //  - 4 vertices (subcell_index 0-3, subcell_dimension=0)
//    // Another example: a Line (element) would have
//    //  - 2 vertices (subcell_index 0-1, subcell_dimension=0)
//    std::vector<GO> subcell_indexes;
//    std::vector<GO> subcell_dimensions;
//
//    // Grab everything
//    panzer_stk::workset_utils::getSideElementCascade(mesh, description.getElementBlock(),side_entities, subcell_dimensions,subcell_indexes,elements);
//
//    // build local cell_ids, mapped by local side id
//    for(GO element_index=0; element_index<elements.size(); ++element_index) {
//      const GO subcell_dimension = subcell_dimensions[element_index];
//      const GO subcell_index = subcell_indexes[element_index];
//      const GO element_local_index = mesh.elementLocalId(elements[element_index]);
//
//      // Add subcell to map
//      local_cell_indexes_by_subcell[std::pair<int,int>(subcell_dimension, subcell_index)].push_back(element_local_index);
//    }
//
//    // We now have a list of all local cells 'touching' the side set, as well as the faces, edges, and vertices that define the 'touch'
//
//    // For our purposes we only really want the cells and faces associated with the side
//
//    std::vector<GO> local_cell_indexes;
//    {
//      // We use a set to avoid duplicates
//      std::set<GO> local_cell_indexes_set;
//      for(const auto & cell_subcell_pair : local_cell_indexes_by_subcell){
//        const std::pair<GO,GO> & subcell_definition = cell_subcell_pair.first;
//        if(subcell_definition.first == face_subcell_dimension){
//          const std::vector<GO> & local_cell_indexes_for_subcell = cell_subcell_pair.second;
//          for(const GO & local_cell_index : local_cell_indexes_for_subcell){
//            local_cell_indexes_set.insert(local_cell_index);
//          }
//        }
//      }
//
//      for(const GO & cell_index : local_cell_indexes_set){
//        local_cell_indexes.push_back(cell_index);
//      }
//    }
//
//
//    const GO num_total_cells = local_cell_indexes.size();
//    const GO num_left_over = num_total_cells % requested_chunk_size;
//    const GO num_chunks = num_total_cells / requested_chunk_size + ((num_left_over > 0)?1:0);
//
//    // We create our chunks
//    for(GO chunk_index=0; chunk_index<num_chunks; ++chunk_index){
//      const GO cell_index_start = chunk_index * requested_chunk_size;
//      const GO cell_index_end = std::min(num_total_cells, cell_index_start + requested_chunk_size);
//      const GO num_real_cells = cell_index_end - cell_index_start;
//      const GO num_virtual_cells = (is_domain_boundary)?num_real_cells:0;
//      const GO num_cells = num_real_cells + num_virtual_cells;
//
//      chunks.push_back(panzer::LocalMeshChunk<LO,GO>());
//      panzer::LocalMeshChunk<LO,GO> & chunk = chunks.back();
//
//      chunk.num_cells = num_cells;
//      chunk.subcell_dimension = face_subcell_dimension;
//      chunk.element_block_name = description.getElementBlock();
//      chunk.sideset_name = sideset_name;
//      chunk.local_cell_global_indexes.assign(local_cell_indexes.begin()+cell_index_start,local_cell_indexes.begin()+cell_index_end);
//
//      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//      // Side has either virtual or ghost cells
//      chunk.virtual_cell_global_indexes;
//
//      // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//      // Side has either virtual or ghost cells
//      chunk.ghost_cell_global_indexes;
//
//      // Must take virtual cells into account
//      {
//        mesh.getElementVertices(chunk.local_cell_global_indexes, description.getElementBlock(), data_allocation);
//
//        const GO num_cells_with_vertices = data_allocation.dimension_0();
//        const GO num_vertices_per_cell = data_allocation.dimension_1();
//        const GO num_dims = data_allocation.dimension_2();
//
//        TEUCHOS_ASSERT(num_real_cells == num_cells_with_vertices);
//
//        chunk.cell_vertices = Kokkos::View<double***>("cell_vertices", num_cells,num_vertices_per_cell,num_dims);
//        Kokkos::deep_copy(chunk.cell_vertices,0.);
//
//        for(GO i=0;i<num_cells_with_vertices;++i){
//          for(GO j=0;j<num_vertices_per_cell;++j){
//            for(GO k=0;k<num_dims;++k){
//              chunk.cell_vertices(i,j,k) = data_allocation(i,j,k);
//            }
//          }
//        }
//      }
//
//      // Each face in this side set connects a real cell to a ghost or virtual cell
//
//      const std::vector<GO> & other_side_cells = (is_domain_boundary) ? chunk.virtual_cell_global_indexes : chunk.ghost_cell_global_indexes;
//
//      chunk.face_to_cells = Kokkos::View<LO*[2],PHX::Device>("face_to_cells", num_real_cells);
//      for(GO i=0; i<num_real_cells; ++i){
//
//      }
//
//
//
//
//
//    }


    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    //===================================================

  } else {

    // We need to partition this chunk of the mesh
    TEUCHOS_ASSERT(description.getWorksetSize()==panzer::WorksetDescriptor::FULL);

    panzer::LocalMeshInfo<LO,GO> mesh_info = panzer_stk::generateLocalMeshInfo<LO,GO>(mesh, description.getElementBlock());

    chunks.push_back(panzer::LocalMeshChunk<LO,GO>());
    panzer::LocalMeshChunk<LO,GO> & chunk = chunks.back();

    const GO num_owned_cells = mesh_info.owned_cells.dimension(0);
    const GO num_ghost_cells = mesh_info.ghstd_cells.dimension(0);
    const GO num_virtual_cells = mesh_info.virtual_cells.dimension(0);
    const GO num_total_cells = num_owned_cells+num_ghost_cells+num_virtual_cells;
    const GO num_vertices_per_cell = mesh_info.owned_vertices.dimension(1);
    const GO num_dims_per_vertex = mesh_info.owned_vertices.dimension(2);

    chunk.element_block_name = description.getElementBlock();

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
