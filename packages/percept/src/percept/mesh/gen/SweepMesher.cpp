// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stdexcept>
#include <iostream>
#include <utility>
#include <set>

#undef NDEBUG

#include <cassert>

#include "SweepMesher.hpp"

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>
#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ionit_Initializer.h>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>

namespace percept
  {
    using namespace util;

    enum { SpatialDim = 3 };
    enum { node_count = 21 };
    enum { number_hex = 3 };
    enum { number_wedge = 3 };
    enum { number_tetra = 3 };
    enum { number_pyramid = 2 };
    enum { number_shell_quad = 3 };
    enum { number_shell_tri = 3 };

    template<> void SweepMesher::breakElement<shards_Quadrilateral_4, shards_Triangle_3>(unsigned elemIndex)
    {
      unsigned* elem = &m_elems[shards_Quadrilateral_4][elemIndex*m_elemInfo[shards_Quadrilateral_4].vertex_count];
      VectorOfInt newElem = VectorOfInt(3);
      newElem[0] = elem[0];
      newElem[1] = elem[1];
      newElem[2] = elem[2];
      push_back(m_elems[shards_Triangle_3], newElem);
      newElem[0] = elem[0];
      newElem[1] = elem[2];
      newElem[2] = elem[3];
      push_back(m_elems[shards_Triangle_3], newElem);
      if (m_deleteAfterBreak)
        {
          VectorOfInt::iterator pos = m_elems[shards_Quadrilateral_4].begin()+(elemIndex*m_elemInfo[shards_Quadrilateral_4].vertex_count);
          m_elems[shards_Quadrilateral_4].erase(pos, pos + m_elemInfo[shards_Quadrilateral_4].vertex_count);
        }
    }

    /// based on UseCase_3 in stk_mesh/use_cases - creates nodes and elements in stk::mesh database
    void SweepMesher::stkMeshCreate(stk::ParallelMachine& comm)
    {
      stkMeshCreateMetaNoCommit(comm);
      m_metaData->commit();
      stkMeshCreateBulkAfterMetaCommit(comm);
    }

    static std::vector<std::string> get_entity_rank_names(unsigned dim)
    {
      std::vector<std::string> names = stk::mesh::entity_rank_names();
#if PERCEPT_USE_FAMILY_TREE
      names.push_back("FAMILY_TREE");
#endif
      return names;
    }

    void SweepMesher::stkMeshCreateMetaNoCommit(stk::ParallelMachine& comm)
    {
      stk::mesh::MeshBuilder builder(comm);
      builder.set_spatial_dimension(3);
      builder.set_entity_rank_names(get_entity_rank_names(3u));
      m_bulkData = builder.create();
      m_metaData = std::shared_ptr<stk::mesh::MetaData>(&m_bulkData->mesh_meta_data(), [](auto ptrWeWontDelete){});

      m_parts.resize(NUM_ELEM_TYPES);

      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          if (m_elems[ieletype].size() > 0)
            {
              m_parts[ieletype] = &(m_metaData->declare_part( std::string("block_").append(std::string(m_elemInfo[ieletype].name)) , stk::topology::ELEMENT_RANK ));
            }
        }
      m_coordinates_field = &m_metaData->declare_field<double>( stk::topology::NODE_RANK, "coordinates" );

      // set cell topology
      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          if (m_elems[ieletype].size() > 0)
            {
              m_elemInfo[ieletype].setCellTopoFptr(*m_parts[ieletype]);
              stk::io::put_io_part_attribute(*m_parts[ieletype]);
            }
        }

      // Field restrictions:
      stk::mesh::Part & universal = m_metaData->universal_part();

      put_field_on_mesh( *m_coordinates_field , universal , m_metaData->spatial_dimension(), nullptr);

      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          if (m_elems[ieletype].size() > 0)
            {
#if 0
              std::cout << "shards::Hexahedron<> ::node_count = " << shards::Hexahedron<> ::node_count
                        << " " << m_elemInfo[ieletype].node_count << std::endl;
              std::cout << "shards::Wedge<> ::node_count = " << shards::Wedge<> ::node_count
                        << " " << m_elemInfo[ieletype].node_count << std::endl;
#endif
            }
        }
    }

    void SweepMesher::stkMeshCreateBulkAfterMetaCommit(stk::ParallelMachine& comm)
    {

      //stk::mesh::BulkData & bulkData = modifiableBulkData();
      stk::mesh::BulkData & bulkData = *m_bulkData;
      bulkData.modification_begin();

      stk::mesh::EntityId elem_id = 1 ;

      for (unsigned ieletype = 0; ieletype < NUM_ELEM_TYPES; ieletype++)
        {
          if (m_elems[ieletype].size() > 0)
            {
              static stk::mesh::EntityIdVector node_ids(27);  // FIXME - do we have more than 27 nodes?
              stk::mesh::Part & part =  *m_parts[ieletype];
              unsigned nodes_per_elem = m_elemInfo[ieletype].vertex_count;
              unsigned numElems = m_elems[ieletype].size()/nodes_per_elem;

              //std::cout << " elems[" << m_elemInfo[ieletype].name << "] = " << numElems << std::endl;

              for (unsigned jelem = 0; jelem < numElems; jelem++, ++elem_id)
                {
                  for (unsigned inode = 0; inode < m_elemInfo[ieletype].vertex_count; inode++)
                    {
                      node_ids[inode] = 1+m_elems[ieletype][jelem*m_elemInfo[ieletype].vertex_count + inode];
                    }

                  //std::cout << "elem_id = " << elem_id << std::endl;
                  stk::mesh::declare_element( bulkData , part , elem_id , node_ids );
                }
            }
        }

      unsigned node_count_1 = m_node_coords.size();
      for ( unsigned i = 0 ; i < node_count_1 ; ++i ) {
        stk::mesh::Entity const node = m_bulkData->get_entity( stk::topology::NODE_RANK , i + 1 );
        double * const coord = stk::mesh::field_data( *m_coordinates_field , node );

        coord[0] = m_node_coords[i][0];
        coord[1] = m_node_coords[i][1];
        coord[2] = m_node_coords[i][2];
      }

      stk::mesh::fixup_ghosted_to_shared_nodes(bulkData);
      bulkData.modification_end();
    }

    void SweepMesher::writeSTKMesh(const char *filename)
    {
      const std::string out_filename(filename);

      stk::io::StkMeshIoBroker mesh;
      mesh.set_bulk_data(m_bulkData);
      size_t result_file_index = mesh.create_output_mesh(out_filename, stk::io::WRITE_RESULTS);
      mesh.write_output_mesh(result_file_index);
    }

  }//namespace percept
