// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


/**
 */

//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>


#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <percept/fixtures/TriQuadSurfaceMesh3D.hpp>

//----------------------------------------------------------------------

  namespace percept {

    typedef stk::topology::topology_type<stk::topology::HEX_8>          Hex8;
    typedef stk::topology::topology_type<stk::topology::WEDGE_6>        Wedge6;
    typedef stk::topology::topology_type<stk::topology::TET_4>          Tet4;
    typedef stk::topology::topology_type<stk::topology::PYRAMID_5>      Pyramid5;

    typedef stk::topology::topology_type<stk::topology::SHELL_QUAD_4>   ShellQuad4;
    typedef stk::topology::topology_type<stk::topology::SHELL_TRI_3>    ShellTriangle3;

    typedef stk::topology::topology_type<stk::topology::QUAD_4>         Quad4;
    typedef stk::topology::topology_type<stk::topology::TRI_3>          Triangle3;

    TriQuadSurfaceMesh3D::TriQuadSurfaceMesh3D( stk::ParallelMachine comm, bool doCommit )
      : m_bulkDataPtr(stk::mesh::MeshBuilder(comm).set_spatial_dimension(3)
                                                  .set_entity_rank_names(stk::mesh::entity_rank_names())
                                                  .create())
      , m_bulkData(*m_bulkDataPtr)
      , m_metaData(m_bulkData.mesh_meta_data())
      , m_sideset(0)

      , m_block_quad_shell( m_metaData.declare_part_with_topology( "block_1", stk::topology::SHELL_QUAD_4 ))
      , m_block_tri_shell(  m_metaData.declare_part_with_topology( "block_2", stk::topology::SHELL_TRI_3 ))


    {
      m_coordinates_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "coordinates" );
//      m_centroid_field    = &m_metaData.declare_field<double>( stk::topology::ELEMENT_RANK, "centroid" );
//      m_temperature_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "temperature" );
//      m_volume_field      = &m_metaData.declare_field<double>( stk::topology::ELEMENT_RANK, "volume" );

      // Define where fields exist on the mesh:
      stk::mesh::Part & universal = m_metaData.universal_part();

      put_field_on_mesh( *m_coordinates_field , universal , m_metaData.spatial_dimension(), nullptr);
      // put_field_on_mesh( *m_centroid_field , universal , m_metaData.spatial_dimension(), nullptr);
      // put_field_on_mesh( *m_temperature_field, universal , nullptr);
      // put_field_on_mesh( *m_volume_field, m_block_hex , nullptr);
      // put_field_on_mesh( *m_volume_field, m_block_wedge , nullptr);
      // put_field_on_mesh( *m_volume_field, m_block_tet , nullptr);
      // put_field_on_mesh( *m_volume_field, m_block_pyramid , nullptr);

      stk::io::put_io_part_attribute(  m_block_tri_shell);
      stk::io::put_io_part_attribute(  m_block_quad_shell);

      if (doCommit)
        m_metaData.commit();
    }

    TriQuadSurfaceMesh3D::~TriQuadSurfaceMesh3D()
    { }

    //------------------------------------------------------------------------------

    void TriQuadSurfaceMesh3D::populate(unsigned nNodes, unsigned nTriFaces, unsigned nQuadFaces, Point *coords, TriIds *triFaces, QuadIds *quadFaces,
                                        stk::mesh::EntityId *nodeIdMap)
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();
      int p_size = m_bulkData.parallel_size();
      int p_rank = m_bulkData.parallel_rank();
      //stk::mesh::EntityId offset = (1ull + p_rank)*(1ull << 50);
      if (1)
        {
          stk::mesh::EntityId curr_elem_id = 1 + p_rank;

          // For each element topology declare elements

          //stk::mesh::Entity tris[nTriFaces];

          for ( unsigned i = 0 ; i < nTriFaces ; ++i , curr_elem_id += p_size ) {
            TriIds & tids = triFaces[i];
            // for (unsigned j=0; j < 3; ++j)
            //   tids[j] += offset;
            stk::mesh::declare_element( m_bulkData, m_block_tri_shell, curr_elem_id,  tids );
          }
          for ( unsigned i = 0 ; i < nQuadFaces ; ++i , curr_elem_id += p_size ) {
            QuadIds & qids = quadFaces[i];
            // for (unsigned j=0; j < 4; ++j)
            //   qids[j] += offset;
            stk::mesh::declare_element( m_bulkData, m_block_quad_shell, curr_elem_id, qids );
            //std::cout << "qids= " << qids[0] << std::endl;
          }

          // For all nodes assign nodal coordinates
          for ( unsigned i = 0 ; i < nNodes ; ++i ) {
            stk::mesh::EntityId id = (nodeIdMap ? nodeIdMap[i] : i + 1);
            stk::mesh::Entity const node = m_bulkData.get_entity( stk::topology::NODE_RANK , id );
            //std::cout << "id= " << id << std::endl;
            if (!m_bulkData.is_valid(node))
              continue;
            double * const coord = stk::mesh::field_data( *m_coordinates_field , node );
            coord[0] = coords[i][0];
            coord[1] = coords[i][1];
            coord[2] = coords[i][2];
            //std::cout << "coord= " << coord[0] << std::endl;
         }
        }
      stk::mesh::fixup_ghosted_to_shared_nodes(m_bulkData);
      m_bulkData.modification_end();
    }


  } //namespace percept
