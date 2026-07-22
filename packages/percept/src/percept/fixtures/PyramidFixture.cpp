// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <percept/fixtures/PyramidFixture.hpp>
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
#include <percept/FieldTypes.hpp>
#include <percept/PerceptUtils.hpp>

//----------------------------------------------------------------------

  namespace percept {

    typedef stk::topology::topology_type<stk::topology::HEX_8>         Hex8;
    typedef stk::topology::topology_type<stk::topology::WEDGE_6>       Wedge6;
    typedef stk::topology::topology_type<stk::topology::TET_4>         Tet4;
    typedef stk::topology::topology_type<stk::topology::PYRAMID_5>     Pyramid5;

    typedef stk::topology::topology_type<stk::topology::SHELL_QUAD_4>  ShellQuad4;
    typedef stk::topology::topology_type<stk::topology::SHELL_TRI_3>   ShellTriangle3;

    typedef stk::topology::topology_type<stk::topology::QUAD_4>        Quad4;
    typedef stk::topology::topology_type<stk::topology::TRI_3>        Triangle3;

    PyramidFixture::PyramidFixture( stk::ParallelMachine comm, bool doCommit, bool do_sidesets ) :
      m_spatial_dimension(3)
      , m_bulkDataPtr(stk::mesh::MeshBuilder(comm).set_spatial_dimension(m_spatial_dimension)
                                                  .set_entity_rank_names(entity_rank_names_and_family_tree())
                                                  .create())
      , m_bulkData(*m_bulkDataPtr)
      , m_metaData(m_bulkData.mesh_meta_data())
      , m_block_pyramid(    m_metaData.declare_part_with_topology( "block_4", stk::topology::PYRAMID_5 ))
      , m_sideset_quad(0), m_sideset_quad_subset(0)
      , m_sideset_tri(0), m_sideset_tri_subset(0)
    {      
      m_coordinates_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "coordinates" );
      m_centroid_field    = &m_metaData.declare_field<double>( stk::topology::ELEMENT_RANK, "centroid" );
      m_temperature_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "temperature" );
      m_volume_field      = &m_metaData.declare_field<double>( stk::topology::ELEMENT_RANK, "volume" );

      // Define where fields exist on the mesh:
      stk::mesh::Part & universal = m_metaData.universal_part();

      if (do_sidesets)
        {
          m_sideset_quad_subset = &m_metaData.declare_part(std::string("surface_pyramid5_quad2d2_1"), stk::topology::FACE_RANK);
          m_sideset_quad =        &m_metaData.declare_part(std::string("surface_1"), stk::topology::FACE_RANK);
          stk::mesh::set_topology(*m_sideset_quad_subset, stk::topology::QUAD_4);
          stk::io::put_io_part_attribute(*m_sideset_quad_subset);
          stk::io::put_io_part_attribute(*m_sideset_quad);
          m_metaData.declare_part_subset(*m_sideset_quad, *m_sideset_quad_subset);

          m_sideset_tri_subset = &m_metaData.declare_part(std::string("surface_pyramid5_tri2d2_1"), stk::topology::FACE_RANK);
          m_sideset_tri =        &m_metaData.declare_part(std::string("surface_2"), stk::topology::FACE_RANK);
          stk::mesh::set_topology(*m_sideset_tri_subset, stk::topology::TRI_3);
          stk::io::put_io_part_attribute(*m_sideset_tri_subset);
          stk::io::put_io_part_attribute(*m_sideset_tri);
          m_metaData.declare_part_subset(*m_sideset_tri, *m_sideset_tri_subset);
        }

      put_field_on_mesh( *m_coordinates_field , universal, m_metaData.spatial_dimension(), nullptr);
      put_field_on_mesh( *m_centroid_field , universal, m_metaData.spatial_dimension(), nullptr);
      put_field_on_mesh( *m_temperature_field, universal, nullptr);

      put_field_on_mesh( *m_volume_field, m_block_pyramid, nullptr);

      stk::io::put_io_part_attribute(  m_block_pyramid );

      if (doCommit)
        m_metaData.commit();
      init();
    }

    PyramidFixture::~PyramidFixture()
    { }

    //------------------------------------------------------------------------------
    // Use case specific mesh generation data:
    //------------------------------------------------------------------------------

    void PyramidFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          stk::mesh::EntityId curr_elem_id = 1;

          // For each element topology declare elements

          stk::mesh::Entity pyramids[2];
          for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++curr_elem_id ) {
            pyramids[i] = stk::mesh::declare_element( m_bulkData, m_block_pyramid, curr_elem_id, pyramid_node_ids[i] );
          }

          if (m_sideset_quad)
            {
              for ( unsigned i = 0 ; i < number_quad ; ++i , ++curr_elem_id ) {
                m_bulkData.declare_element_side( pyramids[i], // element,
                                                 4,            //j_side, // local_side_ord,
                                                 stk::mesh::PartVector{m_sideset_quad_subset});
              }
            }

          if (m_sideset_tri)
            {
              unsigned j_side=0;
              for ( unsigned i = 0 ; i < 3 ; ++i , ++curr_elem_id ) {
                if (i == 2) ++j_side;
                m_bulkData.declare_element_side( pyramids[0], // element,
                                                 j_side,            //j_side, // local_side_ord,
                                                 stk::mesh::PartVector{m_sideset_tri_subset});
                ++j_side;
              }
              j_side=1;
              for ( unsigned i = 0 ; i < 3 ; ++i , ++curr_elem_id ) {
                  m_bulkData.declare_element_side( pyramids[1], // element,
                                                 j_side,            //j_side, // local_side_ord,
                                                 stk::mesh::PartVector{m_sideset_tri_subset});
                ++j_side;
              }
            }

          // For all nodes assign nodal coordinates
          for ( unsigned i = 0 ; i < node_count ; ++i ) {
            stk::mesh::Entity const node = m_bulkData.get_entity( stk::topology::NODE_RANK , i + 1 );
            double * const coord = stk::mesh::field_data( *m_coordinates_field , node );
            coord[0] = node_coord_data[i][0] ;
            coord[1] = node_coord_data[i][1] ;
            coord[2] = node_coord_data[i][2] ;
          }

        }
      stk::mesh::fixup_ghosted_to_shared_nodes(m_bulkData);
      m_bulkData.modification_end();

    }

    // Verify mesh for 6 different parts
    bool verifyMesh( const PyramidFixture & mesh )
    {
      bool result = true;

      const stk::mesh::BulkData & bulkData = mesh.m_bulkData ;
      //const CoordinatesFieldType & node_coord = mesh.m_coordinates_field ;

      stk::mesh::BucketVector element_buckets = bulkData.buckets( stk::topology::ELEMENT_RANK );

      // Create a pair containing Part and matching node_count

      typedef std::pair<stk::mesh::Part*, unsigned> PartNodeCountPair;
      std::vector<PartNodeCountPair> part_and_node_counts;
      const unsigned pyramid5_num_nodes = Pyramid5::num_nodes;
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_pyramid, pyramid5_num_nodes));

      // Verify that entities in each part are set up correctly.
      // Use a PartVector iterator for parts_to_check and call
      // verify_elem_node_coord_by_part in UseCase_Common.cpp for
      // each part in turn.
#if 0
      for( std::vector<PartNodeCountPair>::const_iterator  i = part_and_node_counts.begin() ; i != part_and_node_counts.end() ; ++i )
        {
          result = result &&
            verify_elem_node_coord_by_part(
                                           *(i->first),
                                           element_buckets,
                                           elem_node_coord,
                                           node_coord,
                                           i->second
                                           );
        }
#endif

      // Check that all the nodes were allocated.
      for ( unsigned i = 0 ; i < mesh.node_count ; ++i ) {
        stk::mesh::Entity const node = bulkData.get_entity( stk::topology::NODE_RANK , i + 1 );
        if ( !bulkData.is_valid(node) ) {
          std::cerr << "Error!  Invalid null pointer for node returned from "
                    << "bulkData.get_entity( stk::topology::NODE_RANK, " << i+1 << " ) " << std::endl;
          result = false;
        }
      }

      return result;
    }

  } //namespace percept
