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

#include <percept/fixtures/HeterogeneousFixture.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <percept/PerceptUtils.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>


#include <stk_mesh/base/FEMHelpers.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
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

    HeterogeneousFixture::HeterogeneousFixture( stk::ParallelMachine comm, bool doCommit, bool do_sidesets, bool do_one_sideset ) :
      m_spatial_dimension(3)
      , m_bulkDataPtr(stk::mesh::MeshBuilder(comm).set_spatial_dimension(m_spatial_dimension)
                                                  .set_entity_rank_names(entity_rank_names_and_family_tree())
                                                  .create())
      , m_bulkData(*m_bulkDataPtr)
      , m_metaData(m_bulkData.mesh_meta_data())
      , m_block_hex(        m_metaData.declare_part_with_topology(  "block_1", stk::topology::HEX_8 ))
      , m_block_wedge(      m_metaData.declare_part_with_topology( "block_2", stk::topology::WEDGE_6 ))
      , m_block_tet(        m_metaData.declare_part_with_topology( "block_3", stk::topology::TET_4 ))
      , m_block_pyramid(    m_metaData.declare_part_with_topology( "block_4", stk::topology::PYRAMID_5 ))

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      , m_block_quad_shell( m_metaData.declare_part_with_topology( "block_5", stk::topology::SHELL_QUAD_4 ))
      , m_block_tri_shell(  m_metaData.declare_part_with_topology( "block_6", stk::topology::SHELL_TRI_3 ))
#endif
      , m_sideset_quad(0), m_sideset_quad_subset(0)
      , m_sideset_tri(0), m_sideset_tri_subset(0)
    {
      m_coordinates_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "coordinates" );
      m_centroid_field    = &m_metaData.declare_field<double>( stk::topology::ELEMENT_RANK, "centroid" );
      m_temperature_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "temperature" );
      m_volume_field      = &m_metaData.declare_field<double>( stk::topology::ELEMENT_RANK, "volume" );

      // Define where fields exist on the mesh:
      stk::mesh::Part & universal = m_metaData.universal_part();

      put_field_on_mesh( *m_coordinates_field , universal, m_metaData.spatial_dimension(), nullptr);
      put_field_on_mesh( *m_centroid_field , universal, m_metaData.spatial_dimension(), nullptr);
      put_field_on_mesh( *m_temperature_field, universal, nullptr);
      put_field_on_mesh( *m_volume_field, m_block_hex, nullptr);
      put_field_on_mesh( *m_volume_field, m_block_wedge, nullptr);
      put_field_on_mesh( *m_volume_field, m_block_tet, nullptr);
      put_field_on_mesh( *m_volume_field, m_block_pyramid, nullptr);

      stk::io::put_io_part_attribute(  m_block_hex );
      stk::io::put_io_part_attribute(  m_block_wedge );
      stk::io::put_io_part_attribute(  m_block_tet );
      stk::io::put_io_part_attribute(  m_block_pyramid );

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      stk::io::put_io_part_attribute(  m_block_quad_shell );
      stk::io::put_io_part_attribute(  m_block_tri_shell );
#endif

      if (do_sidesets)
        {
          if (!do_one_sideset)
            {
              m_sideset_quad_subset = &m_metaData.declare_part(std::string("surface_wedge5_quad2d2_1"), stk::topology::FACE_RANK);
              m_sideset_quad =        &m_metaData.declare_part(std::string("surface_1"), stk::topology::FACE_RANK);
              stk::mesh::set_topology(*m_sideset_quad_subset, stk::topology::QUAD_4);
              stk::io::put_io_part_attribute(*m_sideset_quad_subset);
              stk::io::put_io_part_attribute(*m_sideset_quad);
              m_metaData.declare_part_subset(*m_sideset_quad, *m_sideset_quad_subset);

              m_sideset_tri_subset = &m_metaData.declare_part(std::string("surface_wedge5_tri2d2_1"), stk::topology::FACE_RANK);
              m_sideset_tri =        &m_metaData.declare_part(std::string("surface_2"), stk::topology::FACE_RANK);
              stk::mesh::set_topology(*m_sideset_tri_subset, stk::topology::TRI_3);
              stk::io::put_io_part_attribute(*m_sideset_tri_subset);
              stk::io::put_io_part_attribute(*m_sideset_tri);
              m_metaData.declare_part_subset(*m_sideset_tri, *m_sideset_tri_subset);
            }
          else
            {
              m_sideset_quad_subset = &m_metaData.declare_part(std::string("surface_wedge5_quad2d2_1"), stk::topology::FACE_RANK);
              m_sideset_quad =        &m_metaData.declare_part(std::string("surface_1"), stk::topology::FACE_RANK);
              stk::mesh::set_topology(*m_sideset_quad_subset, stk::topology::QUAD_4);
              stk::io::put_io_part_attribute(*m_sideset_quad_subset);
              stk::io::put_io_part_attribute(*m_sideset_quad);
              m_metaData.declare_part_subset(*m_sideset_quad, *m_sideset_quad_subset);

              m_sideset_tri_subset = &m_metaData.declare_part(std::string("surface_wedge5_tri2d2_1"), stk::topology::FACE_RANK);
              m_sideset_tri = m_sideset_quad; //       &m_metaData.declare_part(std::string("surface_1"), stk::topology::FACE_RANK);
              stk::mesh::set_topology(*m_sideset_tri_subset, stk::topology::TRI_3);
              stk::io::put_io_part_attribute(*m_sideset_tri_subset);
              //stk::io::put_io_part_attribute(*m_sideset_tri);
              m_metaData.declare_part_subset(*m_sideset_tri, *m_sideset_tri_subset);

            }
        }
      if (doCommit)
        m_metaData.commit();
      init();
    }

    HeterogeneousFixture::~HeterogeneousFixture()
    { }

    //------------------------------------------------------------------------------
    // Use case specific mesh generation data:

    enum { SpatialDim = 3 };
    enum { node_count = 21 };
    enum { number_hex = 3 };
    enum { number_wedge = 3 };
    enum { number_tetra = 3 };
    enum { number_pyramid = 2 };
    enum { number_shell_quad = 3 };
    enum { number_shell_tri = 3 };
    enum { number_quad = 3 };
    enum { number_tri = 2 };

    namespace {

#if 0
      // Hard coded node coordinate data for all the nodes in the entire mesh
      static const double node_coord_data[ node_count ][ SpatialDim ] = {
        { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,
        { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } , { 3 , 1 , 0 } ,
        { 0 , 2 , 0 } , { 1 , 2 , 0 } ,
        { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
        { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
        { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
        { 1 , 1 , -2 } };
#endif

      // Hard coded hex node ids for all the hex nodes in the entire mesh
      static const stk::mesh::EntityIdVector hex_node_ids[number_hex] {
        { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
        { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
        { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 } };

      // Hard coded wedge node ids for all the wedge nodes in the entire mesh
      static const stk::mesh::EntityIdVector wedge_node_ids[number_wedge] {
        { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
        { 20 , 19 , 16 , 10 ,  9 ,  6 } ,
        { 16 , 17 , 20 ,  6 ,  7 , 10 } };

      // Hard coded tetra node ids for all the tetra nodes in the entire mesh
      static const stk::mesh::EntityIdVector tetra_node_ids[number_tetra] {
        { 15 , 19 , 16 , 21 } ,
        { 19 , 20 , 16 , 21 } ,
        { 16 , 20 , 17 , 21 } };

      // Hard coded pyramid node ids for all the pyramid nodes in the entire mesh
      static const stk::mesh::EntityIdVector pyramid_node_ids[number_pyramid] {
        { 11 , 15 , 16 , 12 , 21 } ,
        { 12 , 16 , 17 , 13 , 21 } };

#if 0
      // Hard coded shell quad node ids for all the shell quad nodes in the entire mesh
      static const stk::mesh::EntityIdVector shell_quad_node_ids[number_shell_quad] {
        { 9 , 6 , 16 , 19 } ,
        { 6 , 7 , 17 , 16 } ,
        { 7 , 8 , 18 , 17 } };

      // Hard coded shell tri node ids for all the shell tri nodes in the entire mesh
      static const stk::mesh::EntityIdVector shell_tri_node_ids[number_shell_tri] {
        { 19 , 16 , 21 } ,
        { 16 , 17 , 21 } ,
        { 17 , 13 , 21 } };

      // NOTE: some quad, tri's for wedge sideset testing
      // Hard coded quad node ids for all the quad nodes in the entire mesh
      static const stk::mesh::EntityIdVector quad_node_ids[number_quad] {
        { 5, 9, 19, 15},
        { 7, 17, 20, 10 },
        { 10, 20, 19, 9}
      };
#endif

      // wedge element id, side id
      static const stk::mesh::EntityIdVector quad_node_side_ids[number_quad] {
        {4, 2},
        {6, 1},
        {5, 0}
      };

#if 0
      // Hard coded tri node ids for all the tri nodes in the entire mesh
      static const stk::mesh::EntityIdVector tri_node_ids[number_tri] {
        { 5, 6, 9},
        { 6, 10, 9}
      };
#endif

      // wedge element id, side id
      static const stk::mesh::EntityIdVector tri_node_side_ids[number_quad] {
        {4, 4},
        {5, 3}
      };

    }

    //------------------------------------------------------------------------------

    void HeterogeneousFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          stk::mesh::EntityId curr_elem_id = 1;

          // For each element topology declare elements

          std::vector<stk::mesh::Entity> wedges(number_wedge);

          for ( unsigned i = 0 ; i < number_hex ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_hex, curr_elem_id, hex_node_ids[i] );
          }

          for ( unsigned i = 0 ; i < number_wedge ; ++i , ++curr_elem_id ) {
            wedges[i] = stk::mesh::declare_element( m_bulkData, m_block_wedge, curr_elem_id, wedge_node_ids[i] );
          }

          for ( unsigned i = 0 ; i < number_tetra ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_tet, curr_elem_id, tetra_node_ids[i] );
          }

          for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_pyramid, curr_elem_id, pyramid_node_ids[i] );
          }

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
          for ( unsigned i = 0 ; i < number_shell_quad ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_quad_shell, curr_elem_id, shell_quad_node_ids[i]);
          }

          for ( unsigned i = 0 ; i < number_shell_tri ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_tri_shell, curr_elem_id, shell_tri_node_ids[i] );
          }
#endif

          if (m_sideset_quad)
            {
              for ( unsigned i = 0 ; i < number_quad ; ++i , ++curr_elem_id ) {
                std::cout << "quad i= " << i << std::endl;
                m_bulkData.declare_element_side( wedges[quad_node_side_ids[i][0] - 4], // element,
                                                 quad_node_side_ids[i][1],   //j_side, // local_side_ord,
                                                 stk::mesh::PartVector{m_sideset_quad_subset});
              }
            }

          if (m_sideset_tri)
            {
              for ( unsigned i = 0 ; i < number_tri ; ++i , ++curr_elem_id ) {
                std::cout << "tri i= " << i << std::endl;
                m_bulkData.declare_element_side( wedges[tri_node_side_ids[i][0] - 4], // element,
                                                 tri_node_side_ids[i][1],   //j_side, // local_side_ord,
                                                 stk::mesh::PartVector{m_sideset_tri_subset});
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
    bool verifyMesh( const HeterogeneousFixture & mesh )
    {
      bool result = true;

      const stk::mesh::BulkData & bulkData = mesh.m_bulkData ;
      //const CoordinatesFieldType & node_coord = mesh.m_coordinates_field ;

      stk::mesh::BucketVector element_buckets = bulkData.buckets( stk::topology::ELEMENT_RANK );

      // Create a pair containing Part and matching node_count

      typedef std::pair<stk::mesh::Part*, unsigned> PartNodeCountPair;
      std::vector<PartNodeCountPair> part_and_node_counts;
      const unsigned hex8_num_nodes = Hex8::num_nodes;
      const unsigned wedge6_num_nodes = Wedge6::num_nodes;
      const unsigned tet4_num_nodes = Tet4::num_nodes;
      const unsigned pyramid5_num_nodes = Pyramid5::num_nodes;
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_hex, hex8_num_nodes));
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_wedge, wedge6_num_nodes));
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_tet, tet4_num_nodes));
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_pyramid, pyramid5_num_nodes));

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_quad_shell, ShellQuad4::num_nodes));
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_tri_shell, ShellTriangle3::num_nodes));
#endif

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
      for ( unsigned i = 0 ; i < node_count ; ++i ) {
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
