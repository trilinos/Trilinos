// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#include <percept/fixtures/TetWedgeFixture.hpp>
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
//----------------------------------------------------------------------

  namespace percept {

    typedef stk::topology::topology_type<stk::topology::WEDGE_6>       Wedge6;
    typedef stk::topology::topology_type<stk::topology::TET_4>         Tet4;

    typedef stk::topology::topology_type<stk::topology::QUAD_4>        Quad4;
    typedef stk::topology::topology_type<stk::topology::TRI_3>         Triangle3;

    TetWedgeFixture::TetWedgeFixture( stk::ParallelMachine comm, bool doCommit, bool do_sidesets ) :
      m_spatial_dimension(3)
      , m_bulkDataPtr(stk::mesh::MeshBuilder(comm).set_spatial_dimension(m_spatial_dimension).create())
      , m_bulkData(*m_bulkDataPtr)
      , m_metaData(m_bulkDataPtr->mesh_meta_data())
      , m_block_wedge(      m_metaData.declare_part_with_topology( "block_2", stk::topology::WEDGE_6 ))
      , m_block_tet(        m_metaData.declare_part_with_topology( "block_3", stk::topology::TET_4 ))

      , m_sideset_quad(0), m_sideset_quad_subset(0)
      , m_sideset_tri(0), m_sideset_tri_subset(0)

    {
      m_coordinates_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "coordinates" );
      m_centroid_field    = &m_metaData.declare_field<double>( stk::topology::ELEMENT_RANK, "centroid" );
      m_temperature_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "temperature" );
      m_volume_field      = &m_metaData.declare_field<double>( stk::topology::ELEMENT_RANK, "volume" );

      // Define where fields exist on the mesh:
      stk::mesh::Part & universal = m_metaData.universal_part();

      put_field_on_mesh( *m_coordinates_field , universal, m_spatial_dimension, nullptr);
      put_field_on_mesh( *m_centroid_field , universal, m_spatial_dimension, nullptr);
      put_field_on_mesh( *m_temperature_field, universal, nullptr);
      put_field_on_mesh( *m_volume_field, m_block_wedge, nullptr);
      put_field_on_mesh( *m_volume_field, m_block_tet, nullptr);

      stk::io::put_io_part_attribute(  m_block_wedge );
      stk::io::put_io_part_attribute(  m_block_tet );

      if (do_sidesets)
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

      if (doCommit)
        m_metaData.commit();
    }

    TetWedgeFixture::~TetWedgeFixture()
    { }

    //------------------------------------------------------------------------------
    // Use case specific mesh generation data:

    enum { SpatialDim = 3 };
    enum { node_count = 21 };
    enum { number_wedge = 3 };
    enum { number_tetra = 3 };
    enum { number_quad = 3 };
    enum { number_tri = 2 };

    namespace {

      // Hard coded node coordinate data for all the nodes in the entire mesh
      static const double node_coord_data[ node_count ][ SpatialDim ] = {
        { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,

        { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } ,

        { 3 , 1 , 0 } ,

        { 0 , 2 , 0 } , { 1 , 2 , 0 } ,

        { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
        { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
        { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
        { 1 , 1 , -2 } };

      // Hard coded wedge node ids for all the wedge nodes in the entire mesh
      static const stk::mesh::EntityIdVector wedge_node_ids[number_wedge] {
        { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
        { 20 , 19 , 16,  10 ,  9 ,  6  } ,
        { 16 , 17 , 20 ,  6 ,  7 , 10 } };

      // Hard coded tetra node ids for all the tetra nodes in the entire mesh
      static const stk::mesh::EntityIdVector tetra_node_ids[number_tetra] {
        { 15 , 19 , 16 , 21 } ,
        { 19 , 20 , 16 , 21 } ,
        { 16 , 20 , 17 , 21 } };

#if 0
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

    void TetWedgeFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          stk::mesh::EntityId curr_elem_id = 1;

          // For each element topology declare elements

          stk::mesh::Entity wedges[number_wedge];

          for ( unsigned i = 0 ; i < number_wedge ; ++i , ++curr_elem_id ) {
            wedges[i] = stk::mesh::declare_element( m_bulkData, m_block_wedge, curr_elem_id, wedge_node_ids[i] );
            for ( unsigned jj = 0 ; jj < Wedge6::num_nodes ; ++jj ) {
              unsigned j = wedge_node_ids[i][jj] - 1;
              stk::mesh::Entity const node = m_bulkData.get_entity( stk::topology::NODE_RANK , j + 1 );
              double * const coord = stk::mesh::field_data( *m_coordinates_field , node );
              coord[0] = node_coord_data[j][0] ;
              coord[1] = node_coord_data[j][1] ;
              coord[2] = node_coord_data[j][2] ;
            }
          }

          for ( unsigned i = 0 ; i < number_tetra ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_tet, curr_elem_id, tetra_node_ids[i] );
            for ( unsigned jj = 0 ; jj < Tet4::num_nodes ; ++jj ) {
              unsigned j = tetra_node_ids[i][jj] - 1;
              stk::mesh::Entity const node = m_bulkData.get_entity( stk::topology::NODE_RANK , j + 1 );
              double * const coord = stk::mesh::field_data( *m_coordinates_field , node );
              coord[0] = node_coord_data[j][0] ;
              coord[1] = node_coord_data[j][1] ;
              coord[2] = node_coord_data[j][2] ;
            }
          }

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

        }
      stk::mesh::fixup_ghosted_to_shared_nodes(m_bulkData);
      m_bulkData.modification_end();

    }


  } //namespace percept
