/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <stk_percept/fixtures/TetWedgeFixture.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>

#include <stk_io/IossBridge.hpp>

//----------------------------------------------------------------------

namespace stk{
  namespace percept {

    typedef shards::Wedge<6>               Wedge6;
    typedef shards::Tetrahedron<4>         Tet4;

    typedef shards::Quadrilateral<4>  Quad4;
    typedef shards::Triangle<3>       Triangle3;

    TetWedgeFixture::TetWedgeFixture( stk::ParallelMachine comm, bool doCommit, bool do_sidesets ) :
      m_spatial_dimension(3)
      , m_metaData(m_spatial_dimension, stk::mesh::entity_rank_names() )
      , m_bulkData(m_metaData, comm )
      , m_block_wedge(      m_metaData.declare_part< Wedge6 >( "block_2" ))
      , m_block_tet(        m_metaData.declare_part< Tet4 >( "block_3" ))

      , m_sideset_quad(0), m_sideset_quad_subset(0)
      , m_sideset_tri(0), m_sideset_tri_subset(0)

      , m_elem_rank( stk::mesh::MetaData::ELEMENT_RANK )
      , m_coordinates_field( m_metaData.declare_field< VectorFieldType >( "coordinates" ))
      , m_centroid_field(    m_metaData.declare_field< VectorFieldType >( "centroid" ))
      , m_temperature_field( m_metaData.declare_field< ScalarFieldType >( "temperature" ))
      , m_volume_field( m_metaData.declare_field< ScalarFieldType >( "volume" ))
      , m_element_node_coordinates_field( m_metaData.declare_field< ElementNodePointerFieldType >( "elem_node_coord" ))
    {
      // Define where fields exist on the mesh:
      stk::mesh::Part & universal = m_metaData.universal_part();

      put_field( m_coordinates_field , stk::mesh::MetaData::NODE_RANK , universal );
      put_field( m_centroid_field , m_elem_rank , universal );
      put_field( m_temperature_field, stk::mesh::MetaData::NODE_RANK, universal );
      put_field( m_volume_field, m_elem_rank, m_block_wedge );
      put_field( m_volume_field, m_elem_rank, m_block_tet );

      // Define the field-relation such that the values of the
      // 'element_node_coordinates_field' are pointers to the
      // element's nodal 'coordinates_field'.
      // I.e., let:
      //   double *const* elem_node_coord =
      //     field_data( m_element_node_coordinates_field , element );
      // then
      //     elem_node_coord[n][0..2] is the coordinates of element node 'n'
      //     that are attached to that node.

      m_metaData.declare_field_relation(
                                        m_element_node_coordinates_field ,
                                        stk::mesh::get_element_node_stencil(3) ,
                                        m_coordinates_field
                                        );

      // Define element node coordinate field for all element parts
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_wedge, Wedge6::node_count );
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_tet, Tet4::node_count );

      stk::io::put_io_part_attribute(  m_block_wedge );
      stk::io::put_io_part_attribute(  m_block_tet );

      if (do_sidesets)
        {
          m_sideset_quad_subset = &m_metaData.declare_part(std::string("surface_wedge5_quad2d2_1"), stk::mesh::MetaData::FACE_RANK);
          m_sideset_quad =        &m_metaData.declare_part(std::string("surface_1"), stk::mesh::MetaData::FACE_RANK);
          stk::mesh::set_cell_topology< Quad4 >(*m_sideset_quad_subset);
          stk::io::put_io_part_attribute(*m_sideset_quad_subset);
          stk::io::put_io_part_attribute(*m_sideset_quad);
          m_metaData.declare_part_subset(*m_sideset_quad, *m_sideset_quad_subset);

          m_sideset_tri_subset = &m_metaData.declare_part(std::string("surface_wedge5_tri2d2_1"), stk::mesh::MetaData::FACE_RANK);
          m_sideset_tri =        &m_metaData.declare_part(std::string("surface_2"), stk::mesh::MetaData::FACE_RANK);
          stk::mesh::set_cell_topology< Triangle3 >(*m_sideset_tri_subset);
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
      static const stk::mesh::EntityId wedge_node_ids[number_wedge][ Wedge6::node_count ] = {
        { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
        { 20 , 19 , 16,  10 ,  9 ,  6  } ,
        { 16 , 17 , 20 ,  6 ,  7 , 10 } };

      // Hard coded tetra node ids for all the tetra nodes in the entire mesh
      static const stk::mesh::EntityId tetra_node_ids[number_tetra][ Tet4::node_count ] = {
        { 15 , 19 , 16 , 21 } ,
        { 19 , 20 , 16 , 21 } ,
        { 16 , 20 , 17 , 21 } };

      // NOTE: some quad, tri's for wedge sideset testing
      // Hard coded quad node ids for all the quad nodes in the entire mesh
      static const stk::mesh::EntityId quad_node_ids[number_quad][ Quad4::node_count ] = {
        { 5, 9, 19, 15},
        { 7, 17, 20, 10 },
        { 10, 20, 19, 9}
      };

      // wedge element id, side id
      static const stk::mesh::EntityId quad_node_side_ids[number_quad][ 2 ] = {
        {4, 2},
        {6, 1},
        {5, 0}
      };

      // Hard coded tri node ids for all the tri nodes in the entire mesh
      static const stk::mesh::EntityId tri_node_ids[number_tri][ Triangle3::node_count ] = {
        { 5, 6, 9},
        { 6, 10, 9}
      };

      // wedge element id, side id
      static const stk::mesh::EntityId tri_node_side_ids[number_quad][ 2 ] = {
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
            for ( unsigned jj = 0 ; jj < Wedge6::node_count ; ++jj ) {
              unsigned j = wedge_node_ids[i][jj] - 1;
              stk::mesh::Entity const node = m_bulkData.get_entity( stk::mesh::MetaData::NODE_RANK , j + 1 );
              double * const coord = field_data( m_coordinates_field , node );
              coord[0] = node_coord_data[j][0] ;
              coord[1] = node_coord_data[j][1] ;
              coord[2] = node_coord_data[j][2] ;
            }
          }

          for ( unsigned i = 0 ; i < number_tetra ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_tet, curr_elem_id, tetra_node_ids[i] );
            for ( unsigned jj = 0 ; jj < Tet4::node_count ; ++jj ) {
              unsigned j = tetra_node_ids[i][jj] - 1;
              stk::mesh::Entity const node = m_bulkData.get_entity( stk::mesh::MetaData::NODE_RANK , j + 1 );
              double * const coord = field_data( m_coordinates_field , node );
              coord[0] = node_coord_data[j][0] ;
              coord[1] = node_coord_data[j][1] ;
              coord[2] = node_coord_data[j][2] ;
            }
          }

          if (m_sideset_quad)
            {
              for ( unsigned i = 0 ; i < number_quad ; ++i , ++curr_elem_id ) {
                std::cout << "quad i= " << i << std::endl;
                stk::mesh::declare_element_side( m_bulkData,
                                                 curr_elem_id, //side_id,
                                                 wedges[quad_node_side_ids[i][0] - 4], // element,
                                                 quad_node_side_ids[i][1],   //j_side, // local_side_ord,
                                                 m_sideset_quad_subset);
              }
            }

          if (m_sideset_tri)
            {
              for ( unsigned i = 0 ; i < number_tri ; ++i , ++curr_elem_id ) {
                std::cout << "tri i= " << i << std::endl;
                stk::mesh::declare_element_side( m_bulkData,
                                                 curr_elem_id, //side_id,
                                                 wedges[tri_node_side_ids[i][0] - 4], // element,
                                                 tri_node_side_ids[i][1],   //j_side, // local_side_ord,
                                                 m_sideset_tri_subset);
              }
            }

        }
      m_bulkData.modification_end();

    }


  } //namespace percept
} //namespace stk
