/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#include <stk_percept/fixtures/HeterogeneousFixture.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

//using stk::mesh::fem::NODE_RANK;

//----------------------------------------------------------------------

namespace stk{
  namespace percept {

    typedef shards::Hexahedron<8>          Hex8;
    typedef shards::Wedge<6>               Wedge6;
    typedef shards::Tetrahedron<4>         Tet4;
    typedef shards::Pyramid<5>             Pyramid4;

    typedef shards::ShellQuadrilateral<4>  ShellQuad4;
    typedef shards::ShellTriangle<3>       ShellTriangle3;

    HeterogeneousFixture::HeterogeneousFixture( stk::ParallelMachine comm, bool doCommit ) :
      m_spatial_dimension(3)
      , m_metaData(m_spatial_dimension, fem::entity_rank_names(m_spatial_dimension) )
      , m_bulkData( stk::mesh::fem::FEMMetaData::get_meta_data(m_metaData) , comm )
      , m_block_hex(        declare_part< Hex8 >( m_metaData, "block_1" ))
      , m_block_wedge(      declare_part< Wedge6 >(m_metaData,  "block_2" ))
      , m_block_tet(        declare_part< Tet4 >(m_metaData,  "block_3" ))

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      , m_block_pyramid(    declare_part< Pyramid4 >(m_metaData,  "block_4" ))
      , m_block_quad_shell( declare_part< ShellQuad4 >(m_metaData,  "block_5" ))
      , m_block_tri_shell(  declare_part< ShellTriangle3 >(m_metaData,  "block_6" ))
#endif

      , m_elem_rank( m_metaData.element_rank() )
      , m_coordinates_field( m_metaData.declare_field< VectorFieldType >( "coordinates" ))
      , m_centroid_field(    m_metaData.declare_field< VectorFieldType >( "centroid" ))
      , m_temperature_field( m_metaData.declare_field< ScalarFieldType >( "temperature" ))
      , m_volume_field( m_metaData.declare_field< ScalarFieldType >( "volume" ))
      , m_element_node_coordinates_field( m_metaData.declare_field< ElementNodePointerFieldType >( "elem_node_coord" ))
    {
      // Define where fields exist on the mesh:
      Part & universal = m_metaData.universal_part();

      put_field( m_coordinates_field , stk::mesh::fem::NODE_RANK , universal );
      put_field( m_centroid_field , m_elem_rank , universal );
      put_field( m_temperature_field, stk::mesh::fem::NODE_RANK, universal );
      put_field( m_volume_field, m_elem_rank, m_block_hex );
      put_field( m_volume_field, m_elem_rank, m_block_wedge );
      put_field( m_volume_field, m_elem_rank, m_block_tet );

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      put_field( m_volume_field, m_elem_rank, m_block_pyramid );
#endif

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
                                        fem::get_element_node_stencil(3) ,
                                        m_coordinates_field
                                        );

      // Define element node coordinate field for all element parts
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_hex, Hex8::node_count );
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_wedge, Wedge6::node_count );
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_tet, Tet4::node_count );

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_pyramid, Pyramid4::node_count );
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_quad_shell, ShellQuad4::node_count);
      put_field( m_element_node_coordinates_field, m_elem_rank, m_block_tri_shell, ShellTriangle3::node_count );
#endif

      if (doCommit)
        m_metaData.commit();
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

    namespace {

      // Hard coded node coordinate data for all the nodes in the entire mesh
      static const double node_coord_data[ node_count ][ SpatialDim ] = {
        { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,
        { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } , { 3 , 1 , 0 } ,
        { 0 , 2 , 0 } , { 1 , 2 , 0 } ,
        { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
        { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
        { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
        { 1 , 1 , -2 } };

      // Hard coded hex node ids for all the hex nodes in the entire mesh
      static const EntityId hex_node_ids[number_hex][ Hex8::node_count ] = {
        { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
        { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
        { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 } };

      // Hard coded wedge node ids for all the wedge nodes in the entire mesh
      static const EntityId wedge_node_ids[number_wedge][ Wedge6::node_count ] = {
        { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
        { 10 ,  9 ,  6 , 20 , 19 , 16 } ,
        { 16 , 17 , 20 ,  6 ,  7 , 10 } };

      // Hard coded tetra node ids for all the tetra nodes in the entire mesh
      static const EntityId tetra_node_ids[number_tetra][ Tet4::node_count ] = {
        { 15 , 19 , 16 , 21 } ,
        { 19 , 20 , 16 , 21 } ,
        { 16 , 20 , 17 , 21 } };

      // Hard coded pyramid node ids for all the pyramid nodes in the entire mesh
      static const EntityId pyramid_node_ids[number_pyramid][ Pyramid4::node_count ] = {
        { 11 , 15 , 16 , 12 , 21 } ,
        { 12 , 16 , 17 , 13 , 21 } };

      // Hard coded shell quad node ids for all the shell quad nodes in the entire mesh
      static const EntityId shell_quad_node_ids[number_shell_quad][ ShellQuad4::node_count ]={
        { 9 , 6 , 16 , 19 } ,
        { 6 , 7 , 17 , 16 } ,
        { 7 , 8 , 18 , 17 } };

      // Hard coded shell tri node ids for all the shell tri nodes in the entire mesh
      static const EntityId shell_tri_node_ids[number_shell_tri][ ShellTriangle3::node_count ] ={
        { 19 , 16 , 21 } ,
        { 16 , 17 , 21 } ,
        { 17 , 13 , 21 } };

    }

    //------------------------------------------------------------------------------

    void HeterogeneousFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          EntityId curr_elem_id = 1;

          // For each element topology declare elements

          for ( unsigned i = 0 ; i < number_hex ; ++i , ++curr_elem_id ) {
            stk::mesh::fem::declare_element( m_bulkData, m_block_hex, curr_elem_id, hex_node_ids[i] );
          }

          for ( unsigned i = 0 ; i < number_wedge ; ++i , ++curr_elem_id ) {
            stk::mesh::fem::declare_element( m_bulkData, m_block_wedge, curr_elem_id, wedge_node_ids[i] );
          }

          for ( unsigned i = 0 ; i < number_tetra ; ++i , ++curr_elem_id ) {
            stk::mesh::fem::declare_element( m_bulkData, m_block_tet, curr_elem_id, tetra_node_ids[i] );
          }

#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
          for ( unsigned i = 0 ; i < number_pyramid ; ++i , ++curr_elem_id ) {
            stk::mesh::fem::declare_element( m_bulkData, m_block_pyramid, curr_elem_id, pyramid_node_ids[i] );
          }

          for ( unsigned i = 0 ; i < number_shell_quad ; ++i , ++curr_elem_id ) {
            stk::mesh::fem::declare_element( m_bulkData, m_block_quad_shell, curr_elem_id, shell_quad_node_ids[i]);
          }

          for ( unsigned i = 0 ; i < number_shell_tri ; ++i , ++curr_elem_id ) {
            stk::mesh::fem::declare_element( m_bulkData, m_block_tri_shell, curr_elem_id, shell_tri_node_ids[i] );
          }
#endif

          // For all nodes assign nodal coordinates
          for ( unsigned i = 0 ; i < node_count ; ++i ) {
            Entity * const node = m_bulkData.get_entity( stk::mesh::fem::NODE_RANK , i + 1 );
            double * const coord = field_data( m_coordinates_field , *node );
            coord[0] = node_coord_data[i][0] ;
            coord[1] = node_coord_data[i][1] ;
            coord[2] = node_coord_data[i][2] ;
          }

        }
      m_bulkData.modification_end();

    }

    // Verify mesh for 6 different parts
    bool verifyMesh( const HeterogeneousFixture & mesh )
    {
      bool result = true;

      const BulkData & bulkData = mesh.m_bulkData ;
      //const VectorFieldType & node_coord = mesh.m_coordinates_field ;
      //const ElementNodePointerFieldType & elem_node_coord  =  mesh.m_element_node_coordinates_field ;

      std::vector<Bucket *> element_buckets = bulkData.buckets( mesh.m_elem_rank );

      // Create a pair containing Part and matching node_count

      typedef std::pair<Part*, unsigned> PartNodeCountPair;
      std::vector<PartNodeCountPair> part_and_node_counts;
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_hex, Hex8::node_count));
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_wedge, Wedge6::node_count));
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_tet, Tet4::node_count));
#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_pyramid, Pyramid4::node_count));
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_quad_shell, ShellQuad4::node_count));
      part_and_node_counts.push_back(PartNodeCountPair(&mesh.m_block_tri_shell, ShellTriangle3::node_count));
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
        Entity * const node = bulkData.get_entity( stk::mesh::fem::NODE_RANK , i + 1 );
        if ( node == NULL ) {
          std::cerr << "Error!  Invalid null pointer for node returned from "
                    << "bulkData.get_entity( stk::mesh::fem::NODE_RANK, " << i+1 << " ) " << std::endl;
          result = false;
        }
      }

      return result;
    }

  } //namespace percept
} //namespace stk
