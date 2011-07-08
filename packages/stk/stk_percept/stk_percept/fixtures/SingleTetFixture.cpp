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

#include <stk_percept/fixtures/SingleTetFixture.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

//----------------------------------------------------------------------

namespace stk{
  namespace percept {

    typedef shards::Tetrahedron<4>         Tet4;

    SingleTetFixture::SingleTetFixture( stk::ParallelMachine comm, bool doCommit ) :
      m_spatial_dimension(3)
      , m_metaData(m_spatial_dimension, stk::mesh::fem::entity_rank_names(m_spatial_dimension) )
      , m_bulkData( stk::mesh::fem::FEMMetaData::get_meta_data(m_metaData) , comm )
      , m_block_tet(        m_metaData.declare_part< Tet4 >( "block_1" ))
      , m_elem_rank( m_metaData.element_rank() )
      , m_coordinates_field( m_metaData.declare_field< VectorFieldType >( "coordinates" ))
    {
      // Define where fields exist on the mesh:
      stk::mesh::Part & universal = m_metaData.universal_part();

      put_field( m_coordinates_field , stk::mesh::fem::FEMMetaData::NODE_RANK , universal );

      if (doCommit)
        m_metaData.commit();
    }

    SingleTetFixture::~SingleTetFixture()
    { }

    //------------------------------------------------------------------------------
    // Use case specific mesh generation data:

    enum { SpatialDim = 3 };
    enum { node_count = 4 };

    namespace {

      // Hard coded node coordinate data for all the nodes in the entire mesh
      static const double node_coord_data[ node_count ][ SpatialDim ] = {
        { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 } };

      // Hard coded tetra node ids for all the tetra nodes in the entire mesh
      static const stk::mesh::EntityId tetra_node_ids[1][ 4 ] = {
        { 1, 2, 3, 4} };

    }

    //------------------------------------------------------------------------------

    void SingleTetFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          stk::mesh::EntityId curr_elem_id = 1;

          // For each element topology declare elements

          for ( unsigned i = 0 ; i < 1 ; ++i , ++curr_elem_id ) {
            stk::mesh::fem::declare_element( m_bulkData, m_block_tet, curr_elem_id, tetra_node_ids[i] );
          }

          // For all nodes assign nodal coordinates
          for ( unsigned i = 0 ; i < node_count ; ++i ) {
            stk::mesh::Entity * const node = m_bulkData.get_entity( stk::mesh::fem::FEMMetaData::NODE_RANK , i + 1 );
            double * const coord = field_data( m_coordinates_field , *node );
            coord[0] = node_coord_data[i][0] ;
            coord[1] = node_coord_data[i][1] ;
            coord[2] = node_coord_data[i][2] ;
          }

        }
      m_bulkData.modification_end();

    }

  } //namespace percept
} //namespace stk
