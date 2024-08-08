// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


/**
 */

#include <percept/fixtures/BeamFixture.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
//----------------------------------------------------------------------

  namespace percept {

    typedef stk::topology::topology_type<stk::topology::BEAM_2> stkBeam2;

    BeamFixture::BeamFixture( stk::ParallelMachine comm, bool doCommit ) :
      m_spatial_dimension(3)
      , m_bulkDataPtr(stk::mesh::MeshBuilder(comm).set_spatial_dimension(m_spatial_dimension)
                                                .create())
      , m_bulkData(*m_bulkDataPtr)
      , m_metaData(m_bulkData.mesh_meta_data())
      , m_block_beam( m_metaData.declare_part_with_topology( "block_2", stk::topology::BEAM_2 ) )
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
      put_field_on_mesh( *m_volume_field, m_block_beam, nullptr);

      if (doCommit)
        m_metaData.commit();
    }

    BeamFixture::~BeamFixture()
    { }

    //------------------------------------------------------------------------------
    // Use case specific mesh generation data:

    enum { SpatialDim = 3 };
    enum { node_count = 4 };
    enum { number_beam = 2 };

    namespace {

      // Hard coded node coordinate data for all the nodes in the entire mesh
      static const double node_coord_data[ node_count ][ SpatialDim ] = {
        { 10 , 10 , 10 } , { 11 , 12 , 13 } , { 24 , 24 , 24 } , { 28 , 29 , 30 }
      };

      // Hard coded beam node ids for all the beam nodes in the entire mesh
      static const stk::mesh::EntityIdVector beam_node_ids[number_beam] {
        { 1, 2 },
        { 3, 4 }
      };

    }

    //------------------------------------------------------------------------------

    void BeamFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          stk::mesh::EntityId curr_elem_id = 1;

          // For each element topology declare elements
          for ( unsigned i = 0 ; i < number_beam ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_beam, curr_elem_id, beam_node_ids[i] );
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
    bool verifyMesh( const BeamFixture & mesh )
    {
      bool result = true;

      return result;
    }

  } //namespace percept
