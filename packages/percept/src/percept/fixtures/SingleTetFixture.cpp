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

#include <percept/fixtures/SingleTetFixture.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>


#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/MeshUtils.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
//----------------------------------------------------------------------

  namespace percept {

    SingleTetFixture::SingleTetFixture( stk::ParallelMachine comm, bool doCommit, unsigned npts, Point *points, unsigned ntets, TetIds *tetIds,
                                        stk::mesh::EntityId elem_id_start
                                        ) :
      m_spatial_dimension(3)
      , m_bulkDataPtr(stk::mesh::MeshBuilder(comm).set_spatial_dimension(m_spatial_dimension)
                                                  .create())
      , m_bulkData(*m_bulkDataPtr)
      , m_metaData(m_bulkData.mesh_meta_data())
      , m_block_tet(        m_metaData.declare_part_with_topology( "block_1", stk::topology::TET_4 ))
      , m_elem_rank( stk::topology::ELEMENT_RANK )
      , m_npts(npts), m_points(points)
      , m_ntets(ntets), m_tetIds(tetIds)
      , m_elem_id_start(elem_id_start)
    {
      m_coordinates_field = &m_metaData.declare_field<double>( stk::topology::NODE_RANK, "coordinates" );

      // Define where fields exist on the mesh:
      stk::mesh::Part & universal = m_metaData.universal_part();

      put_field_on_mesh( *m_coordinates_field , universal, m_spatial_dimension, nullptr);

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
      static  SingleTetFixture::Point node_coord_data[ node_count ] = {
        { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 0 , 1 , 0 } , { 0 , 0 , 1 } };

      // Hard coded tetra node ids for all the tetra nodes in the entire mesh
      static  SingleTetFixture::TetIds tetra_node_ids[] = {
        SingleTetFixture::TetIds { 1, 2, 3, 4} };

    }

    //------------------------------------------------------------------------------

    void SingleTetFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          stk::mesh::EntityId curr_elem_id = (m_elem_id_start ? m_elem_id_start : 1);

          // For each element topology declare elements
          unsigned ntets = 1;
          unsigned npts = 4;
          Point *pts = node_coord_data;
          TetIds *tets = tetra_node_ids;

          if (m_ntets && m_npts)
            {
              ntets = m_ntets;
              npts = m_npts;
              pts = m_points;
              tets = m_tetIds;
            }
          for ( unsigned i = 0 ; i < ntets ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_tet, curr_elem_id, tets[i] );
            //std::cout << "tmp SingleTetFixture::populate tets[i]= " << i << " " << tets[i][0] << " " << tets[i][1] << " " << tets[i][2] << " " << tets[i][3] << std::endl;
          }

          // For all nodes assign nodal coordinates
          for ( unsigned i = 0 ; i < npts ; ++i ) {
            stk::mesh::Entity const node = m_bulkData.get_entity( stk::topology::NODE_RANK , i + 1 );
            double * const coord = stk::mesh::field_data( *m_coordinates_field , node );
            coord[0] = pts[i][0] ;
            coord[1] = pts[i][1] ;
            coord[2] = pts[i][2] ;
            //std::cout << "tmp SingleTetFixture::populate coords= " << i << " " << coord[0] << " " << coord[1] << " "  << coord[2] << std::endl;
          }

        }
      stk::mesh::fixup_ghosted_to_shared_nodes(m_bulkData);
      m_bulkData.modification_end();

    }

  } //namespace percept
