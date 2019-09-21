// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_PyramidFixture_hpp
#define percept_PyramidFixture_hpp

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>

#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>

#include <stk_mesh/base/CoordinateSystems.hpp>

#include <percept/fixtures/HeterogeneousFixture.hpp>
#include <percept/FieldTypes.hpp>

/** stk_mesh Use Case 3 - copied and modified here */

#define HET_FIX_INCLUDE_EXTRA_ELEM_TYPES 0

  namespace percept {

    /** Use case with mixed element topologies and
     *  field relations to provide fast access to node field data
     *  from an element.
     *
     *  copied from stk_mesh and modified
     */

    class PyramidFixture {
    public:


      ~PyramidFixture();

      PyramidFixture( stk::ParallelMachine comm, bool doCommit = true, bool do_sidesets = false);

      void populate();

      unsigned SpatialDim;
      unsigned node_count;
      unsigned number_pyramid;
      unsigned number_quad;
      unsigned number_tri;

      static const unsigned MAX_NPYR = 100;

      double node_coord_data[ MAX_NPYR ][ 3 ];

      stk::mesh::EntityIdVector pyramid_node_ids[MAX_NPYR];

      void init()
      {
        SpatialDim = 3 ;
        node_count = 7 ;
        number_pyramid = 2 ;
        number_quad = 2 ;
        number_tri = 6 ;

        // Hard coded node coordinate data for all the nodes in the entire mesh
        double node_coord_data_0[ 7 ][ 3 ] = {

          { 0 , 0 , -1 } , { 1 , 0 , -1 } ,  { 2 , 0 , -1 } ,

          { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } ,

          { 1 , 1 , -2 }
        };

        // Hard coded pyramid node ids for all the pyramid nodes in the entire mesh
        stk::mesh::EntityIdVector pyramid_node_ids_0[2] {
          { 1 , 4 , 5 , 2 , 7 } ,
            { 2 , 5 , 6 , 3 , 7 } };

        init (node_count, number_pyramid, node_coord_data_0, pyramid_node_ids_0);

      }

      void init(unsigned nnodes, unsigned npyr, double coord[][3], stk::mesh::EntityIdVector *nid)
      {
        number_pyramid = npyr;
        node_count = nnodes;
        for (unsigned ii=0; ii < node_count; ++ii)
          {
            for (int jj=0; jj < 3; ++jj)
              node_coord_data[ii][jj] = coord[ii][jj];
          }
        for (unsigned ii=0; ii < number_pyramid; ++ii)
          {
            pyramid_node_ids[ii].resize(5);
            for (int jj=0; jj < 5; ++jj)
              pyramid_node_ids[ii][jj] = nid[ii][jj];
          }
      }


      const int m_spatial_dimension;
      stk::mesh::MetaData m_metaData;
      stk::mesh::BulkData m_bulkData;

      stk::mesh::Part & m_block_pyramid;
      stk::mesh::Part * m_sideset_quad;
      stk::mesh::Part * m_sideset_quad_subset;
      stk::mesh::Part * m_sideset_tri;
      stk::mesh::Part * m_sideset_tri_subset;

      CoordinatesFieldType & m_coordinates_field;
      CoordinatesFieldType & m_centroid_field;
      ScalarFieldType & m_temperature_field;
      ScalarFieldType & m_volume_field;
    };

    bool verifyMesh( const PyramidFixture & mesh );

  } //namespace percept

#endif // Stk_Mesh_Use_Cases_UseCase_3_hpp
