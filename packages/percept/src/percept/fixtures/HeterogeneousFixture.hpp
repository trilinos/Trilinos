// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_HeterogeneousFixture_hpp
#define percept_HeterogeneousFixture_hpp

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>


#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>

#include <percept/FieldTypes.hpp>
#include <memory>

/** stk_mesh Use Case 3 - copied and modified here */

#define HET_FIX_INCLUDE_EXTRA_ELEM_TYPES 0

  namespace percept {

    /** Use case with mixed element topologies and
     *  field relations to provide fast access to node field data
     *  from an element.
     *
     *  copied from stk_mesh and modified
     */

    class HeterogeneousFixture {
    public:


      ~HeterogeneousFixture();

      HeterogeneousFixture( stk::ParallelMachine comm, bool doCommit = true, bool do_sidesets=false, bool do_one_sideset=false);

      void populate();

      const int m_spatial_dimension;
      std::shared_ptr<stk::mesh::BulkData> m_bulkDataPtr;
      stk::mesh::BulkData& m_bulkData;
      stk::mesh::MetaData& m_metaData;

      stk::mesh::Part & m_block_hex;
      stk::mesh::Part & m_block_wedge;
      stk::mesh::Part & m_block_tet;
      stk::mesh::Part & m_block_pyramid;
#if HET_FIX_INCLUDE_EXTRA_ELEM_TYPES
      stk::mesh::Part & m_block_quad_shell;
      stk::mesh::Part & m_block_tri_shell;
#endif
      stk::mesh::Part * m_sideset_quad;
      stk::mesh::Part * m_sideset_quad_subset;
      stk::mesh::Part * m_sideset_tri;
      stk::mesh::Part * m_sideset_tri_subset;

      CoordinatesFieldType * m_coordinates_field;
      CoordinatesFieldType * m_centroid_field;
      ScalarFieldType * m_temperature_field;
      ScalarFieldType * m_volume_field;

      unsigned SpatialDim ;
      unsigned node_count ;
      unsigned number_hex ;
      unsigned number_wedge ;
      unsigned number_tetra ;
      unsigned number_pyramid ;
      unsigned number_shell_quad ;
      unsigned number_shell_tri ;
      unsigned number_quad ;
      unsigned number_tri ;

      static const unsigned MAX_NODES = 100;
      static const unsigned MAX_ELEMS = 100;

      double node_coord_data[ MAX_NODES ][ 3 ];
      stk::mesh::EntityIdVector hex_node_ids[MAX_ELEMS];
      stk::mesh::EntityIdVector tetra_node_ids[MAX_ELEMS];
      stk::mesh::EntityIdVector wedge_node_ids[MAX_ELEMS];
      stk::mesh::EntityIdVector pyramid_node_ids[MAX_ELEMS];

      void init() {

        SpatialDim = 3 ;
        node_count = 21 ;
        number_hex = 3 ;
        number_wedge = 3 ;
        number_tetra = 3 ;
        number_pyramid = 2 ;
        number_shell_quad = 3 ;
        number_shell_tri = 3 ;
        number_quad = 3 ;
        number_tri = 2 ;

        // Hard coded node coordinate data for all the nodes in the entire mesh
        double node_coord_data_0[ 21 ][ 3 ] = {
          { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2 , 0 , 0 } , { 3 , 0 , 0 } ,
          { 0 , 1 , 0 } , { 1 , 1 , 0 } , { 2 , 1 , 0 } , { 3 , 1 , 0 } ,
          { 0 , 2 , 0 } , { 1 , 2 , 0 } ,
          { 0 , 0 , -1 } , { 1 , 0 , -1 } , { 2 , 0 , -1 } , { 3 , 0 , -1 } ,
          { 0 , 1 , -1 } , { 1 , 1 , -1 } , { 2 , 1 , -1 } , { 3 , 1 , -1 } ,
          { 0 , 2 , -1 } , { 1 , 2 , -1 } ,
          { 1 , 1 , -2 } };


        // Hard coded hex node ids for all the hex nodes in the entire mesh
        stk::mesh::EntityIdVector hex_node_ids_0[3] {
          { 1 , 2 , 12 , 11 , 5 , 6 , 16 , 15 } ,
            { 2 , 3 , 13 , 12 , 6 , 7 , 17 , 16 } ,
              { 3 , 4 , 14 , 13 , 7 , 8 , 18 , 17 } };

        // Hard coded wedge node ids for all the wedge nodes in the entire mesh
        stk::mesh::EntityIdVector wedge_node_ids_0[3] {
          { 15 , 16 , 19 ,  5 ,  6 ,  9 } ,
            { 20 , 19 , 16 , 10 ,  9 ,  6 } ,
              { 16 , 17 , 20 ,  6 ,  7 , 10 } };

        // Hard coded tetra node ids for all the tetra nodes in the entire mesh
        stk::mesh::EntityIdVector tetra_node_ids_0[3] {
          { 15 , 19 , 16 , 21 } ,
            { 19 , 20 , 16 , 21 } ,
              { 16 , 20 , 17 , 21 } };

        // Hard coded pyramid node ids for all the pyramid nodes in the entire mesh
        stk::mesh::EntityIdVector pyramid_node_ids_0[2] {
          { 11 , 15 , 16 , 12 , 21 } ,
            { 12 , 16 , 17 , 13 , 21 } };

#if 0
        stk::mesh::EntityIdVector quad_node_side_ids_0[3] {
          {4, 2},
            {6, 1},
              {5, 0}
        };


        stk::mesh::EntityIdVector tri_node_side_ids_0[2] {
          {4, 4},
            {5, 3}
        };
#endif
        init (node_count, node_coord_data_0,
              number_hex, hex_node_ids_0,
              number_tetra, tetra_node_ids_0,
              number_wedge, wedge_node_ids_0,
              number_pyramid, pyramid_node_ids_0);

      }

      void init(unsigned nnodes, double coord[][3],
                unsigned nhex, stk::mesh::EntityIdVector *hex_id,
                unsigned ntet, stk::mesh::EntityIdVector *tet_id,
                unsigned nwedge, stk::mesh::EntityIdVector *wedge_id,
                unsigned npyr, stk::mesh::EntityIdVector *pyr_id)
      {
        node_count = nnodes;
        for (unsigned ii=0; ii < node_count; ++ii)
          {
            for (int jj=0; jj < 3; ++jj)
              node_coord_data[ii][jj] = coord[ii][jj];
          }
        number_hex = nhex;
        for (unsigned ii=0; ii < number_hex; ++ii)
          {
            hex_node_ids[ii].resize(8);
            for (int jj=0; jj < 8; ++jj)
              hex_node_ids[ii][jj] = hex_id[ii][jj];
          }
        number_tetra = ntet;
        for (unsigned ii=0; ii < number_tetra; ++ii)
          {
            tetra_node_ids[ii].resize(4);
            for (int jj=0; jj < 4; ++jj)
              tetra_node_ids[ii][jj] = tet_id[ii][jj];
          }
        number_wedge = nwedge;
        for (unsigned ii=0; ii < number_wedge; ++ii)
          {
            wedge_node_ids[ii].resize(6);
            for (int jj=0; jj < 6; ++jj)
              wedge_node_ids[ii][jj] = wedge_id[ii][jj];
          }
        number_pyramid = npyr;
        for (unsigned ii=0; ii < number_pyramid; ++ii)
          {
            pyramid_node_ids[ii].resize(5);
            for (int jj=0; jj < 5; ++jj)
              pyramid_node_ids[ii][jj] = pyr_id[ii][jj];
          }
      }


    };

    bool verifyMesh( const HeterogeneousFixture & mesh );

  } //namespace percept

#endif // Stk_Mesh_Use_Cases_UseCase_3_hpp
