// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_TetWedgeFixture_hpp
#define percept_TetWedgeFixture_hpp

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

/** stk_mesh Use Case 3 - copied and modified here */

#define HET_FIX_INCLUDE_EXTRA_ELEM_TYPES 0

  namespace percept {

    /** Use case with mixed element topologies and
     *  field relations to provide fast access to node field data
     *  from an element.
     *
     *  copied from stk_mesh and modified
     */

    class TetWedgeFixture {
    public:


      ~TetWedgeFixture();

      TetWedgeFixture( stk::ParallelMachine comm, bool doCommit = true, bool do_sidesets=false);

      void populate();

      const int m_spatial_dimension;
      std::shared_ptr<stk::mesh::BulkData> m_bulkDataPtr;
      stk::mesh::BulkData& m_bulkData;
      stk::mesh::MetaData& m_metaData;

      stk::mesh::Part & m_block_wedge;
      stk::mesh::Part & m_block_tet;
      stk::mesh::Part * m_sideset_quad;
      stk::mesh::Part * m_sideset_quad_subset;
      stk::mesh::Part * m_sideset_tri;
      stk::mesh::Part * m_sideset_tri_subset;

      CoordinatesFieldType * m_coordinates_field;
      CoordinatesFieldType * m_centroid_field;
      ScalarFieldType * m_temperature_field;
      ScalarFieldType * m_volume_field;
    };

    bool verifyMesh( const TetWedgeFixture & mesh );

  } //namespace percept

#endif // Stk_Mesh_Use_Cases_UseCase_3_hpp
