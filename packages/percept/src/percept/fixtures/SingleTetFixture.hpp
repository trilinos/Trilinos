// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_SingleTetFixture_hpp
#define percept_SingleTetFixture_hpp

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>


#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>


#include <percept/fixtures/SingleTetFixture.hpp>
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

    class SingleTetFixture {
    public:

      typedef double Point[3];
      typedef stk::mesh::EntityIdVector TetIds;


      ~SingleTetFixture();

      SingleTetFixture( stk::ParallelMachine comm, bool doCommit = true, unsigned npts=0, Point *points=0, unsigned ntets=0, TetIds *tetIds=0,
                        stk::mesh::EntityId elem_id_start=0);

      void populate();

      int m_spatial_dimension;
      std::shared_ptr<stk::mesh::BulkData> m_bulkDataPtr;
      stk::mesh::BulkData& m_bulkData;
      stk::mesh::MetaData& m_metaData;

      stk::mesh::Part & m_block_tet;

      const stk::mesh::EntityRank m_elem_rank;

      CoordinatesFieldType * m_coordinates_field;

      unsigned m_npts;
      Point *m_points;
      unsigned m_ntets;
      TetIds *m_tetIds;

      stk::mesh::EntityId m_elem_id_start;

    };

    //bool verifyMesh( const SingleTetFixture & mesh );

  } //namespace percept

#endif // Stk_Mesh_Use_Cases_UseCase_3_hpp
