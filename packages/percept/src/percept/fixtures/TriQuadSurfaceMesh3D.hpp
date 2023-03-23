// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_TriQuadSurfaceMesh3D_hpp
#define percept_TriQuadSurfaceMesh3D_hpp

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>


#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>


#include <percept/PerceptBoostArray.hpp>
#include <percept/FieldTypes.hpp>
#include <memory>

  namespace percept {

    /** Use case with mixed element topologies and
     *  field relations to provide fast access to node field data
     *  from an element.
     *
     *  copied from stk_mesh and modified
     */

    class TriQuadSurfaceMesh3D {
    public:

      //typedef double Point[3];
      typedef std::array<double,3> Point;
      typedef stk::mesh::EntityIdVector QuadIds;
      typedef stk::mesh::EntityIdVector TriIds;

      ~TriQuadSurfaceMesh3D();

      TriQuadSurfaceMesh3D( stk::ParallelMachine comm, bool doCommit = true);

      void populate(unsigned nNodes, unsigned nTriFaces, unsigned nQuadFaces, Point *coords, TriIds *triFaces, QuadIds *quadFaces,
                    stk::mesh::EntityId *nodeIdMap = 0);

      std::shared_ptr<stk::mesh::BulkData> m_bulkDataPtr;
      stk::mesh::BulkData& m_bulkData;
      stk::mesh::MetaData& m_metaData;

      stk::mesh::Part * m_sideset;
      stk::mesh::Part& m_block_quad_shell;
      stk::mesh::Part& m_block_tri_shell;

      CoordinatesFieldType * m_coordinates_field;
      // CoordinatesFieldType * m_centroid_field;
      // ScalarFieldType * m_temperature_field;
      // ScalarFieldType * m_volume_field;
    };


  } //namespace percept

#endif // Stk_Mesh_Use_Cases_UseCase_3_hpp
