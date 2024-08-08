// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stddef.h>                     // for size_t, NULL
#include <string.h>                     // for memcpy
#include <unistd.h>                     // for unlink
#include <iostream>                     // for operator<<, etc
#include <map>                          // for map, etc
#include <set>                          // for set, operator==, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities
#include <string>                       // for string, basic_string, etc
#include <utility>                      // for make_pair, pair
#include <vector>                       // for vector
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include "mpi.h"                        // for MPI_Comm_size, etc
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/IossBridge.hpp"        // for is_part_io_part, etc
#include "stk_io/StkMeshIoBroker.hpp"   // for StkMeshIoBroker
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, field_data, etc
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, put_field
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|, etc
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityProc, etc
#include "stk_mesh/base/FEMHelpers.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_mesh/baseImpl/MeshImplUtils.hpp"
#include <stk_tools/mesh_clone/MeshClone.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { namespace fixtures { class BoxFixture; } } }
namespace stk { namespace mesh { namespace fixtures { class RingFixture; } } }
namespace stk { namespace mesh { struct EntityKey; } }
namespace stk { namespace mesh { class FieldBase; } }

namespace
{

void testSubMesh(stk::unit_test_util::BulkDataTester &oldBulkData, stk::mesh::Selector select, int elementToTestId, size_t goldNumberNodes, size_t goldNumberElements)
{
  const stk::mesh::MetaData &oldMetaData = oldBulkData.mesh_meta_data();
  stk::mesh::MetaData newMetaData;
  stk::unit_test_util::BulkDataTester newBulkData(newMetaData, MPI_COMM_SELF);

  stk::tools::copy_mesh(oldBulkData, select, newBulkData);

  const stk::mesh::PartVector &oldParts = oldMetaData.get_parts();
  const stk::mesh::PartVector &newParts = newMetaData.get_parts();
  ASSERT_EQ(oldParts.size(), newParts.size());
  for(size_t i=0; i<oldParts.size(); i++)
  {
    EXPECT_EQ(oldParts[i]->name(), newParts[i]->name());
    EXPECT_EQ(stk::io::is_part_io_part(*oldParts[i]), stk::io::is_part_io_part(*newParts[i]));
  }

  std::vector<size_t> entityCounts;
  stk::mesh::count_entities(newMetaData.universal_part(), newBulkData, entityCounts);
  EXPECT_EQ(goldNumberElements, entityCounts[stk::topology::ELEMENT_RANK]);
  EXPECT_EQ(goldNumberNodes, entityCounts[stk::topology::NODE_RANK]);

  stk::mesh::Entity elementParallel = oldBulkData.get_entity(stk::topology::ELEMENT_RANK, elementToTestId);
  stk::mesh::Entity elementSerial = newBulkData.get_entity(stk::topology::ELEMENT_RANK, elementToTestId);
  EXPECT_EQ(oldBulkData.identifier(elementParallel), newBulkData.identifier(elementSerial));
  unsigned numNodesParallel = oldBulkData.num_nodes(elementParallel);
  unsigned numNodesSerial = newBulkData.num_nodes(elementSerial);
  ASSERT_EQ(numNodesParallel, numNodesSerial);

  const stk::mesh::Entity *nodesParallel = oldBulkData.begin_nodes(elementParallel);
  const stk::mesh::Entity *nodesSerial = newBulkData.begin_nodes(elementSerial);
  std::set<stk::mesh::EntityId> parallelNodeIds;
  std::set<stk::mesh::EntityId> serialNodeIds;
  for(unsigned i = 0; i < numNodesParallel; i++)
  {
    parallelNodeIds.insert(oldBulkData.identifier(nodesParallel[i]));
    serialNodeIds.insert(newBulkData.identifier(nodesSerial[i]));
  }
  EXPECT_TRUE(parallelNodeIds == serialNodeIds);


  const stk::mesh::FieldVector &oldFields = oldMetaData.get_fields();
  const stk::mesh::FieldVector &newFields = newMetaData.get_fields();
  ASSERT_EQ(oldFields.size(), newFields.size());
  for(size_t i=0; i<oldFields.size(); i++)
  {
    EXPECT_EQ(oldFields[i]->name(), newFields[i]->name());
  }
}

TEST(CloningParallelMesh, destinationHasMpiCommSelf_destinationHasAuraEntities)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = -1;
  MPI_Comm_size(comm, &numProcs);
  if(numProcs == 2)
  {
    std::string exodusFileName = "generated:1x1x4|sideset:xXyYzZ|nodeset:xXyYzZ";
    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, comm);

    stk::io::StkMeshIoBroker exodusFileReader(comm);

    exodusFileReader.set_bulk_data(stkMeshBulkData);
    exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    exodusFileReader.create_input_mesh();
    exodusFileReader.populate_bulk_data();

    {
      stk::mesh::Selector sel = stkMeshMetaData.universal_part();
      int elementToTestId = 3;
      size_t goldNumberNodes = 16;
      size_t goldNumberElements = 3;

      testSubMesh(stkMeshBulkData, sel, elementToTestId,  goldNumberNodes, goldNumberElements);
    }

    {
      stk::mesh::Selector sel = stkMeshMetaData.locally_owned_part() | stkMeshMetaData.globally_shared_part();
      int elementToTestId = 2 + stkMeshBulkData.parallel_rank();
      size_t goldNumberNodes = 12;
      size_t goldNumberElements = 2;

      testSubMesh(stkMeshBulkData, sel, elementToTestId,  goldNumberNodes, goldNumberElements);
    }

    {
      stk::mesh::Part *part = stkMeshMetaData.get_part("block_1");
      stk::mesh::Selector sel = *part;
      int elementToTestId = 2 + stkMeshBulkData.parallel_rank();
      size_t goldNumberNodes = 16;
      size_t goldNumberElements = 3;

      testSubMesh(stkMeshBulkData, sel, elementToTestId,  goldNumberNodes, goldNumberElements);
    }
  }
}

}
