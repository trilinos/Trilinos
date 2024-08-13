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
#include "UnitTestCEOCommonUtils.hpp"   // for add_nodes_to_move, etc
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

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { namespace fixtures { class BoxFixture; } } }
namespace stk { namespace mesh { namespace fixtures { class RingFixture; } } }
namespace stk { namespace mesh { struct EntityKey; } }

namespace stk
{
namespace mesh
{
class FieldBase;
}
}

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fixtures::BoxFixture;

namespace
{

class FieldMgr
{
public:
  FieldMgr(stk::mesh::MetaData &stkMeshMetaData)
    : m_stkMeshMetaData(stkMeshMetaData),
      m_commListNodeFieldName("CommListNode"),
      m_commListNodeField(NULL),
      m_sharingCommMapNodeFieldName("NodeCommInfo"),
      m_sharingCommMapNodeField(NULL),
      m_auraCommMapNodeFieldName("AuraCommMapNode"),
      m_auraCommMapNodeField(NULL),
      m_auraCommMapElementFieldName("ElementCommInfo"),
      m_auraCommMapElementField(NULL)
  {
  }

  ~FieldMgr() {}

  void addCommListNodeField()
  {
    stk::mesh::Field<double> &field = m_stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK,
                                                                              m_commListNodeFieldName, 1);
    m_commListNodeField = &field;
    stk::mesh::put_field_on_mesh(*m_commListNodeField, m_stkMeshMetaData.universal_part(), nullptr);
  }

  void addSharingCommMapNodeField()
  {
    stk::mesh::Field<double> &field = m_stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK,
                                                                              m_sharingCommMapNodeFieldName, 1);
    m_sharingCommMapNodeField = &field;
    stk::mesh::put_field_on_mesh(*m_sharingCommMapNodeField, m_stkMeshMetaData.universal_part(), nullptr);
  }

  void addAuraCommMapNodeField()
  {
    stk::mesh::Field<double> &field = m_stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK,
                                                                              m_auraCommMapNodeFieldName, 1);
    m_auraCommMapNodeField = &field;
    stk::mesh::put_field_on_mesh(*m_auraCommMapNodeField, m_stkMeshMetaData.universal_part(), nullptr);
  }

  void addAuraCommMapElementField()
  {
    stk::mesh::Field<double> &field = m_stkMeshMetaData.declare_field<double>(stk::topology::ELEMENT_RANK,
                                                                              m_auraCommMapElementFieldName, 1);
    m_auraCommMapElementField = &field;
    stk::mesh::put_field_on_mesh(*m_auraCommMapElementField, m_stkMeshMetaData.universal_part(), nullptr);
  }

  stk::mesh::Field<double> * getCommListNodeField() { return m_commListNodeField; }
  stk::mesh::Field<double> * getSharingCommMapNodeField() { return m_sharingCommMapNodeField; }
  stk::mesh::Field<double> * getAuraCommMapNodeField() { return m_auraCommMapNodeField; }
  stk::mesh::Field<double> * getAuraCommMapElementField() { return m_auraCommMapElementField; }

private:
  stk::mesh::MetaData &m_stkMeshMetaData;

  std::string m_commListNodeFieldName;
  stk::mesh::Field<double> *m_commListNodeField;

  std::string m_sharingCommMapNodeFieldName;
  stk::mesh::Field<double> *m_sharingCommMapNodeField;

  std::string m_auraCommMapNodeFieldName;
  stk::mesh::Field<double> *m_auraCommMapNodeField;

  std::string m_auraCommMapElementFieldName;
  stk::mesh::Field<double> *m_auraCommMapElementField;
};

void checkCommMaps(std::string message, stk::unit_test_util::BulkDataTester &stkMeshBulkData, int numElements, bool ownerOfElement[], bool isElementInAuraCommMap[], bool isElementValid[],
                   int numNodes, bool ownerOfNode[], bool isNodeInSharedCommMap[], bool isNodeInAuraCommMap[], bool isNodeValid[]);
void putCommInfoDataOnFields(stk::unit_test_util::BulkDataTester &bulkData, FieldMgr &fieldMgr);

void writeCommInfoFields(stk::unit_test_util::BulkDataTester &bulkData, FieldMgr &fieldMgr, const std::string& filename, double time);

void moveElements2And3ToProc2(stk::unit_test_util::BulkDataTester &stkMeshBulkData);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(UnitTestChangeEntityOwner, changeEntityOwnerCase1)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int numProcs = -1;
  MPI_Comm_size(comm, &numProcs);
  if(numProcs == 3)
  {
    //        std::string exodusFileName = "generated:1x1x6|sideset:xXyYzZ|nodeset:xXyYzZ";
    std::string exodusFileName = "generated:1x1x6";
    const int spatialDim = 3;
    stk::mesh::MetaData stkMeshMetaData(spatialDim);
    stk::unit_test_util::BulkDataTester stkMeshBulkData(stkMeshMetaData, comm);

    stk::io::StkMeshIoBroker exodusFileReader(comm);

    exodusFileReader.set_bulk_data(stkMeshBulkData);
    exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
    exodusFileReader.create_input_mesh();

    FieldMgr fieldMgr(stkMeshMetaData);
    fieldMgr.addCommListNodeField();
    fieldMgr.addSharingCommMapNodeField();
    fieldMgr.addAuraCommMapNodeField();
    fieldMgr.addAuraCommMapElementField();

    exodusFileReader.populate_bulk_data();

    putCommInfoDataOnFields(stkMeshBulkData, fieldMgr);

    double time = 0.5;
    writeCommInfoFields(stkMeshBulkData, fieldMgr, "testBefore.exo",  time);

    moveElements2And3ToProc2(stkMeshBulkData);

    putCommInfoDataOnFields(stkMeshBulkData, fieldMgr);

    writeCommInfoFields(stkMeshBulkData, fieldMgr, "testAfter.exo",  time);

    MPI_Barrier(MPI_COMM_WORLD);

    if(stkMeshBulkData.parallel_rank() == 0)
    {
      unlink("testBefore.exo.3.0");
      unlink("testBefore.exo.3.1");
      unlink("testBefore.exo.3.2");
      unlink("testAfter.exo.3.0");
      unlink("testAfter.exo.3.1");
      unlink("testAfter.exo.3.2");
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void writeCommInfoFields(stk::unit_test_util::BulkDataTester &bulkData, FieldMgr &fieldMgr, const std::string& filename, double time)
{
  stk::io::StkMeshIoBroker ioBroker(bulkData.parallel());
  ioBroker.set_bulk_data(bulkData);

  int indexAfter = ioBroker.create_output_mesh(filename, stk::io::WRITE_RESULTS);
  ioBroker.add_field(indexAfter, *fieldMgr.getCommListNodeField());
  ioBroker.add_field(indexAfter, *fieldMgr.getSharingCommMapNodeField());
  ioBroker.add_field(indexAfter, *fieldMgr.getAuraCommMapNodeField());
  ioBroker.add_field(indexAfter, *fieldMgr.getAuraCommMapElementField());
  ioBroker.write_output_mesh(indexAfter);
  ioBroker.begin_output_step(indexAfter, time);
  ioBroker.write_defined_output_fields(indexAfter);
  ioBroker.end_output_step(indexAfter);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void runProc0(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  {
    int numElements = 6;
    bool ownerOfElement[] = { true, true, false, false, false, false };
    bool isElementInAuraCommMap[] = { false, true, true, false, false, false };
    bool isElementValid[] = { true, true, true, false, false, false };

    int numNodes = 28;
    bool ownerOfNode[] = {
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInSharedCommMap[] = {
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInAuraCommMap[] = {
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeValid[] = {
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false
    };

    checkCommMaps("Before on Proc 0", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                  numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
  }

  stk::mesh::Entity entity = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, 2);
  std::vector<EntityProc> entitiesToMove;
  int proc2=2;
  entitiesToMove.push_back(std::make_pair(entity, proc2));

  CEOUtils::add_nodes_to_move(stkMeshBulkData, entity, proc2, entitiesToMove);

  stkMeshBulkData.change_entity_owner(entitiesToMove);

  {
    int numElements = 6;
    bool ownerOfElement[] = { true, false, false, false, false, false };
    bool isElementInAuraCommMap[] = { true, true, false, false, false, false };
    bool isElementValid[] = { true, true, false, false, false, false };

    int numNodes = 28;
    bool ownerOfNode[] = {
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInSharedCommMap[] = {
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInAuraCommMap[] = {
      true, true, true, true,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeValid[] = {
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false
    };

    checkCommMaps("After on Proc0", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                  numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void runProc1(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  {
    int numElements = 6;
    bool ownerOfElement[] = { false, false, true, true, false, false };
    bool isElementInAuraCommMap[] = { false, true, true, true, true, false };
    bool isElementValid[] = { false, true, true, true, true, false };

    int numNodes = 28;
    bool ownerOfNode[] = {
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInSharedCommMap[] = {
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInAuraCommMap[] = {
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false
    };

    bool isNodeValid[] = {
      false, false, false, false,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false
    };

    checkCommMaps("Before on Proc 1", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                  numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
  }

  stk::mesh::Entity entity = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, 3);
  std::vector<EntityProc> entitiesToMove;
  int proc2=2;
  entitiesToMove.push_back(std::make_pair(entity, proc2));

  CEOUtils::add_nodes_to_move(stkMeshBulkData, entity, proc2, entitiesToMove);

  stkMeshBulkData.change_entity_owner(entitiesToMove);

  {
    int numElements = 6;
    bool ownerOfElement[] = { false, false, false, true, false, false };
    bool isElementInAuraCommMap[] = { false, true, true, true, true, false };
    bool isElementValid[] = { false, false, true, true, true, false };

    int numNodes = 28;
    bool ownerOfNode[] = {
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInSharedCommMap[] = {
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInAuraCommMap[] = {
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false
    };

    bool isNodeValid[] = {
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false
    };

    checkCommMaps("After on Proc 1", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                  numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void runProc2(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  {
    int numElements = 6;
    bool ownerOfElement[] = { false, false, false, false, true, true };
    bool isElementInAuraCommMap[] = { false, false, false, true, true, false };
    bool isElementValid[] = { false, false, false, true, true, true };

    int numNodes = 28;
    bool ownerOfNode[] = {
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      true, true, true, true
    };

    bool isNodeInSharedCommMap[] = {
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInAuraCommMap[] = {
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      true, true, true, true,
      false, false, false, false
    };

    bool isNodeValid[] = {
      false, false, false, false,
      false, false, false, false,
      false, false, false, false,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true
    };

    checkCommMaps("Before on Proc 2", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                  numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
  }

  std::vector<EntityProc> entitiesToMove;

  stkMeshBulkData.change_entity_owner(entitiesToMove);

  {
    int numElements = 6;
    bool ownerOfElement[] = { false, true, true, false, true, true };
    bool isElementInAuraCommMap[] = { true, true, true, true, true, false };
    bool isElementValid[] = { true, true, true, true, true, true };

    int numNodes = 28;
    bool ownerOfNode[] = {
      false, false, false, false,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false,
      true, true, true, true,
      true, true, true, true
    };

    bool isNodeInSharedCommMap[] = {
      false, false, false, false,
      true, true, true, true,
      false, false, false, false,
      true, true, true, true,
      true, true, true, true,
      false, false, false, false,
      false, false, false, false
    };

    bool isNodeInAuraCommMap[] = {
      true, true, true, true,               // left of element 1
      false, false, false, false,           // element 1 : element 2
      true, true, true, true,               // element 2 : element 3
      false, false, false, false,           // element 3 : element 4
      false, false, false, false,           // element 4 : element 5
      true, true, true, true,               // element 5 : element 6
      false, false, false, false            // right of element 6
    };

    bool isNodeValid[] = {
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true,
      true, true, true, true
    };

    checkCommMaps("After on Proc 2", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                  numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void moveElements2And3ToProc2(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  if ( stkMeshBulkData.parallel_rank() == 0 )
  {
    runProc0(stkMeshBulkData);
  }
  else if ( stkMeshBulkData.parallel_rank() == 1 )
  {
    runProc1(stkMeshBulkData);
  }
  else if ( stkMeshBulkData.parallel_rank() == 2 )
  {
    runProc2(stkMeshBulkData);
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void checkCommMaps(std::string message, stk::unit_test_util::BulkDataTester &stkMeshBulkData, int numElements, bool ownerOfElement[], bool isElementInAuraCommMap[], bool isElementValid[],
                   int numNodes, bool ownerOfNode[], bool isNodeInSharedCommMap[], bool isNodeInAuraCommMap[], bool isNodeValid[])
{
  for (int i=0;i<numElements;i++)
  {
    stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, i+1);
    EXPECT_EQ(isElementValid[i], stkMeshBulkData.is_valid(element)) << message << " for element " << i+1;
    if ( isElementValid[i] && stkMeshBulkData.is_valid(element) )
    {
      bool amIOwner = ( stkMeshBulkData.parallel_owner_rank(element) == stkMeshBulkData.parallel_rank() );
      EXPECT_EQ(ownerOfElement[i], amIOwner) << message << " for element " << i+1;
      EXPECT_EQ(isElementInAuraCommMap[i], stkMeshBulkData.is_entity_in_ghosting_comm_map(element)) << message << " for element " << i+1;
    }
  }

  for (int i=0;i<numNodes;i++)
  {
    stk::mesh::Entity node = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, i+1);
    EXPECT_EQ(isNodeValid[i], stkMeshBulkData.is_valid(node)) << message << " for node " << i+1;
    if ( isNodeValid[i] && stkMeshBulkData.is_valid(node) )
    {
      bool amIOwner = ( stkMeshBulkData.parallel_owner_rank(node) == stkMeshBulkData.parallel_rank() );
      EXPECT_EQ(ownerOfNode[i], amIOwner) << message << " for node " << i+1;
      EXPECT_EQ(isNodeInSharedCommMap[i], stkMeshBulkData.my_is_entity_in_sharing_comm_map(node)) << message << " for node " << i+1;
      EXPECT_EQ(isNodeInAuraCommMap[i], stkMeshBulkData.is_entity_in_ghosting_comm_map(node)) << message << " for node " << i+1;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void putCommInfoDataOnFields(stk::unit_test_util::BulkDataTester &bulkData, FieldMgr &fieldMgr)
{
  const stk::mesh::BucketVector &nodeBuckets = bulkData.buckets(stk::topology::NODE_RANK);

  const int ownedValue = 0;
  const int ownedAndSharedValue = 4;
  const int ownedAndGhostedValue = 6;
  const int ownedAndSharedAndGhostedValue = 11;
  const int sharedToThisProcValue = 16;
  const int ghostedToThisProcValue = 20;

  for(size_t i = 0; i < nodeBuckets.size(); i++)
  {
    const stk::mesh::Bucket &bucket = *nodeBuckets[i];
    //(bucket.owned())
    {
      for(size_t j = 0; j < bucket.size(); j++)
      {
        double *commListNode = stk::mesh::field_data(*fieldMgr.getCommListNodeField(), bucket[j]);
        double *sharedCommMapNode = stk::mesh::field_data(*fieldMgr.getSharingCommMapNodeField(), bucket[j]);
        double *auraCommMapNode = stk::mesh::field_data(*fieldMgr.getAuraCommMapNodeField(), bucket[j]);

        // 0: not ghosted
        // 1: ghosted else where
        // 2: ghosted here

        if(bulkData.is_entity_in_ghosting_comm_map(bucket[j]))
        {
          if ( bucket.in_aura() )
          {
            *auraCommMapNode = 2;
          }
          else
          {
            *auraCommMapNode = 1;
          }
        }
        else
        {
          *auraCommMapNode = 0;
        }

        if(bulkData.is_entity_in_ghosting_comm_map(bucket[j]))
        {
          if ( bucket.owned() )
          {
            if(bulkData.is_entity_in_ghosting_comm_map(bucket[j]))
            {
              *sharedCommMapNode = ownedAndSharedAndGhostedValue;
            }
            else
            {
              *sharedCommMapNode = ownedAndSharedValue;
            }
          }
          else
          {
            *sharedCommMapNode = sharedToThisProcValue;
          }
        }
        else if(bulkData.is_entity_in_ghosting_comm_map(bucket[j]))
        {
          if (bucket.in_aura() )
          {
            *sharedCommMapNode = ghostedToThisProcValue;
          }
          else
          {
            *sharedCommMapNode = ownedAndGhostedValue;
          }
        }
        else
        {
          *sharedCommMapNode = ownedValue;
        }

        if(CEOUtils::isEntityValidOnCommList(bulkData, bucket[j]))
        {
          *commListNode = 1;
        }
        else
        {
          *commListNode = 0;
        }
      }
    }
  }

  const stk::mesh::BucketVector &elementBuckets = bulkData.buckets(stk::topology::ELEMENT_RANK);

  for(size_t i = 0; i < elementBuckets.size(); i++)
  {
    const stk::mesh::Bucket &bucket = *elementBuckets[i];
    //if(bucket.owned())
    {
      for(size_t j = 0; j < bucket.size(); j++)
      {
        double *auraCommMap = stk::mesh::field_data(*fieldMgr.getAuraCommMapElementField(), bucket[j]);

        // 0: Not ghosted (owned)
        // 1: ghosted to other proc
        // 2: ghosted here from another proc

        if(bulkData.is_entity_in_ghosting_comm_map(bucket[j]))
        {
          if ( bucket.in_aura() )
          {
            *auraCommMap = ghostedToThisProcValue;
          }
          else
          {
            *auraCommMap = ownedAndGhostedValue;
          }
        }
        else
        {
          *auraCommMap = ownedValue;
        }
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
