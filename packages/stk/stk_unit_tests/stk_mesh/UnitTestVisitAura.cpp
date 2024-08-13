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

#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/baseImpl/Visitors.hpp>
#include <stk_topology/topology.hpp>

namespace
{

using EntitySet = std::set<stk::mesh::Entity>;

class VisitAura : public stk::unit_test_util::MeshFixture
{
public:
  VisitAura()
    : MeshFixture(3, {"node","edge","face","elem","constraint"}),
      thisProc(stk::parallel_machine_rank(MPI_COMM_WORLD)) {}

  void check_for_expected_entities(const EntitySet& entities,
                                   stk::mesh::EntityRank rank,
                                   const stk::mesh::EntityIdVector& expected)
  {
    for(stk::mesh::EntityId id : expected) {
      stk::mesh::Entity expectedEntity = get_bulk().get_entity(rank, id);
      EXPECT_TRUE(get_bulk().is_valid(expectedEntity));
      EXPECT_TRUE(entities.find(expectedEntity) != entities.end())<<"Failed to find "<<get_bulk().entity_key(expectedEntity);
    }
  }

  EntitySet visit_aura_closure(stk::mesh::EntityRank rank, stk::mesh::EntityId id)
  {
    stk::mesh::Entity entity = get_bulk().get_entity(rank, id);
    EXPECT_TRUE(get_bulk().is_valid(entity));

    stk::mesh::impl::OnlyVisitGhostsOnce onlyVisitGhostsOnce(get_bulk());
    EntitySet entities;
    stk::mesh::impl::StoreInSet<EntitySet> storeInSet(entities);
    stk::mesh::impl::VisitAuraClosureGeneral(get_bulk(), entity, storeInSet, onlyVisitGhostsOnce);
    return entities;
  }

  stk::mesh::Entity create_constraint_entity(stk::mesh::EntityId constraintId,
                                             const stk::mesh::EntityVector& connectedEntities)
  {
    get_bulk().modification_begin();
    stk::mesh::Entity constraintEntity = get_bulk().declare_constraint(constraintId, stk::mesh::ConstPartVector{});
    for(unsigned i=0; i<connectedEntities.size(); ++i) {
      get_bulk().declare_relation(constraintEntity, connectedEntities[i], i);
    }
    get_bulk().modification_end();
    return constraintEntity;
  }

  const int thisProc;
};

static constexpr int numProcs = 2;

TEST_F(VisitAura, auraClosureOfOwnedNode_isEmpty)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != numProcs) { return; }
  setup_mesh("generated:1x1x2|sideset:xXyYzZ", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::EntityId ownedNodeId[numProcs] = {1, 12};

  EntitySet entities = visit_aura_closure(stk::topology::NODE_RANK, ownedNodeId[thisProc]);

  EXPECT_EQ(0u, entities.size());
  check_for_expected_entities(entities, stk::topology::ELEM_RANK, {});
  check_for_expected_entities(entities, stk::topology::FACE_RANK, {});
  check_for_expected_entities(entities, stk::topology::NODE_RANK, {});
}

TEST_F(VisitAura, auraClosureOfSharedNode)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != numProcs) { return; }
  setup_mesh("generated:1x1x2|sideset:xXyYzZ", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::EntityId sharedNodeId[numProcs] = {5, 5};

  EntitySet entities = visit_aura_closure(stk::topology::NODE_RANK, sharedNodeId[thisProc]);

  stk::mesh::EntityIdVector expectedElems[numProcs] = {{2}, {1}};
  stk::mesh::EntityIdVector expectedFaces[numProcs] = {{21,22,23,24,26}, {11,12,13,14,15}};
  stk::mesh::EntityIdVector expectedNodes[numProcs] = {{9,10,11,12}, {1,2,3,4}};

  EXPECT_EQ(10u, entities.size());
  check_for_expected_entities(entities, stk::topology::ELEM_RANK, expectedElems[thisProc]);
  check_for_expected_entities(entities, stk::topology::FACE_RANK, expectedFaces[thisProc]);
  check_for_expected_entities(entities, stk::topology::NODE_RANK, expectedNodes[thisProc]);
}

TEST_F(VisitAura, auraClosureOfAuraElem)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != numProcs) { return; }
  setup_mesh("generated:1x1x2|sideset:xXyYzZ", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::EntityId recvAuraElemId[numProcs] = {2, 1};

  EntitySet entities = visit_aura_closure(stk::topology::ELEM_RANK, recvAuraElemId[thisProc]);

  stk::mesh::EntityIdVector expectedElems[numProcs] = {{2}, {1}};
  stk::mesh::EntityIdVector expectedFaces[numProcs] = {{21,22,23,24,26}, {11,12,13,14,15}};
  stk::mesh::EntityIdVector expectedNodes[numProcs] = {{9,10,11,12}, {1,2,3,4}};

  EXPECT_EQ(10u, entities.size());
  check_for_expected_entities(entities, stk::topology::ELEM_RANK, expectedElems[thisProc]);
  check_for_expected_entities(entities, stk::topology::FACE_RANK, expectedFaces[thisProc]);
  check_for_expected_entities(entities, stk::topology::NODE_RANK, expectedNodes[thisProc]);
}

TEST_F(VisitAura, auraClosureOfAuraFace)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != numProcs) { return; }
  setup_mesh("generated:1x1x2|sideset:xXyYzZ", stk::mesh::BulkData::AUTO_AURA);

  stk::mesh::EntityId recvAuraFaceId[numProcs] = {22, 13};

  EntitySet entities = visit_aura_closure(stk::topology::FACE_RANK, recvAuraFaceId[thisProc]);

  stk::mesh::EntityIdVector expectedElems[numProcs] = {{2}, {1}};
  stk::mesh::EntityIdVector expectedFaces[numProcs] = {{21,22,23,24,26}, {11,12,13,14,15}};
  stk::mesh::EntityIdVector expectedNodes[numProcs] = {{9,10,11,12}, {1,2,3,4}};

  EXPECT_EQ(10u, entities.size());
  check_for_expected_entities(entities, stk::topology::ELEM_RANK, expectedElems[thisProc]);
  check_for_expected_entities(entities, stk::topology::FACE_RANK, expectedFaces[thisProc]);
  check_for_expected_entities(entities, stk::topology::NODE_RANK, expectedNodes[thisProc]);
}

TEST_F(VisitAura, auraClosureOfOwnedConstraintWithElem)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != numProcs) { return; }
  setup_mesh("generated:1x1x2|sideset:xXyYzZ", stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::EntityId constraintIds[numProcs] = {1, 2};
  stk::mesh::EntityId elemIds[numProcs] = {1, 2};
  stk::mesh::Entity ownedElem = get_bulk().get_entity(stk::topology::ELEM_RANK, elemIds[thisProc]);
  create_constraint_entity(constraintIds[thisProc], {ownedElem});
  EntitySet entities = visit_aura_closure(stk::topology::CONSTRAINT_RANK, constraintIds[thisProc]);

  stk::mesh::EntityIdVector expectedConstraints[numProcs] = {{2}, {1}};
  stk::mesh::EntityIdVector expectedElems[numProcs] = {{2}, {1}};
  stk::mesh::EntityIdVector expectedFaces[numProcs] = {{21,22,23,24,26}, {11,12,13,14,15}};
  stk::mesh::EntityIdVector expectedNodes[numProcs] = {{9,10,11,12}, {1,2,3,4}};

  EXPECT_EQ(11u, entities.size());
  check_for_expected_entities(entities, stk::topology::CONSTRAINT_RANK, expectedConstraints[thisProc]);
  check_for_expected_entities(entities, stk::topology::ELEM_RANK, expectedElems[thisProc]);
  check_for_expected_entities(entities, stk::topology::FACE_RANK, expectedFaces[thisProc]);
  check_for_expected_entities(entities, stk::topology::NODE_RANK, expectedNodes[thisProc]);
}

TEST_F(VisitAura, auraClosureOfOwnedConstraintWithNodes)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != numProcs) { return; }
  setup_mesh("generated:1x1x2|sideset:xXyYzZ", stk::mesh::BulkData::AUTO_AURA);
  stk::mesh::EntityId constraintIds[numProcs] = {1, 2};
  stk::mesh::EntityId ownedNodeIds[numProcs] = {1, 12};
  stk::mesh::EntityId sharedNodeIds[numProcs] = {5, 5};
  stk::mesh::Entity ownedNode = get_bulk().get_entity(stk::topology::NODE_RANK, ownedNodeIds[thisProc]);
  stk::mesh::Entity sharedNode = get_bulk().get_entity(stk::topology::NODE_RANK, sharedNodeIds[thisProc]);
  create_constraint_entity(constraintIds[thisProc], {ownedNode, sharedNode});
  EntitySet entities = visit_aura_closure(stk::topology::CONSTRAINT_RANK, constraintIds[thisProc]);

  stk::mesh::EntityIdVector expectedConstraints[numProcs] = {{2}, {1}};
  stk::mesh::EntityIdVector expectedElems[numProcs] = {{2}, {1}};
  stk::mesh::EntityIdVector expectedFaces[numProcs] = {{21,22,23,24,26}, {11,12,13,14,15}};
  stk::mesh::EntityIdVector expectedNodes[numProcs] = {{9,10,11,12}, {1,2,3,4}};

  EXPECT_EQ(11u, entities.size());
  check_for_expected_entities(entities, stk::topology::CONSTRAINT_RANK, expectedConstraints[thisProc]);
  check_for_expected_entities(entities, stk::topology::ELEM_RANK, expectedElems[thisProc]);
  check_for_expected_entities(entities, stk::topology::FACE_RANK, expectedFaces[thisProc]);
  check_for_expected_entities(entities, stk::topology::NODE_RANK, expectedNodes[thisProc]);
}

} // empty namespace

