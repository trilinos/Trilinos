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

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/EntityProcMapping.hpp>
#include <stk_io/FillMesh.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>

#include <vector>
#include <set>

using stk::unit_test_util::build_mesh;

TEST(EntityProcMapping, basic)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) > 3) {
    return;
  }

  const unsigned spatialDim = 3;
  std::shared_ptr<stk::mesh::BulkData> bulkPtr = build_mesh(spatialDim, MPI_COMM_WORLD);
  stk::mesh::BulkData& bulk = *bulkPtr;
  stk::io::fill_mesh("generated:1x2x3",bulk);

  if (stk::parallel_machine_rank(MPI_COMM_WORLD) != 0) {
    return;
  }

  stk::mesh::EntityProcMapping entProcMapping(bulk.get_size_of_entity_index_space());
  stk::mesh::Entity elem1 = bulk.get_entity(stk::topology::ELEM_RANK, 1);
  int proc = 0;
  entProcMapping.addEntityProc(elem1, proc);
  proc = 1;
  entProcMapping.addEntityProc(elem1, proc);
  proc = 0;
  entProcMapping.addEntityProc(elem1, proc);

  stk::mesh::Entity elem2 = bulk.get_entity(stk::topology::ELEM_RANK, 2);
  proc = 0;
  entProcMapping.addEntityProc(elem2, proc);
  proc = 1;
  entProcMapping.addEntityProc(elem2, proc);
  proc = 0;
  entProcMapping.addEntityProc(elem2, proc);

  EXPECT_EQ(2u, entProcMapping.get_num_entities());
  EXPECT_EQ(2u, entProcMapping.get_num_procs(elem1));
  EXPECT_EQ(2u, entProcMapping.get_num_procs(elem2));

  EXPECT_TRUE(entProcMapping.find(stk::mesh::EntityProc(elem1, proc)));
  proc = 5;
  EXPECT_FALSE(entProcMapping.find(stk::mesh::EntityProc(elem1, proc)));

  stk::mesh::EntityLess entLess(bulk);
  std::set<stk::mesh::EntityProc, stk::mesh::EntityLess> entityProcSet(entLess);
  entProcMapping.fill_set(entityProcSet);
  EXPECT_EQ(4u, entityProcSet.size());

  std::vector<stk::mesh::EntityProc> entityProcVec;
  entProcMapping.fill_vec(entityProcVec);
  EXPECT_EQ(4u, entityProcVec.size());

  proc = 1;
  entProcMapping.eraseEntityProc(elem1, proc);

  entProcMapping.fill_set(entityProcSet);
  EXPECT_EQ(3u, entityProcSet.size());

  entProcMapping.fill_vec(entityProcVec);
  EXPECT_EQ(3u, entityProcVec.size());
}

void test_add_two_remove_one_then_other_still_found(stk::mesh::EntityProcMapping& mapping, stk::mesh::Entity entity)
{
  mapping.addEntityProc(entity, 0);
  mapping.addEntityProc(entity, 2);
  EXPECT_TRUE(mapping.find(entity,0));
  EXPECT_TRUE(mapping.find(entity,2));

  mapping.eraseEntityProc(entity,2);
  EXPECT_TRUE(mapping.find(entity,0));
  EXPECT_FALSE(mapping.find(entity,2));
}

TEST(EntityProcMapping, add_two_remove_one_then_other_still_found)
{
  stk::mesh::Entity entity(1);
  const unsigned arbitraryMaxNumEntities = 10;
  stk::mesh::EntityProcMapping mapping(arbitraryMaxNumEntities);
  test_add_two_remove_one_then_other_still_found(mapping, entity);
}

TEST(EntityProcMapping, add_two_remove_one_then_other_still_found_with_reset)
{
  stk::mesh::Entity entity(1);
  const unsigned arbitraryMaxNumEntities = 10;
  stk::mesh::EntityProcMapping mapping(arbitraryMaxNumEntities);
  test_add_two_remove_one_then_other_still_found(mapping, entity);

  const unsigned largerMaxNumEntities = 128;
  mapping.reset(largerMaxNumEntities);
  EXPECT_FALSE(mapping.find(entity,0));
  EXPECT_FALSE(mapping.find(entity,2));
  EXPECT_FALSE(mapping.find(entity));

  test_add_two_remove_one_then_other_still_found(mapping, entity);
}

TEST(EntityProcMapping, erase_nonexisting_then_previous_proc_still_found)
{
  stk::mesh::Entity entity(1);
  const unsigned arbitraryMaxNumEntities = 10;
  stk::mesh::EntityProcMapping mapping(arbitraryMaxNumEntities);
  mapping.addEntityProc(entity, 0);
  EXPECT_TRUE(mapping.find(entity,0));

  mapping.eraseEntityProc(entity,2);
  EXPECT_TRUE(mapping.find(entity,0));
}

TEST(EntityProcMapping, visitEntityProcs)
{
  stk::mesh::Entity entity1(1), entity2(2);
  const unsigned arbitraryMaxNumEntities = 10;
  stk::mesh::EntityProcMapping mapping(arbitraryMaxNumEntities);
  mapping.addEntityProc(entity1, 2);
  mapping.addEntityProc(entity2, 1);
  mapping.addEntityProc(entity2, 3);

  std::vector<stk::mesh::EntityProc> gold = {stk::mesh::EntityProc(entity1,2),stk::mesh::EntityProc(entity2,1),stk::mesh::EntityProc(entity2,3)};
  std::vector<stk::mesh::EntityProc> entityProcs;
  mapping.visit_entity_procs([&](stk::mesh::Entity entity, int proc){entityProcs.push_back(stk::mesh::EntityProc(entity,proc));});
  EXPECT_EQ(gold, entityProcs);
}

