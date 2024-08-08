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
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_tools/mesh_clone/MeshClone.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_io/WriteMesh.hpp>

#include <algorithm>

namespace
{
using stk::unit_test_util::build_mesh;

std::vector<size_t> get_entity_counts(const stk::mesh::BulkData &bulk, const stk::mesh::Selector &selector)
{
  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(bulk, counts, &selector);
  return counts;
}

struct EntityConnectivity
{
  EntityConnectivity(const stk::mesh::BulkData &bulk, stk::mesh::Entity argEntity, stk::mesh::EntityRank relationRank_)
    : entity(argEntity),
      relationRank(relationRank_),
      numConnected(bulk.num_connectivity(entity, relationRank)),
      connected(bulk.begin(entity, relationRank), bulk.begin(entity, relationRank)+numConnected),
      ords(bulk.begin_ordinals(entity, relationRank), bulk.begin_ordinals(entity, relationRank)+numConnected),
      perms()
  {
    if (nullptr != bulk.begin_permutations(entity, relationRank))
    {
      perms.assign(bulk.begin_permutations(entity, relationRank), bulk.begin_permutations(entity, relationRank)+numConnected);
    }
    const bool is_upward_connectivity = relationRank > bulk.entity_rank(entity);
    if (is_upward_connectivity) {
      sort_relations_by_identifier(entity, bulk);
    }
  }

  void sort_relations_by_identifier(stk::mesh::Entity entityArg, const stk::mesh::BulkData& bulk)
  {
    std::sort(connected.begin(), connected.end(), stk::mesh::EntityLess(bulk));

    const stk::mesh::Entity* unsortedEntities = bulk.begin(entityArg, relationRank);
    const stk::mesh::ConnectivityOrdinal* unsortedOrds = bulk.begin_ordinals(entityArg, relationRank);
    const stk::mesh::Permutation* unsortedPerms = bulk.begin_permutations(entityArg, relationRank);

    for(size_t i=0; i<connected.size(); ++i)
    {
      stk::mesh::Entity ent = connected[i];
      for(unsigned j=0; j<numConnected; ++j)
      {
        if (i != j && ent == unsortedEntities[j])
        {
          ords[i] = unsortedOrds[j];
          if (!perms.empty())
          {
            perms[i] = unsortedPerms[j];
          }
          break;
        }
      }
    }
  }

  void remove_unselected(const stk::mesh::Selector& selector, const stk::mesh::BulkData& bulk)
  {
    for(unsigned i=0; i<numConnected; ++i)
    {
      unsigned backwardsIdx = numConnected - i - 1;
      if (!selector(bulk.bucket(connected[backwardsIdx])))
      {
        connected.erase(connected.begin()+backwardsIdx);
        ords.erase(ords.begin()+backwardsIdx);
        if (!perms.empty())
        {
          perms.erase(perms.begin()+backwardsIdx);
        }
      }
    }
    numConnected = connected.size();
  }

  stk::mesh::Entity entity;
  stk::mesh::EntityRank relationRank;
  unsigned numConnected;
  std::vector<stk::mesh::Entity> connected;
  std::vector<stk::mesh::ConnectivityOrdinal> ords;
  std::vector<stk::mesh::Permutation> perms;
};

void expect_equal_entity_connectivity(const stk::mesh::BulkData &oldBulk,
                                      const stk::mesh::BulkData &newBulk,
                                      const EntityConnectivity &oldConnectivity,
                                      const EntityConnectivity &newConnectivity,
                                      stk::mesh::EntityRank rank)
{
  ASSERT_EQ(oldConnectivity.numConnected, newConnectivity.numConnected) << "rank: "<<rank << " oldEntity = " << oldBulk.entity_key(oldConnectivity.entity) << " newEntity = " << newBulk.entity_key(newConnectivity.entity);
  unsigned newConIndex = 0;
  for(unsigned oldConIndex = 0; oldConIndex < oldConnectivity.numConnected; oldConIndex++)
  {
    EXPECT_EQ(oldBulk.identifier(oldConnectivity.connected[oldConIndex]), newBulk.identifier(newConnectivity.connected[newConIndex])) << "rank: "<<rank<<", relation_rank: "<<oldConnectivity.relationRank;
    EXPECT_EQ(oldConnectivity.ords[oldConIndex], newConnectivity.ords[newConIndex]) << "oldConIndex: "<<oldConIndex<<", newConIndex: "<<newConIndex
                                                                                    <<", rank: " << rank << " relation_rank: " << oldConnectivity.relationRank;
    if(!oldConnectivity.perms.empty() && !newConnectivity.perms.empty())
      EXPECT_EQ(oldConnectivity.perms[oldConIndex], newConnectivity.perms[newConIndex]) << "rank: " << rank << " relation_rank: " << oldConnectivity.relationRank;
    else
      EXPECT_TRUE(oldConnectivity.perms.empty() && newConnectivity.perms.empty());

    ++newConIndex;
  }
}

void expect_equal_relations(const stk::mesh::Selector &oldSelector,
                            const stk::mesh::BulkData &oldBulk, stk::mesh::Entity oldEntity,
                            const stk::mesh::BulkData &newBulk, stk::mesh::Entity newEntity,
                            stk::mesh::EntityRank endRank)
{
  stk::mesh::EntityRank rank = oldBulk.entity_rank(oldEntity);
  for(stk::mesh::EntityRank relationRank = stk::topology::NODE_RANK; relationRank < endRank; relationRank++)
  {
    if(relationRank != rank)
    {
      EntityConnectivity oldConnectivity(oldBulk, oldEntity, relationRank);
      oldConnectivity.remove_unselected(oldSelector, oldBulk);
      EntityConnectivity newConnectivity(newBulk, newEntity, relationRank);
      expect_equal_entity_connectivity(oldBulk, newBulk, oldConnectivity, newConnectivity, rank);
    }
  }
}

void expect_equal_field_data(const stk::mesh::BulkData &oldBulk, stk::mesh::Entity oldEntity,
                             const stk::mesh::BulkData &newBulk, stk::mesh::Entity newEntity,
                             const stk::mesh::EntityRank rank)
{
  const stk::mesh::FieldVector &oldFields = oldBulk.mesh_meta_data().get_fields(rank);
  const stk::mesh::FieldVector &newFields = newBulk.mesh_meta_data().get_fields(rank);

  EXPECT_EQ(oldFields.size(), newFields.size());

  for(unsigned i=0; i<oldFields.size(); ++i)
  {
    unsigned oldBytesPerEntity = stk::mesh::field_bytes_per_entity(*oldFields[i], oldEntity);
    unsigned newBytesPerEntity = stk::mesh::field_bytes_per_entity(*newFields[i], newEntity);
    EXPECT_EQ(oldBytesPerEntity, newBytesPerEntity);
    unsigned char* oldData = static_cast<unsigned char*>(stk::mesh::field_data(*oldFields[i], oldEntity));
    unsigned char* newData = static_cast<unsigned char*>(stk::mesh::field_data(*newFields[i], newEntity));
    for(unsigned j=0; j<oldBytesPerEntity; ++j)
      EXPECT_EQ(oldData[j], newData[j]);
  }
}

void expect_equal_entity_counts(const stk::mesh::BulkData &oldBulk, const stk::mesh::Selector &oldSelector, const stk::mesh::BulkData &newBulk)
{
  const stk::mesh::MetaData &oldMeta = oldBulk.mesh_meta_data();
  const stk::mesh::MetaData &newMeta = newBulk.mesh_meta_data();
  stk::mesh::Selector newSelector = newBulk.mesh_meta_data().universal_part();

  ASSERT_EQ(oldBulk.is_automatic_aura_on(), newBulk.is_automatic_aura_on()) << "Aura needs to be the same for this test.";

  std::vector<size_t> oldCounts = get_entity_counts(oldBulk, (oldSelector & oldMeta.locally_owned_part()));
  std::vector<size_t> newCounts = get_entity_counts(newBulk, (newSelector & newMeta.locally_owned_part()));
  ASSERT_EQ(oldCounts.size(), newCounts.size());
  stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(oldMeta.entity_rank_count());
  for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; rank++)
    EXPECT_EQ(oldCounts[rank], newCounts[rank]) << rank;
}

void expect_equal_sharing(const stk::mesh::BulkData& oldBulk, stk::mesh::Entity oldEntity,
                          const stk::mesh::BulkData& newBulk, stk::mesh::Entity newEntity)
{
  std::vector<int> oldSharingProcs;
  std::vector<int> newSharingProcs;
  oldBulk.comm_shared_procs(oldBulk.entity_key(oldEntity), oldSharingProcs);
  newBulk.comm_shared_procs(newBulk.entity_key(newEntity), newSharingProcs);
  ASSERT_EQ(oldSharingProcs.size(), newSharingProcs.size());
  for(size_t i=0; i<oldSharingProcs.size(); ++i) {
    EXPECT_EQ(oldSharingProcs[i], newSharingProcs[i]);
  }
}

void expect_superset_sharing(const stk::mesh::BulkData& oldBulk, stk::mesh::Entity oldEntity,
                             const stk::mesh::BulkData& newBulk, stk::mesh::Entity newEntity)
{
  std::vector<int> oldSharingProcs;
  std::vector<int> newSharingProcs;
  oldBulk.comm_shared_procs(oldBulk.entity_key(oldEntity), oldSharingProcs);
  newBulk.comm_shared_procs(newBulk.entity_key(newEntity), newSharingProcs);
  ASSERT_GE(oldSharingProcs.size(), newSharingProcs.size());

  std::set<int> oldSharingProcsSet(oldSharingProcs.begin(), oldSharingProcs.end());
  std::set<int> newSharingProcsSet(newSharingProcs.begin(), newSharingProcs.end());
  ASSERT_TRUE(
        std::includes(oldSharingProcsSet.begin(), oldSharingProcsSet.end(),
                      newSharingProcsSet.begin(), newSharingProcsSet.end())
        );
}

class CloningMesh : public stk::unit_test_util::MeshFixture
{
protected:
  const char *get_mesh_spec_for_1x1x8_with_sideset() const {return "generated:1x1x8|sideset:x";}
  stk::mesh::BulkData::AutomaticAuraOption get_no_auto_aura_option() const {return stk::mesh::BulkData::NO_AUTO_AURA;}

  CloningMesh()
    : stk::unit_test_util::MeshFixture(3, {"node", "edge", "face", "element", "constraint"})
  { }

  void create_constraints()
  {
    size_t numConstraintsPerProc = 10;
    size_t localConstraintId = 1;
    get_bulk().modification_begin();
    stk::mesh::Entity constraint = get_bulk().declare_constraint(get_bulk().parallel_rank() * numConstraintsPerProc + localConstraintId);
    stk::mesh::Part *submesh = get_meta().get_part("submesh");
    get_bulk().change_entity_parts(constraint, stk::mesh::ConstPartVector{submesh}, {});
    declare_constraint_relations_to_block1_elements(constraint);
    get_bulk().modification_end();
  }

  void declare_constraint_relations_to_block1_elements(stk::mesh::Entity constraint)
  {
    stk::mesh::EntityVector block1Elements;
    stk::mesh::Part *block1Part = get_meta().get_part("block_1");
    stk::mesh::get_selected_entities(get_meta().locally_owned_part() & *block1Part, get_bulk().buckets(stk::topology::ELEMENT_RANK), block1Elements);
    for(size_t i=0; i < block1Elements.size(); i++)
      get_bulk().declare_relation(constraint, block1Elements[i], i);
  }

  void declare_parts()
  {
    get_meta().declare_part_with_topology("block_2", stk::topology::HEX_8);
    get_meta().declare_part("submesh", stk::topology::CONSTRAINT_RANK);
  }

  void move_half_of_mesh_to_block2()
  {
    stk::mesh::Part *block1 = get_meta().get_part("block_1");
    stk::mesh::Part *block2 = get_meta().get_part("block_2");
    stk::mesh::PartVector remove{block1};
    stk::mesh::PartVector add{block2};

    stk::mesh::EntityVector localElements;
    stk::mesh::get_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEMENT_RANK), localElements);

    get_bulk().modification_begin();
    for(unsigned i = localElements.size() / 2; i < localElements.size(); i++)
    {
      get_bulk().change_entity_parts(localElements[i], add, remove);
    }
    get_bulk().modification_end();
  }

  void do_clone_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    stk::mesh::BulkData & oldBulk = get_bulk();
    std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(oldBulk.parallel(), auraOption);
    stk::mesh::BulkData& newBulk = *newBulkPtr;
    stk::mesh::Selector copySelector = get_meta().universal_part();
    stk::tools::copy_mesh(oldBulk, copySelector, newBulk);
    test_cloned_mesh(oldBulk, newBulk);
    expect_equal_entity_part_ordinals(oldBulk, copySelector, newBulk);
  }

  void do_clone_test_with_submesh_selector(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    stk::mesh::Part *submesh = get_meta().get_part("submesh");
    stk::mesh::Part *block1 = get_meta().get_part("block_1");
    stk::mesh::Selector copySelector = *submesh | *block1;
    stk::mesh::BulkData & oldBulk = get_bulk();
    std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(oldBulk.parallel(), auraOption);
    stk::mesh::BulkData& newBulk = *newBulkPtr;
    stk::tools::copy_mesh(oldBulk, copySelector, newBulk);
    test_cloned_submesh(oldBulk, copySelector, newBulk);
  }

  void test_cloned_mesh(stk::mesh::BulkData &oldBulk, const stk::mesh::BulkData &newBulk)
  {
    stk::mesh::Part& universalPart = oldBulk.mesh_meta_data().universal_part();

    expect_equal_entity_counts(oldBulk, universalPart, newBulk);

    const stk::mesh::MetaData &oldMeta = oldBulk.mesh_meta_data();
    ASSERT_EQ(5u, static_cast<unsigned>(oldMeta.entity_rank_count()));
    stk::mesh::Selector localOrSharedSelector = oldMeta.locally_owned_part() | oldMeta.globally_shared_part();
    stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(oldMeta.entity_rank_count());
    stk::mesh::EntityVector oldEntities;
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; rank++)
    {
      stk::mesh::get_selected_entities(localOrSharedSelector, oldBulk.buckets(rank), oldEntities);
      for(stk::mesh::Entity oldEntity : oldEntities)
      {
        stk::mesh::Entity newEntity = newBulk.get_entity(rank, oldBulk.identifier(oldEntity));
        EXPECT_TRUE(newBulk.is_valid(newEntity));
        if (rank < stk::topology::ELEM_RANK)
        {
          expect_equal_sharing(oldBulk, oldEntity, newBulk, newEntity);
        }
        expect_equal_relations(universalPart, oldBulk, oldEntity, newBulk, newEntity, endRank);
        expect_equal_field_data(oldBulk, oldEntity, newBulk, newEntity, rank);
      }
    }
    ASSERT_EQ(oldBulk.has_face_adjacent_element_graph(), newBulk.has_face_adjacent_element_graph());
  }

  bool isEntityConnectedToOwnedSelectedElement(stk::mesh::Entity entity, const stk::mesh::Selector& selector, stk::mesh::BulkData& bulk)
  {
    stk::mesh::Selector alsoOwnedSelector = selector & bulk.mesh_meta_data().locally_owned_part();
    bool foundSelectedElement = false;
    const stk::mesh::Entity* elemIt = bulk.begin_elements(entity);
    const stk::mesh::Entity* elemEndIt = bulk.end_elements(entity);
    for (; (elemIt != elemEndIt) && !foundSelectedElement; ++elemIt)
    {
      if (alsoOwnedSelector(bulk.bucket(*elemIt)))
      {
        foundSelectedElement= true;
      }
    }
    return foundSelectedElement;
  }

  void test_cloned_submesh(stk::mesh::BulkData &oldBulk, const stk::mesh::Selector &copySelector, const stk::mesh::BulkData &newBulk)
  {
    expect_equal_entity_counts(oldBulk, copySelector, newBulk);

    const stk::mesh::MetaData &oldMeta = oldBulk.mesh_meta_data();
    ASSERT_EQ(5u, static_cast<unsigned>(oldMeta.entity_rank_count()));
    stk::mesh::Selector localOrSharedSelector = copySelector & (oldMeta.locally_owned_part() | oldMeta.globally_shared_part());
    stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(oldMeta.entity_rank_count());
    stk::mesh::EntityVector oldEntities;
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; rank++)
    {
      stk::mesh::get_selected_entities(localOrSharedSelector, oldBulk.buckets(rank), oldEntities);
      for(stk::mesh::Entity oldEntity : oldEntities)
      {
        bool shouldBeValidInNewBulk = isEntityConnectedToOwnedSelectedElement(oldEntity, copySelector, oldBulk);
        if (shouldBeValidInNewBulk)
        {
          stk::mesh::Entity newEntity = newBulk.get_entity(rank, oldBulk.identifier(oldEntity));
          ASSERT_TRUE(newBulk.is_valid(newEntity)) << "Entity = " << oldBulk.entity_key(oldEntity);
          if (rank < stk::topology::ELEM_RANK)
          {
            expect_superset_sharing(oldBulk, oldEntity, newBulk, newEntity);
          }
          //expect_equal_relations(copySelector, oldBulk, oldEntity, newBulk, newEntity, endRank); // Faces are different
          expect_equal_field_data(oldBulk, oldEntity, newBulk, newEntity, rank);
        }
      }
    }
    ASSERT_EQ(oldBulk.has_face_adjacent_element_graph(), newBulk.has_face_adjacent_element_graph());
  }

  void expect_equal_entity_part_ordinals(stk::mesh::BulkData &oldBulk, const stk::mesh::Selector &oldSelector, stk::mesh::BulkData &newBulk)
  {
    stk::mesh::Selector newSelector = newBulk.mesh_meta_data().universal_part();
    stk::mesh::EntityRank endRank = static_cast<stk::mesh::EntityRank>(oldBulk.mesh_meta_data().entity_rank_count());
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank < endRank; rank++)
    {
      const stk::mesh::BucketVector &oldBuckets = oldBulk.get_buckets(rank, oldSelector);
      const stk::mesh::BucketVector &newBuckets = newBulk.get_buckets(rank, newSelector);
      ASSERT_EQ(oldBuckets.size(), newBuckets.size());
      for(size_t b = 0; b < oldBuckets.size(); b++)
      {
        const stk::mesh::Bucket &oldBucket = *oldBuckets[b];
        const stk::mesh::Bucket &newBucket = *newBuckets[b];
        ASSERT_EQ(oldBucket.size(), newBucket.size());

        const stk::mesh::PartVector &oldParts = oldBucket.supersets();
        const stk::mesh::PartVector &newParts = newBucket.supersets();
        ASSERT_EQ(oldParts.size(), newParts.size());
        for(size_t p = 0; p < oldParts.size(); p++)
          EXPECT_EQ(oldParts[p]->mesh_meta_data_ordinal(), newParts[p]->mesh_meta_data_ordinal());
      }
    }
  }
};

TEST_F(CloningMesh, cloningBulkDataNoAura_selector)
{
  setup_empty_mesh(get_no_auto_aura_option());
  declare_parts();
  stk::io::fill_mesh(get_mesh_spec_for_1x1x8_with_sideset(), get_bulk());
  move_half_of_mesh_to_block2();
  create_constraints();
  do_clone_test_with_submesh_selector(get_no_auto_aura_option());
}

TEST_F(CloningMesh, cloningBulkDataWithAura_selector)
{
  stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA;
  setup_empty_mesh(auraOption);
  declare_parts();
  stk::io::fill_mesh(get_mesh_spec_for_1x1x8_with_sideset(), get_bulk());
  move_half_of_mesh_to_block2();
  create_constraints();
  do_clone_test_with_submesh_selector(auraOption);
}

TEST_F(CloningMesh, cloningBulkDataWithFacesAndAura_selector)
{
  stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA;
  setup_empty_mesh(auraOption);
  declare_parts();
  stk::io::fill_mesh(get_mesh_spec_for_1x1x8_with_sideset(), get_bulk());
  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), stk::mesh::PartVector{}, true);
  move_half_of_mesh_to_block2();
  create_constraints();
  do_clone_test_with_submesh_selector(auraOption);
}

TEST_F(CloningMesh, cloningBulkDataNoAura_same)
{
  setup_mesh(get_mesh_spec_for_1x1x8_with_sideset(), get_no_auto_aura_option());
  do_clone_test(get_no_auto_aura_option());
}

TEST_F(CloningMesh, cloningBulkDataWithAura_same)
{
  stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA;
  setup_mesh(get_mesh_spec_for_1x1x8_with_sideset(), auraOption);
  do_clone_test(auraOption);
}

TEST_F(CloningMesh, cloningBulkDataWithFacesAndAura_same)
{
  stk::mesh::BulkData::AutomaticAuraOption auraOption = stk::mesh::BulkData::AUTO_AURA;
  setup_mesh(get_mesh_spec_for_1x1x8_with_sideset(), auraOption);
  stk::mesh::create_all_sides(get_bulk(), get_meta().universal_part(), stk::mesh::PartVector{}, true);
  do_clone_test(auraOption);
}

TEST_F(CloningMesh, cloningBulkData_modifiableStateThrows)
{
  setup_mesh(get_mesh_spec_for_1x1x8_with_sideset(), get_no_auto_aura_option());

  stk::mesh::BulkData & oldBulk = get_bulk();
  std::shared_ptr<stk::mesh::BulkData> newBulkPtr = build_mesh(oldBulk.parallel(), get_no_auto_aura_option());
  stk::mesh::BulkData& newBulk = *newBulkPtr;
  oldBulk.modification_begin();

  ASSERT_THROW(stk::tools::copy_mesh(oldBulk, get_meta().universal_part(), newBulk), std::logic_error);
}

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
//portability issue: bulk-data size seems to be different (48 bytes smaller) on intel 15.0
//running this test on gcc should be sufficient to detect addition/deletion of class members.

TEST(BulkDataSize, sizeChanges_needToUpdateCopyMesh)
{
  // KHP: different compilers (and different version of a compiler) give different values for the sizeof bulk.
  // Only test on gcc 4.9.3 for now.
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 9))
  std::shared_ptr<stk::mesh::BulkData> bulk = build_mesh(3, MPI_COMM_WORLD);
  EXPECT_TRUE(1176u >= sizeof(*bulk)) << "Size of BulkData changed.  Does mesh copying capability need to be updated?";
#endif
}

#endif

}
