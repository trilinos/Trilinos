#include "gtest/gtest.h"
#include "stk_mesh/baseImpl/DeletedEntityCache.hpp"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_unit_test_utils/TextMesh.hpp"

namespace {

class DeletedEntityCacheTester : public stk::unit_test_util::MeshFixture
{
  protected:
    DeletedEntityCacheTester() :
      MeshFixture(3),
      cache()
    {
      std::string meshDesc =
        "0,1,HEX_8,1,2,3,4,5,6,7,8\n\
         0,2,HEX_8,2,9,10,3,6,11,12,7";
      setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
      stk::unit_test_util::setup_text_mesh(get_bulk(), meshDesc);

      stk::mesh::BucketVector const& buckets = bulkData->get_buckets(stk::topology::NODE_RANK, metaData->universal_part());
      for (auto& bucket : buckets) 
      {
        for (const stk::mesh::Entity& entity : *bucket) 
        {
          nodes.push_back(entity);
          max_local_offset = std::max(entity.local_offset(), max_local_offset);
        }
      }

      cache = std::make_shared<stk::mesh::impl::DeletedEntityCache>(*bulkData);
    }

    std::vector<stk::mesh::Entity::entity_value_type> get_ghost_entity_counts()
    {
      auto& ghostReuseMap = cache->get_ghost_reuse_map();
      std::vector<stk::mesh::Entity::entity_value_type> usedValuesCount(max_local_offset+1, 0);
      for (auto& keyOffsetPair : ghostReuseMap)
      {
        usedValuesCount[keyOffsetPair.second] += 1;
      }

      return usedValuesCount;
    }

  std::shared_ptr<stk::mesh::impl::DeletedEntityCache> cache;
  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::Entity::entity_value_type max_local_offset = 0;

};

}

TEST_F(DeletedEntityCacheTester, mark_entity_as_deleted_nonghost)
{
  if (get_parallel_size() != 1) 
  {
    GTEST_SKIP();
  }

  for (int i=0; i < 3; ++i)
  {
    cache->mark_entity_as_deleted(nodes[i], false);
  }

  auto& ghost_reuse_map = cache->get_ghost_reuse_map();
  EXPECT_EQ(ghost_reuse_map.size(), 0u);

  auto& deleted_entities = cache->get_deleted_entities_current_mod_cycle();
  EXPECT_EQ(deleted_entities.size(), 3u);
  for (int i=0; i < 3; ++i)
  {
    EXPECT_EQ(deleted_entities[i], nodes[i].local_offset());
  }
}

TEST_F(DeletedEntityCacheTester, mark_entity_as_deleted_ghost)
{
  if (get_parallel_size() != 1) 
  {
    GTEST_SKIP();
  }

  for (int i=0; i < 3; ++i)
  {
    cache->mark_entity_as_deleted(nodes[i], true);
  }

  EXPECT_EQ(cache->get_deleted_entities_current_mod_cycle().size(), 0u);
  EXPECT_EQ(cache->get_ghost_reuse_map().size(), 3u);

  auto usedGhosts = get_ghost_entity_counts();
  for (size_t i=0; i < nodes.size(); ++i)
  {
    size_t expected_val = i < 3 ? 1 : 0;
    EXPECT_EQ(usedGhosts[nodes[i].local_offset()], expected_val);
  }
}

TEST_F(DeletedEntityCacheTester, get_entity_for_reuse_initial)
{
  if (get_parallel_size() != 1) 
  {
    GTEST_SKIP();
  }

  EXPECT_EQ(cache->get_entity_for_reuse(), stk::mesh::Entity::InvalidEntity);
}

TEST_F(DeletedEntityCacheTester, update_deleted_entities_container)
{
  if (get_parallel_size() != 1) 
  {
    GTEST_SKIP();
  }
  
  std::vector<stk::mesh::Entity::entity_value_type> destroyedEntities;
  for (int i=0; i < 5; ++i)
  {
    bool isGhost = i < 3;
    cache->mark_entity_as_deleted(nodes[i], isGhost);
    destroyedEntities.push_back(nodes[i].local_offset());
  }
  std::sort(destroyedEntities.begin(), destroyedEntities.end());

  EXPECT_EQ(cache->get_entity_for_reuse(), stk::mesh::Entity::InvalidEntity);
  cache->update_deleted_entities_container();

  std::vector<stk::mesh::Entity::entity_value_type> reusedEntities;
  for (int i=0; i < 5; ++i)
  {
    reusedEntities.push_back(cache->get_entity_for_reuse());
  }
  std::sort(reusedEntities.begin(), reusedEntities.end());

  for (int i=0; i < 5; ++i)
  {
    EXPECT_EQ(destroyedEntities[i], reusedEntities[i]);
  }

  for (int i=0; i < 10; ++i)
  {
    EXPECT_EQ(cache->get_entity_for_reuse(), stk::mesh::Entity::InvalidEntity);
  }
}





