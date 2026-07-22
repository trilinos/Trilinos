#include "gtest/gtest.h"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/EntitySorterBase.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>

namespace {

class EntityReverseSorter : public stk::mesh::EntitySorterBase
{
public:
  virtual void sort(stk::mesh::BulkData &bulk, stk::mesh::EntityVector& entityVector) const
  {
    std::sort(entityVector.begin(), entityVector.end(),
              [&bulk](stk::mesh::Entity a, stk::mesh::Entity b) { return bulk.identifier(a) > bulk.identifier(b); });
  }
};

class HowToSortEntities : public stk::unit_test_util::MeshFixture
{
protected:
  void sort_and_check()
  {
    if(stk::parallel_machine_size(get_comm()) == 1)
    {
      setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
      get_bulk().sort_entities(EntityReverseSorter());
      expect_entities_in_reverse_order();
    }
  }
  void expect_entities_in_reverse_order()
  {
    const stk::mesh::BucketVector buckets = get_bulk().buckets(stk::topology::NODE_RANK);
    ASSERT_EQ(1u, buckets.size());
    expect_bucket_in_reverse_order(*buckets[0]);
  }
  void expect_bucket_in_reverse_order(const stk::mesh::Bucket &bucket)
  {
    ASSERT_EQ(20u, bucket.size());
    for(size_t i=1; i<bucket.size(); i++)
      EXPECT_GT(get_bulk().identifier(bucket[i-1]), get_bulk().identifier(bucket[i]));
  }
};
TEST_F(HowToSortEntities, example_reverse)
{
  sort_and_check();
}

}
