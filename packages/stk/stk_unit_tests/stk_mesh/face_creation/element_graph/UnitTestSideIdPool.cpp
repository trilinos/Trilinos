#include <gtest/gtest.h>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/SideIdPool.hpp>
#include <stk_util/parallel/ParallelVectorConcat.hpp>
#include <stk_util/util/SortAndUnique.hpp>

class SideIdPoolTester : public stk::mesh::SideIdPool
{
public:
  SideIdPoolTester(stk::mesh::BulkData &bulk)
    : SideIdPool(bulk)
  {

  }

  const stk::mesh::EntityIdVector &my_get_all_ids()
  {
    return get_all_ids();
  }
};


class SideIdPoolInitialIdsTest : public stk::unit_test_util::MeshFixture
{
protected:
  SideIdPoolInitialIdsTest()
  {
    setup_mesh("generated:1x1x4", stk::mesh::BulkData::AUTO_AURA);
    sideIdPool = new SideIdPoolTester(get_bulk());
  }

  virtual ~SideIdPoolInitialIdsTest()
  {
    delete sideIdPool;
  }

  void generate_1000_ids_on_all_procs()
  {
    sideIdPool->generate_initial_ids(numInitialIdsToRequestPerProc);
    const stk::mesh::EntityIdVector &idsThisProc = sideIdPool->my_get_all_ids();
    stk::parallel_vector_concat(get_comm(), idsThisProc, ids);
  }

  void expect_ids_are_parallel_unique(const unsigned expectedParallelSize)
  {
    size_t sizeBeforeSort = ids.size();
    stk::util::sort_and_unique(ids);
    EXPECT_EQ(sizeBeforeSort, ids.size());
    EXPECT_EQ(expectedParallelSize*get_bulk().parallel_size(), ids.size());
  }

  const unsigned numInitialIdsToRequestPerProc = 1000;
  stk::mesh::EntityIdVector ids;
  SideIdPoolTester *sideIdPool;
};

TEST_F(SideIdPoolInitialIdsTest, ids_are_unique)
{
  generate_1000_ids_on_all_procs();
  expect_ids_are_parallel_unique(2*numInitialIdsToRequestPerProc);
}


class SideIdPoolAdditionalIdTest : public SideIdPoolInitialIdsTest
{
protected:
  void generate_5_additional_ids_on_all_procs()
  {
    sideIdPool->generate_additional_ids_collective(numAdditionalIdsToRequestPerProc);
    const stk::mesh::EntityIdVector &idsThisProc = sideIdPool->my_get_all_ids();
    stk::parallel_vector_concat(get_comm(), idsThisProc, ids);
  }

  void expect_5_additional_ids_are_parallel_unique()
  {
    // when asking for 5 more, it will double the requested size if enough ids are not available
    expect_ids_are_parallel_unique(2*numInitialIdsToRequestPerProc);
  }

  const unsigned numAdditionalIdsToRequestPerProc = 5;
};
TEST_F(SideIdPoolAdditionalIdTest, additional_ids_are_unique)
{
  generate_1000_ids_on_all_procs();
  generate_5_additional_ids_on_all_procs();
  expect_5_additional_ids_are_parallel_unique();
}
