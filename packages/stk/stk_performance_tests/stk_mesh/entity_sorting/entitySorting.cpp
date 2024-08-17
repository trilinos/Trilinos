#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_util/diag/StringUtil.hpp>
#include <stk_mesh/baseImpl/Partition.hpp>
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_mesh/baseImpl/GlobalIdEntitySorter.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>  // for BulkDataTester
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>

namespace
{

class EntitySortingPerformance : public stk::unit_test_util::PerformanceTester
{
public:
    EntitySortingPerformance(stk::unit_test_util::BulkDataTester &bulk)
      : stk::unit_test_util::PerformanceTester(bulk.parallel()),
        bulkData(bulk)
    {
    }

protected:
    virtual void run_algorithm_to_time()
    {
        stk::mesh::impl::BucketRepository& bucketRepository = bulkData.my_get_bucket_repository();
        for (size_t i=0 ; i<numTimesToSort ; ++i)
            sort_bucket_repository(bucketRepository);
    }
    virtual size_t get_value_to_output_as_iteration_count()
    {
        return numTimesToSort;
    }


private:
    stk::unit_test_util::BulkDataTester &bulkData;
    size_t numTimesToSort = 1000;

    void sort_bucket_repository(stk::mesh::impl::BucketRepository& bucketRepository)
    {
        for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEM_RANK; ++rank)
            for(stk::mesh::impl::Partition* partition : bucketRepository.get_partitions(rank))
                partition->sort(stk::mesh::impl::GlobalIdEntitySorter());
    }
};

class SortEntitiesCustomLess : public ::testing::Test
{
protected:
    void run_entity_sort_performance_test(stk::unit_test_util::BulkDataTester& bulk)
    {
        EntitySortingPerformance perfTester(bulk);
        perfTester.run_performance_test();
    }
};

TEST_F(SortEntitiesCustomLess, test_entity_sorting_performance)
{
    stk::mesh::MetaData meta;
    stk::unit_test_util::BulkDataTester bulk(meta, MPI_COMM_WORLD);
    stk::io::fill_mesh("generated:100x100x100",bulk);
    run_entity_sort_performance_test(bulk);
}

} // namespace

