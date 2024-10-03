#include "gtest/gtest.h"                // for AssertHelper, EXPECT_EQ, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_util/diag/StringUtil.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/PerformanceTester.hpp>

namespace
{

class GetPartsByNamePerformance : public stk::unit_test_util::PerformanceTester
{
public:
  GetPartsByNamePerformance(stk::mesh::BulkData &bulk)
    : stk::unit_test_util::PerformanceTester(bulk.parallel()),
      bulkData(bulk)
  {
  }

protected:
  virtual void run_algorithm_to_time()
  {
    const unsigned numParts = 1e6;
    for(unsigned i=0; i<numParts; i++)
    {
      const std::string name = "part_" + sierra::to_string(i);
      bulkData.mesh_meta_data().declare_part(name);
      const stk::mesh::Part *part = bulkData.mesh_meta_data().get_part(name);
      EXPECT_TRUE(part != nullptr);
    }
  }
  virtual size_t get_value_to_output_as_iteration_count()
  {
    return bulkData.mesh_meta_data().get_parts().size();
  }

  stk::mesh::BulkData &bulkData;
};

class GetPartsByName : public stk::unit_test_util::MeshFixture
{
protected:
  void run_get_parts_perf_test()
  {
    GetPartsByNamePerformance perfTester(get_bulk());
    perfTester.run_performance_test();
  }
};

TEST_F(GetPartsByName, test_get_parts)
{

  setup_mesh("generated:1x1x16", stk::mesh::BulkData::AUTO_AURA);
  run_get_parts_perf_test();
}

}

