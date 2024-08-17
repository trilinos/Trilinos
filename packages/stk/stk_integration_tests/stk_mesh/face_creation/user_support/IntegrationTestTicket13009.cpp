#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/CreateAdjacentEntities.hpp>

class Ticket13009 : public stk::unit_test_util::MeshFixture
{
protected:
  Ticket13009() : MeshFixture(2) {};
  void setup_mesh_with_auto_decomp(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    allocate_bulk(auraOption);
    stk::io::fill_mesh_with_auto_decomp(meshSpecification, get_bulk());
  }

  void run_test_on_procs(int num_procs)
  {
    if(stk::parallel_machine_size(get_comm()) == num_procs)
    {
      setup_mesh_with_auto_decomp("basic.e", stk::mesh::BulkData::AUTO_AURA);

      stk::mesh::PartVector emptyPartVector;
      EXPECT_NO_THROW(stk::mesh::create_adjacent_entities(get_bulk(), emptyPartVector));
    }
  }

};

TEST_F(Ticket13009, check_create_adjacent_entities_serial)
{
  run_test_on_procs(1);
}

TEST_F(Ticket13009, check_create_adjacent_entities_parallel)
{
  run_test_on_procs(2);
}
