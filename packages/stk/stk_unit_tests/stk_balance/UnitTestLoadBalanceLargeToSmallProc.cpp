#include <test_utils/MeshFixtureMxNRebalance.hpp>

namespace
{

class TestBalanceBalanceLargeToSmall : public MeshFixtureMxNRebalance
{
protected:
    TestBalanceBalanceLargeToSmall() : MeshFixtureMxNRebalance() {}

    virtual unsigned get_x() const { return 3; }
    virtual unsigned get_y() const { return 3; }
    virtual unsigned get_z() const { return 4; }

    virtual unsigned get_num_procs_initial_decomp() const { return 3; }
    virtual unsigned get_num_procs_target_decomp()  const { return 2; }
};

TEST_F(TestBalanceBalanceLargeToSmall, MxN_decompositionWithAura)
{
    stk::parallel_machine_barrier(get_comm());
    if(stk::parallel_machine_size(get_comm()) == static_cast<int>(get_num_procs_initial_decomp()))
        setup_and_test_balance_of_mesh(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(TestBalanceBalanceLargeToSmall, MxN_decompositionWithoutAura)
{
    stk::parallel_machine_barrier(get_comm());
    if(stk::parallel_machine_size(get_comm()) == static_cast<int>(get_num_procs_initial_decomp()))
        setup_and_test_balance_of_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
}

}
