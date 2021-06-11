#include <test_utils/MeshFixtureMxNRebalance.hpp>
#include <stk_balance/m2n/MtoNRebalancer.hpp>
#include <stk_balance/m2n/M2NDecomposer.hpp>

#include "stk_unit_test_utils/getOption.h"
#include "stk_mesh/base/Comm.hpp"
#include "stk_util/parallel/ParallelVectorConcat.hpp"
#include "stk_mesh/base/DestroyElements.hpp"

namespace
{

class TestBalanceBalanceSmallToLarge : public MeshFixtureMxNRebalance
{
protected:
    TestBalanceBalanceSmallToLarge() : MeshFixtureMxNRebalance() {}

    virtual unsigned get_x() const { return 3; }
    virtual unsigned get_y() const { return 3; }
    virtual unsigned get_z() const { return 3; }

    virtual unsigned get_num_procs_initial_decomp() const { return 2; }
    virtual unsigned get_num_procs_target_decomp()  const { return 3; }
};

TEST_F(TestBalanceBalanceSmallToLarge, MxN_decompositionWithAura)
{
    if(stk::parallel_machine_size(get_comm()) == static_cast<int>(get_num_procs_initial_decomp()))
        setup_and_test_balance_of_mesh(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(TestBalanceBalanceSmallToLarge, MxN_decompositionWithoutAura)
{
    if(stk::parallel_machine_size(get_comm()) == static_cast<int>(get_num_procs_initial_decomp()))
        setup_and_test_balance_of_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
}

class TestBalanceMtoM : public MeshFixtureMxNRebalance
{
protected:
    TestBalanceMtoM() : MeshFixtureMxNRebalance() {}

    virtual unsigned get_x() const { return 3; }
    virtual unsigned get_y() const { return 3; }
    virtual unsigned get_z() const { return 3; }

    virtual unsigned get_num_procs_initial_decomp() const { return 2; }
    virtual unsigned get_num_procs_target_decomp()  const { return 2; }
};

TEST_F(TestBalanceMtoM, MxM_decompositionWithoutAura)
{
    if(stk::parallel_machine_size(get_comm()) == static_cast<int>(get_num_procs_initial_decomp())) {
        setup_initial_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        stk::balance::M2NBalanceSettings balanceSettings(get_output_filename(), get_num_procs_target_decomp());
        EXPECT_NO_THROW(stk::balance::m2n::rebalanceMtoN(m_ioBroker, balanceSettings));
    }
}

}
