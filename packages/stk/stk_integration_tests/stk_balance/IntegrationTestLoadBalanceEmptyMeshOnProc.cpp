#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>

namespace
{

////////////////////////////////////////////////////////////////////////////////////////////

void balance_mesh(stk::mesh::BulkData& bulkData);
void test_that_mesh_is_balanced(const stk::mesh::BulkData& stkMeshBulkData, const unsigned gold_num_elements_per_proc);

////////////////////////////////////////////////////////////////////////////////////////////

void test_load_balancing_when_one_proc_has_no_mesh(stk::mesh::BulkData& stkMeshBulkData);
void move_elements_from_proc_2_to_proc_0(stk::mesh::BulkData& bulkData);
void test_that_empty_mesh_exists_on_proc_2(const stk::mesh::BulkData& stkMeshBulkData);

class EmptyMeshOnProc : public stk::unit_test_util::MeshFixture {};

TEST_F(EmptyMeshOnProc, testEmptyMeshOnProcNoAura)
{
    if (get_parallel_size() != 3) return;
    setup_mesh("generated:1x1x6", stk::mesh::BulkData::NO_AUTO_AURA);
    test_load_balancing_when_one_proc_has_no_mesh(get_bulk());
}

TEST_F(EmptyMeshOnProc, testEmptyMeshOnProcWithAura)
{
    if (get_parallel_size() != 3) return;
    setup_mesh("generated:1x1x6", stk::mesh::BulkData::AUTO_AURA);
    test_load_balancing_when_one_proc_has_no_mesh(get_bulk());
}

void test_load_balancing_when_one_proc_has_no_mesh(stk::mesh::BulkData& stkMeshBulkData)
{
    move_elements_from_proc_2_to_proc_0(stkMeshBulkData);
    test_that_empty_mesh_exists_on_proc_2(stkMeshBulkData);
    balance_mesh(stkMeshBulkData);
    test_that_mesh_is_balanced(stkMeshBulkData, 2);
}

void move_elements_from_proc_2_to_proc_0(stk::mesh::BulkData& bulkData)
{
    stk::mesh::EntityVector local_elements;
    stk::mesh::get_selected_entities(bulkData.mesh_meta_data().locally_owned_part(),
            bulkData.buckets(stk::topology::ELEM_RANK), local_elements);

    int dest_proc = bulkData.parallel_rank() == 2 ? 0 : bulkData.parallel_rank();
    stk::mesh::EntityProcVec decomp(local_elements.size());
    for(size_t i=0;i<local_elements.size();++i)
    {
        decomp[i].first = local_elements[i];
        decomp[i].second = dest_proc;
    }
    stk::balance::internal::rebalance(bulkData, decomp);
}

void test_that_empty_mesh_exists_on_proc_2(const stk::mesh::BulkData& stkMeshBulkData)
{
    std::vector<unsigned> num_elements_per_proc = { 4, 2, 0 };
    unsigned num_elements = stk::mesh::count_selected_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(),
            stkMeshBulkData.buckets(stk::topology::ELEM_RANK));
    EXPECT_EQ(num_elements_per_proc[stkMeshBulkData.parallel_rank()], num_elements);
}

void balance_mesh(stk::mesh::BulkData& bulkData)
{
    stk::balance::BasicZoltan2Settings graphSettings;
    stk::balance::balanceStkMesh(graphSettings, bulkData);
}

void test_that_mesh_is_balanced(const stk::mesh::BulkData& stkMeshBulkData, const unsigned gold_num_elements_per_proc)
{
    std::vector<size_t> counts;
    stk::mesh::count_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData, counts);
    EXPECT_EQ(gold_num_elements_per_proc, counts[stk::topology::ELEM_RANK]);
}

}
