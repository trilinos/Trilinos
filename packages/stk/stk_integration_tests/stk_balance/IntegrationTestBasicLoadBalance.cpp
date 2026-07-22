#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>

namespace
{

////////////////////////////////////////////////////////////////////////////////////////////

void test_simple_load_balance(stk::mesh::BulkData& stkMeshBulkData);
void move_element_to_other_proc(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityId id);
void balance_mesh(stk::mesh::BulkData& bulkData);
void test_that_mesh_is_balanced(const stk::mesh::BulkData& stkMeshBulkData);
stk::mesh::EntityProcVec get_only_valid_entity_proc_vec(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element);

////////////////////////////////////////////////////////////////////////////////////////////

class BasicLoadBalance : public stk::unit_test_util::MeshFixture {};

TEST_F(BasicLoadBalance, testWithAura)
{
  if (get_parallel_size() != 2) return;
  setup_mesh("generated:1x1x6", stk::mesh::BulkData::AUTO_AURA);
  test_simple_load_balance(get_bulk());
}

TEST_F(BasicLoadBalance, testWithoutAura)
{
  if (get_parallel_size() != 2) return;
  setup_mesh("generated:1x1x6", stk::mesh::BulkData::NO_AUTO_AURA);
  test_simple_load_balance(get_bulk());
}

void test_simple_load_balance(stk::mesh::BulkData& stkMeshBulkData)
{
  move_element_to_other_proc(stkMeshBulkData, 4);
  balance_mesh(stkMeshBulkData);
  test_that_mesh_is_balanced(stkMeshBulkData);
}

void move_element_to_other_proc(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityId id)
{
  stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEM_RANK, id);
  stk::mesh::EntityProcVec element_to_proc = get_only_valid_entity_proc_vec(stkMeshBulkData, element);
  stkMeshBulkData.change_entity_owner(element_to_proc);
}

void balance_mesh(stk::mesh::BulkData& bulkData)
{
  stk::balance::GraphCreationSettingsForZoltan2 graphSettings;
  stk::balance::balanceStkMesh(graphSettings, bulkData);
}

void test_that_mesh_is_balanced(const stk::mesh::BulkData& stkMeshBulkData)
{
  const unsigned gold_num_elements_per_proc = 3;
  std::vector<size_t> counts;
  stk::mesh::count_entities(stkMeshBulkData.mesh_meta_data().locally_owned_part(), stkMeshBulkData, counts);
  EXPECT_EQ(gold_num_elements_per_proc, counts[stk::topology::ELEM_RANK]);
}

stk::mesh::EntityProcVec get_only_valid_entity_proc_vec(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element)
{
  stk::mesh::EntityProcVec element_to_proc;
  if(stkMeshBulkData.is_valid(element) && stkMeshBulkData.bucket(element).owned())
    element_to_proc.push_back(std::make_pair(element, 1 - stkMeshBulkData.parallel_rank()));
  return element_to_proc;
}

}
