#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>

namespace
{

class TestBalanceBalanceMultiPhysics : public stk::unit_test_util::MeshFixture
{
protected:
  TestBalanceBalanceMultiPhysics()
    : MeshFixture(), physics1(nullptr), physics2(nullptr) {}

  void setup_and_test_balance_of_active_only(stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    physics1 = &(get_meta().declare_part("physics1"));
    physics2 = &(get_meta().declare_part("physics2"));
    stk::io::fill_mesh("generated:1x1x30", get_bulk());
    put_elements_1_thru_12_into_physics1();
    put_elements_13_thru_30_into_physics2();
    test_balance_of_multiphysics();
  }

  void put_elements_1_thru_12_into_physics1()
  {
    stk::mesh::EntityId begin_id = 1;
    stk::mesh::EntityId ending_id = 12;
    move_element_range_into_part(begin_id, ending_id, *physics1);
  }

  void put_elements_13_thru_30_into_physics2()
  {
    stk::mesh::EntityId begin_id = 13;
    stk::mesh::EntityId ending_id = 30;
    move_element_range_into_part(begin_id, ending_id, *physics2);
  }

  void move_element_range_into_part(stk::mesh::EntityId begin_id, stk::mesh::EntityId ending_id, stk::mesh::Part& part)
  {
    get_bulk().modification_begin();
    for(stk::mesh::EntityId id=begin_id;id<=ending_id;id++)
      move_only_valid_element_into_part(id, part);
    get_bulk().modification_end();
  }

  void move_only_valid_element_into_part(stk::mesh::EntityId id, stk::mesh::Part& part)
  {
    stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK,id);
    if(get_bulk().is_valid(element) && get_bulk().bucket(element).owned())
      get_bulk().change_entity_parts(element, stk::mesh::ConstPartVector{&part});
  }

  void test_balance_of_multiphysics()
  {
    balance_multi_physics_mesh();
    test_balance_of_mesh();
  }

  void balance_multi_physics_mesh()
  {
    std::vector<stk::mesh::Selector> selectors = {*physics1, *physics2};
    stk::balance::GraphCreationSettingsForZoltan2 graphSettings;
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), selectors);
  }

  void test_balance_of_mesh()
  {
    unsigned numBalanceElements = stk::mesh::count_selected_entities(get_meta().locally_owned_part(), get_bulk().buckets(stk::topology::ELEM_RANK));
    EXPECT_EQ(10u, numBalanceElements);
    test_balance_of_physics(4, *physics1);
    test_balance_of_physics(6, *physics2);
  }

  void test_balance_of_physics(unsigned num_gold_elements_per_proc, stk::mesh::Part& physics)
  {
    unsigned numBalancedElementsPerPhysics = stk::mesh::count_selected_entities(get_meta().locally_owned_part() & physics, get_bulk().buckets(stk::topology::ELEM_RANK));
    EXPECT_EQ(num_gold_elements_per_proc, numBalancedElementsPerPhysics);
  }

  stk::mesh::Part *physics1 = nullptr;
  stk::mesh::Part *physics2 = nullptr;
};

TEST_F(TestBalanceBalanceMultiPhysics, NoAura)
{
  if (stk::parallel_machine_size(get_comm()) == 3) {
    setup_and_test_balance_of_active_only(stk::mesh::BulkData::NO_AUTO_AURA);
    stk::io::write_mesh("final_mesh.g", get_bulk());
  }
}

TEST_F(TestBalanceBalanceMultiPhysics, WithAura)
{
  if(stk::parallel_machine_size(get_comm()) == 3)
    setup_and_test_balance_of_active_only(stk::mesh::BulkData::AUTO_AURA);
}

}


