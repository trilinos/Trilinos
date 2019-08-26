#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>

namespace
{

class TestBalanceBalanceActiveEntities : public stk::unit_test_util::MeshFixture
{
protected:
    TestBalanceBalanceActiveEntities()
    : MeshFixture(), activePart(get_meta().declare_part("active")) {}

    void setup_and_test_balance_of_active_only(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        setup_mesh("generated:1x1x6", auraOption);
        stk::unit_test_util::put_mesh_into_part(get_bulk(), activePart);
        make_elements_4_and_5_inactive();
        test_balance_of_active_only();
    }

    void make_elements_4_and_5_inactive()
    {
        mark_element_for_deactivation(4);
        mark_element_for_deactivation(5);
        deactivate_marked_elements();
    }

    void test_balance_of_active_only()
    {
        balance_only_active_mesh();
        test_only_active_is_balanced();
    }

    void balance_only_active_mesh()
    {
        std::vector<stk::mesh::Selector> selectors = {activePart};
        stk::balance::GraphCreationSettingsForZoltan2 graphSettings;
        stk::balance::balanceStkMesh(graphSettings, get_bulk(), selectors);
    }

    void test_only_active_is_balanced()
    {
        unsigned numOwnedActiveElements = stk::mesh::count_selected_entities(get_meta().locally_owned_part() & activePart, get_bulk().buckets(stk::topology::ELEM_RANK));
        EXPECT_EQ(2u, numOwnedActiveElements);
        test_marked_entities_did_not_move();
    }

    void test_marked_entities_did_not_move()
    {
        for(size_t i = 0; i < elementsMarkedForDeactivation.size(); ++i)
            if(get_bulk().parallel_rank() == originalOwningProc[i])
            {
                EXPECT_TRUE(get_bulk().bucket(elementsMarkedForDeactivation[i]).owned());
            }
    }

    void mark_element_for_deactivation(stk::mesh::EntityId id)
    {
        stk::mesh::Entity element = get_bulk().get_entity(stk::topology::ELEM_RANK, id);
        if(get_bulk().is_valid(element) && get_bulk().bucket(element).owned())
        {
            elementsMarkedForDeactivation.push_back(element);
            originalOwningProc.push_back(get_bulk().parallel_rank());
        }
    }

    void deactivate_marked_elements()
    {
        get_bulk().modification_begin();
        for(stk::mesh::Entity& element : elementsMarkedForDeactivation)
            get_bulk().change_entity_parts(element, stk::mesh::ConstPartVector{}, stk::mesh::ConstPartVector{&activePart});
        get_bulk().modification_end();
    }

    stk::mesh::Part &activePart;
    stk::mesh::EntityVector elementsMarkedForDeactivation;
    std::vector<int> originalOwningProc;
};

TEST_F(TestBalanceBalanceActiveEntities, testBasicLoadBalanceWithInactiveNoAura)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
        setup_and_test_balance_of_active_only(stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(TestBalanceBalanceActiveEntities, testBasicLoadBalanceWithInactiveWithAura)
{
    if(stk::parallel_machine_size(get_comm()) == 2)
        setup_and_test_balance_of_active_only(stk::mesh::BulkData::AUTO_AURA);
}

}
