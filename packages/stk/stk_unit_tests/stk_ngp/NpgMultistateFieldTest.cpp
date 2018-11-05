#include <stk_ngp/Ngp.hpp>
#include <stk_ngp/NgpMultistateField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/stk_config.h>


template <typename Field>
void assign_value_to_field(stk::mesh::BulkData &bulk, Field ngpField, double value)
{
    ngp::Mesh ngpMesh(bulk);
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part(), KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngpField.get(entity, 0) = value;
    });
}
void assign_value_to_statedfield(stk::mesh::BulkData &bulk, ngp::ConvenientMultistateField<double> ngpStatedField, double value)
{
    ngp::Mesh ngpMesh(bulk);
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part(), KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        ngpStatedField.get_new(entity, 0) = value;
    });
}

namespace {

class StatedFields : public stk::unit_test_util::MeshFixture
{
protected:
    ngp::ConvenientMultistateField<double> create_field_with_num_states(unsigned numStates)
    {
        make_mesh_with_field(numStates);
        return ngp::ConvenientMultistateField<double>(get_bulk(), *stkField);
    }
    void make_mesh_with_field(unsigned numStates)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        stkField = &get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "myField", numStates);
        double init = -1.0;
        stk::mesh::put_field_on_mesh(*stkField, get_meta().universal_part(), &init);
        stk::io::fill_mesh("generated:1x1x4", get_bulk());
    }
    void test_n_states(unsigned numStates)
    {
        ngp::MultistateField<double> ngpField = create_field_with_num_states(numStates);
        EXPECT_EQ(numStates, ngpField.get_num_states());
    }
    void verify_initial_value_set_per_state(ngp::MultistateField<double> &ngpStatedField)
    {
        ngp::Field<double> ngpField = ngpStatedField.get_field_new_state();
        test_field_has_value(ngpField, 0, -1.0);
        for(unsigned stateCount = 1; stateCount < stkField->number_of_states(); stateCount++)
        {
            ngp::ConstField<double> constNgpField = ngpStatedField.get_field_old_state(static_cast<stk::mesh::FieldState>(stateCount));
            test_field_has_value(constNgpField, stateCount, -1.0);
        }
    }
    void verify_can_assign_values_per_state(ngp::MultistateField<double> &ngpStatedField)
    {
        assign_sequential_values_by_state(ngpStatedField);
        expect_sequential_values_by_state(ngpStatedField);
    }
    template <typename Field>
    void test_field_has_value(Field ngpField, unsigned stateCount, double expectedValue)
    {
        stk::mesh::EntityVector elems;
        stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elems);

        stk::mesh::FieldBase * fieldOfState = stkField->field_state(static_cast<stk::mesh::FieldState>(stateCount));
        ngpField.copy_device_to_host(get_bulk(), *fieldOfState);

        for(stk::mesh::Entity elem : elems)
            EXPECT_EQ(expectedValue, *static_cast<double*>(stk::mesh::field_data(*fieldOfState, elem)));
    }
    stk::mesh::Field<double> * stkField;

private:
    template <typename Field>
    void expect_field_has_value(const Field &ngpField, unsigned stateCount)
    {
        double fieldValue = stateCount;
        test_field_has_value(ngpField, stateCount, fieldValue);
    }

    void assign_sequential_values_by_state(ngp::MultistateField<double>& ngpStatedField)
    {
        assign_value_to_state_new(ngpStatedField, stkField->number_of_states() - 1);
        for(unsigned stateCount = 1; stateCount < stkField->number_of_states(); stateCount++)
        {
            get_bulk().update_field_data_states();
            ngpStatedField.increment_state();
            assign_value_to_state_new(ngpStatedField, stkField->number_of_states() - stateCount - 1);
        }
    }

    void assign_value_to_state_new(ngp::MultistateField<double>& ngpStatedField, double value)
    {
        ngp::Field<double> ngpField = ngpStatedField.get_field_new_state();
        assign_value_to_field(get_bulk(), ngpField, value);
    }

    void expect_sequential_values_by_state(ngp::MultistateField<double>& ngpStatedField)
    {
        ngp::Field<double> ngpField = ngpStatedField.get_field_new_state();
        expect_field_has_value(ngpField, 0);
        for(unsigned stateCount = 1; stateCount < stkField->number_of_states(); stateCount++)
        {
            stk::mesh::FieldState oldState = static_cast<stk::mesh::FieldState>(stateCount);
            ngp::ConstField<double> constNgpField = ngpStatedField.get_field_old_state(oldState);
            expect_field_has_value(constNgpField, stateCount);
        }
    }
};

TEST_F(StatedFields, creatingFromSingleStateStkField_hasOneState)
{
    test_n_states(1);
}
TEST_F(StatedFields, creatingFromDualStateStkField_hasTwoStates)
{
    test_n_states(2);
}
TEST_F(StatedFields, creatingFromSixStateStkField_hasSixStates)
{
    test_n_states(6);
}
TEST_F(StatedFields, creatingFromStatedStkField_ngpStatedFieldHasSameValues)
{
    ngp::MultistateField<double> ngpStatedField = create_field_with_num_states(3);
    verify_initial_value_set_per_state(ngpStatedField);
    verify_can_assign_values_per_state(ngpStatedField);
}
TEST_F(StatedFields, incrementingState_fieldsShiftDown)
{
    ngp::MultistateField<double> ngpStatedField = create_field_with_num_states(3);
    verify_can_assign_values_per_state(ngpStatedField);

    get_bulk().update_field_data_states();
    ngpStatedField.increment_state();

    ngp::Field<double> ngpField = ngpStatedField.get_field_new_state();
    test_field_has_value(ngpField, 0, stkField->number_of_states()-1);
    for(unsigned stateCount = 1; stateCount < stkField->number_of_states(); stateCount++)
    {
        ngp::ConstField<double> constNgpField = ngpStatedField.get_field_old_state(static_cast<stk::mesh::FieldState>(stateCount));
        test_field_has_value(constNgpField, stateCount, stateCount-1);
    }
}
TEST_F(StatedFields, gettingFieldData_direct)
{
    ngp::ConvenientMultistateField<double> ngpStatedField = create_field_with_num_states(3);
    verify_initial_value_set_per_state(ngpStatedField);

    assign_value_to_statedfield(get_bulk(), ngpStatedField, 2);
    get_bulk().update_field_data_states();
    ngpStatedField.increment_state();
    assign_value_to_statedfield(get_bulk(), ngpStatedField, 1);
    get_bulk().update_field_data_states();
    ngpStatedField.increment_state();
    assign_value_to_statedfield(get_bulk(), ngpStatedField, 0);

    ngp::Field<double> ngpField = ngpStatedField.get_field_new_state();
    test_field_has_value(ngpField, 0, 0);
    for(unsigned stateCount = 1; stateCount < stkField->number_of_states(); stateCount++)
    {
        ngp::ConstField<double> constNgpField = ngpStatedField.get_field_old_state(static_cast<stk::mesh::FieldState>(stateCount));
        test_field_has_value(constNgpField, stateCount, stateCount);
    }
}
TEST_F(StatedFields, updateFromHostStkField)
{
    ngp::ConvenientMultistateField<double> ngpStatedField = create_field_with_num_states(1);
    verify_initial_value_set_per_state(ngpStatedField);
    ngp::Field<double> ngpField = ngpStatedField.get_field_new_state();

    const double new_value = 999;
    stk::mesh::field_fill(new_value, *stkField);
    ngpField.copy_host_to_device(*bulkData, *stkField);
    test_field_has_value(ngpField, 0, new_value);
}
} //namespace
