#include <ngp/Ngp.hpp>
#include <ngp/NgpStatedField.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/stk_config.h>


void assign_value_to_field(stk::mesh::BulkData &bulk, ngp::StkNgpField<double> ngpField, double value)
{
    ngp::StkNgpMesh ngpMesh(bulk);
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part(), KOKKOS_LAMBDA(ngp::StkNgpMesh::MeshIndex entity)
    {
        ngpField.get(entity, 0) = value;
    });
}
void assign_value_to_statedfield(stk::mesh::BulkData &bulk, ngp::NgpStatedField<double> ngpStatedField, unsigned state, double value)
{
    ngp::StkNgpMesh ngpMesh(bulk);
    ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, bulk.mesh_meta_data().universal_part(), KOKKOS_LAMBDA(ngp::StkNgpMesh::MeshIndex entity)
    {
        ngpStatedField.get(static_cast<stk::mesh::FieldState>(state), entity, 0) = value;
    });
}

namespace {

class StatedFields : public stk::unit_test_util::MeshFixture
{
protected:
    ngp::NgpStatedField<double> create_field_with_num_states(unsigned numStates)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        stkField = &get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "myField", numStates);
        double init = -1.0;
        stk::mesh::put_field(*stkField, get_meta().universal_part(), &init);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", get_bulk());

        return ngp::NgpStatedField<double>(get_bulk(), *stkField);
    }
    void test_n_states(unsigned numStates)
    {
        ngp::NgpStatedField<double> ngpField = create_field_with_num_states(numStates);
        EXPECT_EQ(numStates, ngpField.get_num_states());
    }
    void verify_initial_value_set_per_state(ngp::NgpStatedField<double> &ngpStatedField)
    {
        for(unsigned stateCount = 0; stateCount < stkField->number_of_states(); stateCount++)
            test_field_has_value(ngpStatedField, stateCount, -1.0);
    }
    void verify_can_assign_values_per_state(ngp::NgpStatedField<double> &ngpStatedField)
    {
        for(unsigned stateCount = 0; stateCount < stkField->number_of_states(); stateCount++)
        {
            ngp::StkNgpField<double> ngpField = ngpStatedField.get_field_of_state(static_cast<stk::mesh::FieldState>(stateCount));
            double fieldValue = stateCount;
            assign_value_to_field(get_bulk(), ngpField, fieldValue);
            test_field_has_value(ngpStatedField, stateCount, fieldValue);
        }
    }
    void test_field_has_value(ngp::NgpStatedField<double> ngpStatedField, unsigned stateCount, double expectedValue)
    {
        stk::mesh::EntityVector elems;
        stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elems);

        ngp::StkNgpField<double> ngpField = ngpStatedField.get_field_of_state(static_cast<stk::mesh::FieldState>(stateCount));
        stk::mesh::FieldBase * fieldOfState = stkField->field_state(static_cast<stk::mesh::FieldState>(stateCount));
        ngpField.copy_device_to_host(get_bulk(), *fieldOfState);

        for(stk::mesh::Entity elem : elems)
        {
            double * value = static_cast<double*>(stk::mesh::field_data(*fieldOfState, elem));
            EXPECT_EQ(expectedValue, *value);
        }
    }
    stk::mesh::Field<double> * stkField;
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
    ngp::NgpStatedField<double> ngpStatedField = create_field_with_num_states(3);
    verify_initial_value_set_per_state(ngpStatedField);
    verify_can_assign_values_per_state(ngpStatedField);
}
TEST_F(StatedFields, incrementingState_fieldsShiftDown)
{
    ngp::NgpStatedField<double> ngpStatedField = create_field_with_num_states(3);
    verify_can_assign_values_per_state(ngpStatedField);

    get_bulk().update_field_data_states();
    ngpStatedField.increment_state();

    test_field_has_value(ngpStatedField, 0, stkField->number_of_states()-1);
    for(unsigned stateCount = 1; stateCount < stkField->number_of_states(); stateCount++)
    {
        test_field_has_value(ngpStatedField, stateCount, stateCount-1);
    }
}
TEST_F(StatedFields, gettingFieldData_direct)
{
    ngp::NgpStatedField<double> ngpStatedField = create_field_with_num_states(3);
    verify_initial_value_set_per_state(ngpStatedField);
    for(unsigned stateCount = 0; stateCount < stkField->number_of_states(); stateCount++)
    {
        double fieldValue = stateCount;
        assign_value_to_statedfield(get_bulk(), ngpStatedField, stateCount, fieldValue);
        test_field_has_value(ngpStatedField, stateCount, fieldValue);
    }
}
} //namespace
