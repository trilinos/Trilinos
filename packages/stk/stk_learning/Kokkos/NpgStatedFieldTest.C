#include <ngp/Ngp.hpp>
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


namespace ngp
{

const unsigned MAX_NUM_STATES = 6;

class NgpStatedField
{
public:
    NgpStatedField(const stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField) :
            numStates(stkField.number_of_states())
    {
        for(unsigned i = 0; i < numStates; i++)
        {
            stk::mesh::FieldBase * fieldOfState = stkField.field_state(static_cast<stk::mesh::FieldState>(i));
            fields[i] = ngp::StkNgpField(bulk, *fieldOfState);
        }
    }
    unsigned get_num_states()
    {
        return numStates;
    }
    ngp::StkNgpField get_field_of_state(stk::mesh::FieldState fs)
    {
        return fields[fs];
    }
private:
    unsigned numStates;
    ngp::StkNgpField fields[MAX_NUM_STATES];
};

}

namespace {


class StatedFields : public stk::unit_test_util::MeshFixture
{
protected:
    ngp::NgpStatedField create_field_with_num_states(unsigned numStates)
    {
        setup_empty_mesh(stk::mesh::BulkData::NO_AUTO_AURA);
        stkField = &get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEM_RANK, "myField", numStates);
        double init = -1.0;
        stk::mesh::put_field(*stkField, get_meta().universal_part(), &init);
        stk::unit_test_util::fill_mesh_using_stk_io("generated:1x1x4", get_bulk());

//        stk::mesh::EntityVector elems;
//        stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elems);
//        for(unsigned stateCount = 0; stateCount < numStates; stateCount++)
//        {
//            stk::mesh::FieldBase * fieldOfState = stkField->field_state(static_cast<stk::mesh::FieldState>(stateCount));
//            for(stk::mesh::Entity elem : elems)
//            {
//                double * value = static_cast<double*>(stk::mesh::field_data(*fieldOfState, elem));
//                *value = stateCount;
//            }
//        }

        return ngp::NgpStatedField(get_bulk(), *stkField);
    }
    void test_n_states(unsigned numStates)
    {
        ngp::NgpStatedField ngpField = create_field_with_num_states(numStates);
        EXPECT_EQ(numStates, ngpField.get_num_states());
    }
    void test_field_has_value(ngp::StkNgpField ngpField, unsigned stateCount, double expectedValue)
    {
        stk::mesh::EntityVector elems;
        stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, elems);

        stk::mesh::FieldBase * fieldOfState = stkField->field_state(static_cast<stk::mesh::FieldState>(stateCount));
        ngpField.copy_device_to_host(get_bulk(), *fieldOfState);

        for(stk::mesh::Entity elem : elems)
        {
            double * value = static_cast<double*>(stk::mesh::field_data(*fieldOfState, elem));
            EXPECT_EQ(expectedValue, *value);
        }
    }
    void verify_initial_value_set_per_state(ngp::NgpStatedField &ngpStatedField)
    {
        for(unsigned stateCount = 0; stateCount < stkField->number_of_states(); stateCount++)
        {
            ngp::StkNgpField ngpField = ngpStatedField.get_field_of_state(static_cast<stk::mesh::FieldState>(stateCount));
            test_field_has_value(ngpField, stateCount, -1.0);
        }
    }
    void assign_value_to_field(ngp::StkNgpField ngpField, double value)
    {
        ngp::StkNgpMesh ngpMesh(get_bulk());
        ngp::for_each_entity_run(ngpMesh, stk::topology::ELEM_RANK, get_meta().universal_part(), KOKKOS_LAMBDA(ngp::StkNgpMesh::MeshIndex entity)
        {
            ngpField.get(entity, 0) = value;
        });
    }
    void verify_can_assign_values_per_state(ngp::NgpStatedField &ngpStatedField)
    {
        for(unsigned stateCount = 0; stateCount < stkField->number_of_states(); stateCount++)
        {
            ngp::StkNgpField ngpField = ngpStatedField.get_field_of_state(static_cast<stk::mesh::FieldState>(stateCount));
            double fieldValue = stateCount;
            assign_value_to_field(ngpField, fieldValue);
            test_field_has_value(ngpField, stateCount, fieldValue);
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
    ngp::NgpStatedField ngpStatedField = create_field_with_num_states(3);
    verify_initial_value_set_per_state(ngpStatedField);
    verify_can_assign_values_per_state(ngpStatedField);
}
} //namespace
