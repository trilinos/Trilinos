#ifndef PACKAGES_STK_STK_LEARNING_KOKKOS_NGPSTATEDFIELD_H_
#define PACKAGES_STK_STK_LEARNING_KOKKOS_NGPSTATEDFIELD_H_

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

const unsigned MAX_NUM_FIELD_STATES = 6;

template<typename T>
class MultistateField
{
public:
    MultistateField() {}
    MultistateField(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField) :
        numOldStates(stkField.number_of_states() - 1)
    {
        stk::mesh::FieldBase * fieldOfState = stkField.field_state(stk::mesh::StateNew);
        fieldNewState = ngp::Field<T>(bulk, *fieldOfState);
        for(unsigned i = 0; i < numOldStates; i++)
        {
            fieldOfState = stkField.field_state(get_old_state_for_index(i));
            fieldOldStates[i] = ngp::ConstField<T>(bulk, *fieldOfState);
        }
    }
    STK_FUNCTION
    virtual ~MultistateField() {}
    STK_FUNCTION
    unsigned get_num_states() const
    {
        return numOldStates + 1;
    }
    STK_FUNCTION
    ngp::Field<T> get_field_new_state() const
    {
        return fieldNewState;
    }
    STK_FUNCTION
    ngp::ConstField<T> get_field_old_state(stk::mesh::FieldState state) const
    {
        ThrowRequire(state != stk::mesh::StateNew);
        return fieldOldStates[state-1];
    }
    STK_FUNCTION
    void increment_state()
    {
        for(unsigned i=numOldStates-1; i>0; i--)
            fieldOldStates[i].swap_data(fieldOldStates[i-1]);
        fieldOldStates[0].swap_data(fieldNewState);
    }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field)
    {
        fieldNewState.copy_device_to_host(bulk, *field.field_state(stk::mesh::StateNew));
        for(unsigned i=0; i<numOldStates; i++)
            fieldOldStates[i].copy_device_to_host(bulk, *field.field_state(get_old_state_for_index(i)));
    }

private:
    stk::mesh::FieldState get_old_state_for_index(unsigned i)
    {
        return static_cast<stk::mesh::FieldState>(i+1);
    }
    unsigned numOldStates;
    ngp::Field<T> fieldNewState;
    ngp::ConstField<T> fieldOldStates[MAX_NUM_FIELD_STATES];
};

template<typename T>
class ConvenientMultistateField : public MultistateField<T>
{
public:
    ConvenientMultistateField(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField) :
            MultistateField<T>(bulk, stkField)
    {}

    STK_FUNCTION
    T& get_new(const ngp::Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        return MultistateField<T>::get_field_new_state().get(ngpMesh, entity, component);
    }

    STK_FUNCTION
    const T& get_old(stk::mesh::FieldState state, const ngp::Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        return MultistateField<T>::get_field_old_state(state).get(ngpMesh, entity, component);
    }

    STK_FUNCTION
    T& get_new(ngp::Mesh::MeshIndex entity, int component) const
    {
        return MultistateField<T>::get_field_new_state().get(entity, component);
    }

    STK_FUNCTION
    const T& get_old(stk::mesh::FieldState state, ngp::Mesh::MeshIndex entity, int component) const
    {
        return MultistateField<T>::get_field_old_state(state).get(entity, component);
    }
};

}

#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPSTATEDFIELD_H_ */
