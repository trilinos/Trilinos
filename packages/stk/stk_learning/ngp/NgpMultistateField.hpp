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
            numStates(stkField.number_of_states())
    {
        for(unsigned i = 0; i < numStates; i++)
        {
            stk::mesh::FieldBase * fieldOfState = stkField.field_state(static_cast<stk::mesh::FieldState>(i));
            fields[i] = ngp::Field<T>(bulk, *fieldOfState);
        }
    }
    virtual ~MultistateField() {}
    STK_FUNCTION
    unsigned get_num_states() const
    {
        return numStates;
    }
    STK_FUNCTION
    ngp::Field<T> get_field_of_state(stk::mesh::FieldState state) const
    {
        return fields[state];
    }
    STK_FUNCTION
    void increment_state()
    {
        ngp::Field<T> oldLast = fields[numStates-1];
        for(unsigned i=numStates-1; i>0; i--)
            fields[i] = fields[i-1];
        fields[0] = oldLast;
    }

    void copy_device_to_host(const stk::mesh::BulkData& bulk, stk::mesh::FieldBase &field)
    {
        for(unsigned i=0; i<numStates; i++)
            fields[i].copy_device_to_host(bulk, *field.field_state(static_cast<stk::mesh::FieldState>(i)));
    }

private:
    unsigned numStates;
    ngp::Field<T> fields[MAX_NUM_FIELD_STATES];
};

template<typename T>
class ConvenientMultistateField : public MultistateField<T>
{
public:
    ConvenientMultistateField(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField) :
            MultistateField<T>(bulk, stkField)
    {}

    STK_FUNCTION
    T& get(stk::mesh::FieldState state, const ngp::Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        return MultistateField<T>::get_field_of_state(state).get(ngpMesh, entity, component);
    }

    STK_FUNCTION
    const T& const_get(stk::mesh::FieldState state, const ngp::Mesh& ngpMesh, stk::mesh::Entity entity, int component) const
    {
        return MultistateField<T>::get_field_of_state(state).const_get(ngpMesh, entity, component);
    }

    STK_FUNCTION
    T& get(stk::mesh::FieldState state, ngp::Mesh::MeshIndex entity, int component) const
    {
        return MultistateField<T>::get_field_of_state(state).get(entity, component);
    }

    STK_FUNCTION
    const T& const_get(stk::mesh::FieldState state, ngp::Mesh::MeshIndex entity, int component) const
    {
        return MultistateField<T>::get_field_of_state(state).const_get(entity, component);
    }
};

}

#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPSTATEDFIELD_H_ */
