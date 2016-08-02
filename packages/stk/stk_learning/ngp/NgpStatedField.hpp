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

class NgpStatedField
{
public:
    NgpStatedField(stk::mesh::BulkData& bulk, const stk::mesh::FieldBase &stkField) :
            numStates(stkField.number_of_states())
    {
        for(unsigned i = 0; i < numStates; i++)
        {
            stk::mesh::FieldBase * fieldOfState = stkField.field_state(static_cast<stk::mesh::FieldState>(i));
            fields[i] = ngp::StkNgpField(bulk, *fieldOfState);
        }
    }
    STK_FUNCTION
    unsigned get_num_states()
    {
        return numStates;
    }
    STK_FUNCTION
    ngp::StkNgpField get_field_of_state(stk::mesh::FieldState fs)
    {
        return fields[fs];
    }
    STK_FUNCTION
    void increment_state()
    {
        ngp::StkNgpField oldLast = fields[numStates-1];
        for(unsigned i=numStates-1; i>0; i--)
            fields[i] = fields[i-1];
        fields[0] = oldLast;
    }
private:
    unsigned numStates;
    ngp::StkNgpField fields[MAX_NUM_FIELD_STATES];
};

}

#endif /* PACKAGES_STK_STK_LEARNING_KOKKOS_NGPSTATEDFIELD_H_ */
