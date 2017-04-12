#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/FieldBase.hpp>

namespace stk {
namespace balance {
namespace internal {

inline void put_entity_data_to_field(const stk::mesh::EntityProcVec& entityDataToBePutIntoField, stk::mesh::FieldBase *field)
{
    for(const stk::mesh::EntityProc& entityProc : entityDataToBePutIntoField)
    {
        double *fieldData = static_cast<double*>(stk::mesh::field_data(*field, entityProc.first));
        *fieldData = static_cast<double>(entityProc.second);
    }
}

}
}
}
