#include <stk_io/AttributeFields.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "IossBridge.hpp"

namespace stk {
namespace io {

stk::mesh::FieldVector get_attribute_fields_for_part(const stk::mesh::MetaData &meta, const stk::mesh::Part *ioPart)
{
    stk::mesh::FieldVector attributes;

    for(stk::mesh::FieldBase *field : meta.get_fields())
    {
        const Ioss::Field::RoleType *fieldRole = stk::io::get_field_role(*field);
        if(fieldRole != nullptr && *fieldRole == Ioss::Field::ATTRIBUTE)
        {
            for(const stk::mesh::FieldBase::Restriction &restriction : field->restrictions())
            {
                const stk::mesh::Selector &selector = restriction.selector();
                if(selector(ioPart))
                    attributes.push_back(field);
            }
        }
    }
    return attributes;
}

void mark_field_as_attribute(stk::mesh::FieldBase &field)
{
    stk::io::set_field_role(field, Ioss::Field::ATTRIBUTE);
}

} // namespace io
} // namespace stk
