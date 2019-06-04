// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_io/AttributeFields.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "IossBridge.hpp"               // for get_field_role, etc
#include "Ioss_Field.h"                 // for Field, Field::RoleType, etc
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, etc
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for FieldVector
#include "stk_topology/topology.hpp"    // for topology, etc
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
namespace io {

Ioss::Field::RoleType get_field_role_for_non_outputted_attribute()
{
    return Ioss::Field::INTERNAL;
}

stk::mesh::FieldVector get_internal_fields_for_part(const stk::mesh::MetaData &meta, const stk::mesh::Part *ioPart)
{
    stk::mesh::FieldVector attributes;

    for(stk::mesh::FieldBase *field : meta.get_fields())
    {
        const Ioss::Field::RoleType *fieldRole = stk::io::get_field_role(*field);
        if(fieldRole != nullptr && *fieldRole == get_field_role_for_non_outputted_attribute())
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

stk::mesh::FieldBase* get_named_internal_field_for_part(const stk::mesh::MetaData &meta, const stk::mesh::Part *ioPart, const std::string &fieldName)
{
    stk::mesh::FieldBase *field = meta.get_field(stk::topology::ELEM_RANK, fieldName);

    if (field != nullptr) {
    	const Ioss::Field::RoleType *fieldRole = stk::io::get_field_role(*field);
    	if(fieldRole != nullptr && *fieldRole == get_field_role_for_non_outputted_attribute())
    	{
    		for(const stk::mesh::FieldBase::Restriction &restriction : field->restrictions())
    		{
    			const stk::mesh::Selector &selector = restriction.selector();
    			if(selector(ioPart))
    				return field;
    		}
    	}
    }

    return nullptr;
}

void mark_field_as_internal(stk::mesh::FieldBase &field)
{
    stk::io::set_field_role(field, get_field_role_for_non_outputted_attribute());
}

} // namespace io
} // namespace stk
