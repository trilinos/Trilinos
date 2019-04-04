#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/FieldBase.hpp"
#include "stk_topology/topology.hpp"
#include "stk_util/util/ReportHandler.hpp"
#include "MeshCloneUtils.hpp"

namespace stk
{
namespace tools
{

void copy_dist_fact_field(stk::mesh::MetaData &newMeta, const stk::mesh::Part &oldPart, stk::mesh::Part &newPart)
{
    const stk::mesh::FieldBase *oldDistFactField = stk::io::get_distribution_factor_field(oldPart);
    if(oldDistFactField != nullptr)
    {
        unsigned fieldOrd = oldDistFactField->mesh_meta_data_ordinal();
        const stk::mesh::FieldBase *newDistFactField = newMeta.get_fields()[fieldOrd];
        stk::io::set_distribution_factor_field(newPart, *newDistFactField);
    }
}

void copy_io_part_attributes(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::PartVector &oldParts = oldMeta.get_parts();
    for(size_t i = 0; i < oldParts.size(); i++)
    {
        stk::mesh::Part* newPart = get_corresponding_part(newMeta, oldParts[i]);
        if(newPart != nullptr)
        {
            if(stk::io::is_part_io_part(*oldParts[i]))
                stk::io::put_io_part_attribute(*newPart);
            copy_dist_fact_field(newMeta, *oldParts[i], *newPart);
        }
    }
}

void copy_field_roles(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    const stk::mesh::FieldVector &oldFields = oldMeta.get_fields();
    const stk::mesh::FieldVector &newFields = newMeta.get_fields();
    ThrowRequire(oldFields.size() == newFields.size());
    for(size_t i = 0; i < oldFields.size(); i++)
    {
        const Ioss::Field::RoleType* role = stk::io::get_field_role(*oldFields[i]);
        if(role != nullptr)
            stk::io::set_field_role(*newFields[i], *role);
    }
}

void copy_io_attributes(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
    copy_io_part_attributes(oldMeta, newMeta);
    copy_field_roles(oldMeta, newMeta);
}
}
}
