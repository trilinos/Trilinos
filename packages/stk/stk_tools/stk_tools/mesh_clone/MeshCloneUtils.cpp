#include "MeshCloneUtils.hpp"
#include <stk_mesh/base/MetaData.hpp>
#include "stk_mesh/base/Part.hpp"

namespace stk {
namespace tools {

stk::mesh::Part * get_corresponding_part(stk::mesh::MetaData &newMeta, const stk::mesh::Part *oldPart)
{
  return newMeta.get_part(oldPart->name());
}

stk::mesh::OrdinalVector get_part_supersets(const stk::mesh::Part &part)
{
  const stk::mesh::PartVector &supersetParts = part.supersets();
  stk::mesh::OrdinalVector supersetOrdinals(supersetParts.size());
  for(size_t i=0; i<supersetParts.size(); i++)
    supersetOrdinals[i] = supersetParts[i]->mesh_meta_data_ordinal();
  return supersetOrdinals;
}

void copy_part_supersets(const stk::mesh::Part &oldPart, stk::mesh::Part &newPart, const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta)
{
  stk::mesh::OrdinalVector oldSupersets = get_part_supersets(oldPart);
  for(stk::mesh::PartOrdinal partOrd : oldSupersets)
  {
    const std::string &oldName = oldMeta.get_part(partOrd).name();
    newMeta.declare_part_subset(*newMeta.get_part(oldName), newPart);
  }
}

}
}
