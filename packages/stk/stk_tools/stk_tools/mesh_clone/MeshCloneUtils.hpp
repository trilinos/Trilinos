#ifndef PACKAGES_STK_STK_MESH_STK_MESH_BASE_MESHCLONE_UTILS_HPP_
#define PACKAGES_STK_STK_MESH_STK_MESH_BASE_MESHCLONE_UTILS_HPP_

#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class MetaData; } }

namespace stk {
namespace tools {

stk::mesh::OrdinalVector get_part_supersets(const stk::mesh::Part &part);
void copy_part_supersets(const stk::mesh::Part &oldPart, stk::mesh::Part &newPart, const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta);
stk::mesh::Part * get_corresponding_part(stk::mesh::MetaData &newMeta, const stk::mesh::Part *oldPart);

}
}


#endif
