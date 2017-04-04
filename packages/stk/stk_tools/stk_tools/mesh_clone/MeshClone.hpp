#ifndef PACKAGES_STK_STK_MESH_STK_MESH_BASE_MESHCLONE_HPP_
#define PACKAGES_STK_STK_MESH_STK_MESH_BASE_MESHCLONE_HPP_

#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Selector; } }

namespace stk {
namespace tools {

void copy_meta_with_io_attributes(const stk::mesh::MetaData &inputMeta, stk::mesh::MetaData &outputMeta);
void copy_mesh(const stk::mesh::BulkData &inputBulk, stk::mesh::Selector selector, stk::mesh::BulkData &outputBulk);
void copy_bulk(const stk::mesh::BulkData &inputBulk, stk::mesh::Selector selector, stk::mesh::BulkData &outputBulk);

}
}

#endif /* PACKAGES_STK_STK_MESH_STK_MESH_BASE_MESHCLONE_HPP_ */
