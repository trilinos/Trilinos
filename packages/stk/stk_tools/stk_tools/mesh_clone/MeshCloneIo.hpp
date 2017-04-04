#ifndef PACKAGES_STK_STK_IO_STK_IO_MESHCLONEIO_HPP_
#define PACKAGES_STK_STK_IO_STK_IO_MESHCLONEIO_HPP_

namespace stk { namespace mesh { class MetaData; } }

namespace stk {
namespace tools {

void copy_io_attributes(const stk::mesh::MetaData &oldMeta, stk::mesh::MetaData &newMeta);

}}

#endif /* PACKAGES_STK_STK_IO_STK_IO_MESHCLONEIO_HPP_ */
