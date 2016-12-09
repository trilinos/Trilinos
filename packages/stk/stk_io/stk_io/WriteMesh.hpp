#ifndef PACKAGES_STK_STK_IO_STK_IO_WRITEMESH_HPP_
#define PACKAGES_STK_STK_IO_STK_IO_WRITEMESH_HPP_

#include <string>
#include <stk_io/DatabasePurpose.hpp> // for stk::io::DatabasePurpose
namespace stk { namespace mesh { class BulkData; }}

namespace stk
{
namespace io
{

void write_mesh(const std::string &filename,
                stk::mesh::BulkData &bulkData,
                stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

} // namespace unit_test_util
} // namespace stk

#endif /* PACKAGES_STK_STK_IO_STK_IO_WRITEMESH_HPP_ */
