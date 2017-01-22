#ifndef PACKAGES_STK_STK_IO_STK_IO_FILLMESH_HPP_
#define PACKAGES_STK_STK_IO_STK_IO_FILLMESH_HPP_

#include <string>
namespace stk { namespace mesh { class BulkData; }}
namespace stk { namespace io { class StkMeshIoBroker; } }

namespace stk
{
namespace io
{

void fill_mesh(const std::string &meshSpec, stk::mesh::BulkData &bulkData);
void fill_mesh_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData);
void fill_mesh_preexisting(stk::io::StkMeshIoBroker & stkIo, const std::string& meshSpec, stk::mesh::BulkData& bulkData);

} // namespace unit_test_util
} // namespace stk

#endif /* PACKAGES_STK_STK_IO_STK_IO_FILLMESH_HPP_ */
