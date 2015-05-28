#ifndef ioUtils_hpp
#define ioUtils_hpp

#include <string>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_io/DatabasePurpose.hpp> // for stk::io::DatabasePurpose
namespace stk { namespace mesh { class BulkData; }}

namespace stk
{
namespace unit_test_util
{

void fill_mesh_using_stk_io(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::ParallelMachine communicator);
void fill_mesh_using_stk_io_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::ParallelMachine communicator);

void write_mesh_using_stk_io(const std::string &filename,
                             stk::mesh::BulkData &bulkData,
                             stk::ParallelMachine communicator,
                             stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

} // namespace unit_test_util
} // namespace stk


#endif
