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

void fill_mesh_using_stk_io(const std::string &meshSpec, stk::mesh::BulkData &bulkData);
void fill_mesh_using_stk_io_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData);

void write_mesh_using_stk_io(const std::string &filename,
                             stk::mesh::BulkData &bulkData,
                             stk::io::DatabasePurpose databasePurpose = stk::io::WRITE_RESULTS);

// Example:  meshSizeSpec = "2x2x1"
void generated_mesh_to_file_in_serial(const std::string &meshSizeSpec, const std::string &fileName);

void read_from_serial_file_and_decompose(const std::string& fileName, stk::mesh::BulkData &mesh, const std::string &decompositionMethod);

// This avoid the GeneratedMesh limitation on the z-dimension >= num_processors,
// and allows cases like 2x2x1 with cyclic decomposition on four processors.
void generate_mesh_from_serial_spec_and_load_in_parallel_with_auto_decomp(const std::string &meshSizeSpec, stk::mesh::BulkData & mesh, const std::string &decompositionMethod);

} // namespace unit_test_util
} // namespace stk


#endif
