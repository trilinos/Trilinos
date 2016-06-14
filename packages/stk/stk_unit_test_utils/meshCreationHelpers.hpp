#ifndef meshcreationhelpers_hpp
#define meshcreationhelpers_hpp

#include <string>
#include <stk_util/parallel/Parallel.hpp>

namespace stk
{
namespace unit_test_util
{
size_t create_mesh_without_time_steps(const std::string & filename, MPI_Comm communicator, stk::mesh::BulkData & bulk, stk::io::StkMeshIoBroker & stkIo);
void create_mesh_with__field_1__field_2__field_3(const std::string & filename, MPI_Comm communicator);
}
}

#endif
