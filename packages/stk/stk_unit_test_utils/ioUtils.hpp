#ifndef ioUtils_hpp
#define ioUtils_hpp

#include <string>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc

namespace stk { namespace mesh { class BulkData; }}

namespace stk
{
namespace unit_test_util
{

void fill_mesh_using_stk_io(const std::string &meshSpec, stk::mesh::BulkData &bulkData, MPI_Comm communicator);

} // namespace unit_test_util
} // namespace stk


#endif
