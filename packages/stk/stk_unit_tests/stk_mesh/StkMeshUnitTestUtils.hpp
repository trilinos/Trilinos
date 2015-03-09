#ifndef stkMeshUnitTestUtils_hpp
#define stkMeshUnitTestUtils_hpp

#include <string>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc

namespace stk
{
namespace mesh
{

class BulkData;

namespace unit_test
{

void readMesh(const std::string &meshSpec, stk::mesh::BulkData &bulkData, MPI_Comm communicator);

} // namespace unit_test
} // namespace mesh
} // namespace stk


#endif
