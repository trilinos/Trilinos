#include "StkMeshUnitTestUtils.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData

namespace stk
{
namespace mesh
{
namespace unit_test
{

void readMesh(const std::string &meshSpec, stk::mesh::BulkData &bulkData, MPI_Comm communicator)
{
    stk::io::StkMeshIoBroker exodusFileReader(communicator);
    exodusFileReader.set_bulk_data(bulkData);
    exodusFileReader.add_mesh_database(meshSpec, stk::io::READ_MESH);
    exodusFileReader.create_input_mesh();
    exodusFileReader.populate_bulk_data();
}

} // namespace unit_test
} // namespace mesh
} // namespace stk

