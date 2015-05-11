#include "ioUtils.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData

namespace stk
{
namespace unit_test_util
{

void fill_mesh_using_stk_io_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData, MPI_Comm communicator)
{
    stk::io::StkMeshIoBroker exodusFileReader(communicator);
    exodusFileReader.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));

    exodusFileReader.set_bulk_data(bulkData);
    exodusFileReader.add_mesh_database(meshSpec, stk::io::READ_MESH);
    exodusFileReader.create_input_mesh();
    exodusFileReader.populate_bulk_data();
}

void fill_mesh_using_stk_io(const std::string &meshSpec, stk::mesh::BulkData &bulkData, MPI_Comm communicator)
{
    stk::io::StkMeshIoBroker exodusFileReader(communicator);

    exodusFileReader.set_bulk_data(bulkData);
    exodusFileReader.add_mesh_database(meshSpec, stk::io::READ_MESH);
    exodusFileReader.create_input_mesh();
    exodusFileReader.populate_bulk_data();
}

} // namespace unit_test_util
} // namespace stk

