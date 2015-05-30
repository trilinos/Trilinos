#include "ioUtils.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"   // for BulkData

namespace stk
{
namespace unit_test_util
{

void fill_mesh_using_stk_io_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::ParallelMachine communicator)
{
    stk::io::StkMeshIoBroker exodusFileReader(communicator);
    exodusFileReader.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));

    exodusFileReader.set_bulk_data(bulkData);
    exodusFileReader.add_mesh_database(meshSpec, stk::io::READ_MESH);
    exodusFileReader.create_input_mesh();
    exodusFileReader.populate_bulk_data();
}

void fill_mesh_using_stk_io(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::ParallelMachine communicator)
{
    stk::io::StkMeshIoBroker exodusFileReader(communicator);

    exodusFileReader.set_bulk_data(bulkData);
    exodusFileReader.add_mesh_database(meshSpec, stk::io::READ_MESH);
    exodusFileReader.create_input_mesh();
    exodusFileReader.populate_bulk_data();
}

void write_mesh_using_stk_io(const std::string &filename,
                             stk::mesh::BulkData &bulkData,
                             stk::ParallelMachine communicator,
                             stk::io::DatabasePurpose databasePurpose)
{
    stk::io::StkMeshIoBroker exodusFileWriter(communicator);

    exodusFileWriter.set_bulk_data(bulkData);
    size_t output_file_index = exodusFileWriter.create_output_mesh(filename, databasePurpose);
    exodusFileWriter.write_output_mesh(output_file_index);
}

} // namespace unit_test_util
} // namespace stk

