#include "stk_io/FillMesh.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"

namespace stk
{
namespace io
{

void fill_mesh_preexisting(stk::io::StkMeshIoBroker & stkIo, const std::string& meshSpec, stk::mesh::BulkData& bulkData)
{
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(meshSpec, stk::io::READ_MESH);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();
    stkIo.populate_bulk_data();
}

void fill_mesh_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData)
{
    stk::io::StkMeshIoBroker stkIo;
    stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    fill_mesh_preexisting(stkIo, meshSpec, bulkData);
}

void fill_mesh(const std::string &meshSpec, stk::mesh::BulkData &bulkData)
{
    stk::io::StkMeshIoBroker stkIo;
    fill_mesh_preexisting(stkIo, meshSpec, bulkData);
}

}
}
