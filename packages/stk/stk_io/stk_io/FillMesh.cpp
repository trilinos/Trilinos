// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_io/FillMesh.hpp"
#include "Ioss_Property.h"             // for Property
#include "stk_io/StkMeshIoBroker.hpp"  // for StkMeshIoBroker
#include "stk_mesh/base/BulkData.hpp"  // for BulkData
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace io
{

void fill_mesh_preexisting(stk::io::StkMeshIoBroker & stkIo,
                           const std::string& meshSpec,
                           stk::mesh::BulkData& bulkData,
                           stk::io::DatabasePurpose purpose)
{
    stkIo.set_bulk_data(bulkData);
    stkIo.add_mesh_database(meshSpec, purpose);
    stkIo.create_input_mesh();
    stkIo.add_all_mesh_fields_as_input_fields();
    const bool delayFieldDataAllocation = true;
    stkIo.populate_mesh(delayFieldDataAllocation);
    stkIo.populate_field_data();
}

void fill_mesh_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::io::StkMeshIoBroker &stkIo)
{
    stkIo.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    fill_mesh_preexisting(stkIo, meshSpec, bulkData);
}

void fill_mesh_with_auto_decomp(const std::string &meshSpec, stk::mesh::BulkData &bulkData)
{
    stk::io::StkMeshIoBroker stkIo;
    fill_mesh_with_auto_decomp(meshSpec, bulkData, stkIo);
}

void fill_mesh(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::io::StkMeshIoBroker &stkIo)
{
    fill_mesh_preexisting(stkIo, meshSpec, bulkData);
}

void fill_mesh(const std::string &meshSpec, stk::mesh::BulkData &bulkData)
{
    stk::io::StkMeshIoBroker stkIo;
    fill_mesh_preexisting(stkIo, meshSpec, bulkData);
}

void save_step_info(stk::io::StkMeshIoBroker &stkIo, int &numSteps, double &maxTime)
{
    numSteps = stkIo.get_num_time_steps();
    if(numSteps>0)
    {
        stkIo.read_defined_input_fields(numSteps);
        maxTime = stkIo.get_max_time();
    }
}

void fill_mesh_save_step_info(const std::string& inFile, stk::mesh::BulkData& inBulk, int &numSteps, double &maxTime)
{
    stk::io::StkMeshIoBroker stkIo;
    stk::io::fill_mesh_preexisting(stkIo, inFile, inBulk);
    save_step_info(stkIo, numSteps, maxTime);
}

void fill_mesh_with_fields(const std::string& inFile, stk::mesh::BulkData& bulk, stk::io::DatabasePurpose purpose)
{
    stk::io::StkMeshIoBroker stkIo;
    fill_mesh_with_fields(inFile, stkIo, bulk, purpose);
}

void fill_mesh_with_fields(const std::string& inFile, stk::io::StkMeshIoBroker& ioBroker,
                           stk::mesh::BulkData& bulk, stk::io::DatabasePurpose purpose)
{
    int numSteps;
    double maxTime;
    stk::io::fill_mesh_preexisting(ioBroker, inFile, bulk, purpose);
    save_step_info(ioBroker, numSteps, maxTime);
}

}
}
