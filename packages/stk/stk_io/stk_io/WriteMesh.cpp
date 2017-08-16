#include "stk_io/WriteMesh.hpp"
#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"

namespace stk
{
namespace io
{

void write_mesh(const std::string &filename,
                stk::mesh::BulkData &bulkData,
                stk::io::DatabasePurpose databasePurpose)
{
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(bulkData);
    size_t outputFileIndex = stkIo.create_output_mesh(filename, databasePurpose);
    stkIo.write_output_mesh(outputFileIndex);
}

void write_mesh_subset(const std::string &filename,
                stk::mesh::BulkData &bulkData,
                stk::mesh::Selector& subsetSelector,
                stk::io::DatabasePurpose databasePurpose)
{
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(bulkData);
    size_t outputFileIndex = stkIo.create_output_mesh(filename, databasePurpose);
    stkIo.set_subset_selector(outputFileIndex, subsetSelector);
    stkIo.write_output_mesh(outputFileIndex);
}

void write_mesh_with_fields(const std::string& filename, stk::mesh::BulkData &bulkData, int step, double time)
{
    stk::io::StkMeshIoBroker outStkIo;
    outStkIo.set_bulk_data(bulkData);

    size_t outputFileIndex = outStkIo.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    const stk::mesh::FieldVector fields = bulkData.mesh_meta_data().get_fields();
    for(stk::mesh::FieldBase* field : fields)
    {
        const Ioss::Field::RoleType* fieldRole = stk::io::get_field_role(*field);
        if(fieldRole == nullptr || *fieldRole == Ioss::Field::TRANSIENT)
            outStkIo.add_field(outputFileIndex, *field);
    }

    outStkIo.write_output_mesh(outputFileIndex);
    if(step>0)
    {
        outStkIo.begin_output_step(outputFileIndex, time);
        outStkIo.write_defined_output_fields(outputFileIndex);
        outStkIo.end_output_step(outputFileIndex);
    }
}

}
}
