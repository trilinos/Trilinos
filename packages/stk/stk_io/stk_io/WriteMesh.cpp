// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "stk_io/WriteMesh.hpp"
#include <stddef.h>                    // for size_t
#include "Ioss_Field.h"                // for Field, Field::RoleType, etc
#include "stk_io/DatabasePurpose.hpp"  // for DatabasePurpose, etc
#include "stk_io/IossBridge.hpp"       // for get_field_role
#include "stk_io/StkMeshIoBroker.hpp"  // for StkMeshIoBroker
#include "stk_mesh/base/BulkData.hpp"  // for BulkData
#include "stk_mesh/base/MetaData.hpp"  // for MetaData
#include "stk_mesh/base/Types.hpp"     // for FieldVector
namespace stk { namespace mesh { class FieldBase; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

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

void write_mesh_with_large_ids(const std::string &filename,
                stk::mesh::BulkData &bulkData,
                stk::io::DatabasePurpose databasePurpose)
{
    stk::io::StkMeshIoBroker stkIo;
    stkIo.set_bulk_data(bulkData);
    stkIo.property_add(Ioss::Property("INTEGER_SIZE_API" , 8));
    stkIo.property_add(Ioss::Property("INTEGER_SIZE_DB" , 8));

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

void write_mesh_with_fields(const std::string& filename, stk::io::StkMeshIoBroker &outStkIo, int step, double time, stk::io::DatabasePurpose databasePurpose)
{
    size_t outputFileIndex = outStkIo.create_output_mesh(filename, databasePurpose);

    if (step>0)
    {
        const stk::mesh::FieldVector fields = outStkIo.bulk_data().mesh_meta_data().get_fields();
        for(stk::mesh::FieldBase* field : fields)
        {
            const Ioss::Field::RoleType* fieldRole = stk::io::get_field_role(*field);
            if(fieldRole == nullptr || *fieldRole == Ioss::Field::TRANSIENT)
                outStkIo.add_field(outputFileIndex, *field);
        }
    }

    outStkIo.write_output_mesh(outputFileIndex);

    if (step>0)
    {
        outStkIo.begin_output_step(outputFileIndex, time);
        outStkIo.write_defined_output_fields(outputFileIndex);
        outStkIo.end_output_step(outputFileIndex);
    }
}

void write_mesh_with_fields(const std::string& filename, stk::mesh::BulkData &bulkData, int step, double time, stk::io::DatabasePurpose databasePurpose)
{
    stk::io::StkMeshIoBroker outStkIo;
    outStkIo.set_bulk_data(bulkData);
    write_mesh_with_fields(filename, outStkIo, step, time, databasePurpose);
}

void write_mesh_with_large_ids_and_fields(const std::string& filename, stk::mesh::BulkData &bulkData, int step, double time, stk::io::DatabasePurpose databasePurpose)
{
    stk::io::StkMeshIoBroker outStkIo;
    outStkIo.property_add(Ioss::Property("INTEGER_SIZE_API" , 8));
    outStkIo.property_add(Ioss::Property("INTEGER_SIZE_DB" , 8));
    outStkIo.set_bulk_data(bulkData);
    write_mesh_with_fields(filename, outStkIo, step, time, databasePurpose);
}


}
}
