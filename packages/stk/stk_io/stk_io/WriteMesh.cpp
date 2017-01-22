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

}
}
