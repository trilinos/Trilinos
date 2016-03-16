#include <stk_unit_test_utils/WriteSelectedEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>                     // for BulkData
#include <stddef.h>                                  // for size_t
#include <fstream>                                   // for operator<<, etc
#include <stk_mesh/base/DataTraits.hpp>              // for DataTraits
#include <stk_mesh/base/FieldTraits.hpp>
#include <stk_mesh/base/FieldBase.hpp>               // for FieldBase, etc
#include <stk_mesh/base/Types.hpp>                   // for EntityVector, etc
#include <stk_util/environment/Env.hpp>              // for parallel_rank, etc
#include <vector>                                    // for vector
#include "stk_mesh/base/Bucket.hpp"                  // for Bucket
#include "stk_mesh/base/BulkData.hpp"                // for BulkData
#include "stk_mesh/base/BulkDataInlinedMethods.hpp"
#include "stk_mesh/base/Entity.hpp"                  // for Entity
#include "stk_mesh/base/EntityKey.hpp"               // for operator<<
#include "stk_mesh/base/GetEntities.hpp"
#include "stk_mesh/base/Part.hpp"                    // for Part
#include "stk_mesh/base/Selector.hpp"                // for Selector
#include "stk_topology/topology.hpp"                 // for operator<<, etc
#include "stk_util/diag/String.hpp"                  // for operator<<
#include <stk_mesh/baseImpl/DebugWriter.hpp>

namespace stk {
namespace debug {

void write_int_value(const std::string &tagString, int scalar, const std::string& callingFile, int lineNumber, int numProcs, int localProc)
{
    DebugWriter writer(callingFile, lineNumber, numProcs, localProc);
    writer.write_to_file(tagString, scalar);
}

void write_real_value(const std::string &tagString, double scalar, const std::string& callingFile, int lineNumber, int numProcs, int localProc)
{
    DebugWriter writer(callingFile, lineNumber, numProcs, localProc);
    writer.write_to_file(tagString, scalar);
}

void write_string_value(const std::string &tagString, const std::string &str, const std::string& callingFile, int lineNumber, int numProcs, int localProc)
{
    DebugWriter writer(callingFile, lineNumber, numProcs, localProc);
    writer.write_to_file(tagString, str);
}


void write_selected_entities(const stk::mesh::BulkData& meshBulk,
                             stk::mesh::Selector selector,
                             const stk::mesh::FieldVector &fields,
                             const std::string& callingFile, int lineNumber, const std::string& basename)
{
    DebugWriter writer(callingFile, lineNumber, meshBulk.parallel_size(), meshBulk.parallel_rank(), basename);
    writer.write_selected_entities(meshBulk, selector, fields);
}

void write_meta_data(const stk::mesh::BulkData& meshBulk,
                     const std::string& callingFile, int lineNumber, const std::string& basename)
{
    DebugWriterWithInternalParts writer(callingFile, lineNumber, meshBulk.parallel_size(), meshBulk.parallel_rank(), basename);
    writer.write_meta_data(meshBulk);
}

void write_selected_entities_with_internal_parts(const stk::mesh::BulkData& meshBulk,
                             stk::mesh::Selector selector,
                             const stk::mesh::FieldVector &fields,
                             const std::string& callingFile, int lineNumber, const std::string& basename)
{
    DebugWriterWithInternalParts writer(callingFile, lineNumber, meshBulk.parallel_size(), meshBulk.parallel_rank(), basename);
    writer.write_selected_entities(meshBulk, selector, fields);
}

}
}

