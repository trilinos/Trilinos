// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_unit_test_utils/WriteSelectedEntities.hpp>
#include <stk_mesh/base/BulkData.hpp>         // for BulkData
#include <stk_mesh/base/Types.hpp>            // for FieldVector
#include <stk_mesh/baseImpl/DebugWriter.hpp>  // for DebugWriter, etc
#include <string>                             // for string, basic_string
#include "stk_mesh/base/Selector.hpp"         // for Selector
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

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

namespace simple_fields {

void write_selected_entities(const stk::mesh::BulkData& mesh_bulk,
                             stk::mesh::Selector selector,
                             const stk::mesh::FieldVector &fields,
                             const std::string& callingFile, int lineNumber, const std::string& basename) {
  stk::debug::write_selected_entities(mesh_bulk, selector, fields, callingFile, lineNumber, basename);
}

void write_selected_entities_with_internal_parts(const stk::mesh::BulkData& mesh_bulk,
                                                 stk::mesh::Selector selector,
                                                 const stk::mesh::FieldVector &fields,
                                                 const std::string& calling_file, int line_number, const std::string& basename) {
  stk::debug::write_selected_entities_with_internal_parts(mesh_bulk, selector, fields, calling_file, line_number, basename);
}

void write_int_value(const std::string &tagString, int scalar, const std::string& calling_file, int numProcs, int localProc, int line_number) {
  stk::debug::write_int_value(tagString, scalar, calling_file, numProcs, localProc, line_number);
}

void write_real_value(const std::string &tagString, double scalar, const std::string& calling_file, int numProcs, int localProc, int line_number) {
  stk::debug::write_real_value(tagString, scalar, calling_file, numProcs, localProc, line_number);
}

void write_string_value(const std::string &tagString, const std::string &str, const std::string& callingFile, int numProcs, int localProc, int lineNumber) {
  stk::debug::write_string_value(tagString, str, callingFile, numProcs, localProc, lineNumber);
}

void write_meta_data(const stk::mesh::BulkData& meshBulk,
                     const std::string& callingFile, int lineNumber, const std::string& basename) {
  stk::debug::write_meta_data(meshBulk, callingFile, lineNumber, basename);
}

} // namespace simple_fields

}
}

