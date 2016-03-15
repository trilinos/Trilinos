#ifndef WRITESELECTEDENTITIES_HPP_
#define WRITESELECTEDENTITIES_HPP_

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <string>

namespace stk {
namespace debug {

void write_selected_entities(const stk::mesh::BulkData& mesh_bulk,
                             stk::mesh::Selector selector,
                             const stk::mesh::FieldVector &fields,
                             const std::string& callingFile, int lineNumber, const std::string& basename = "selected");

void write_selected_entities_with_internal_parts(const stk::mesh::BulkData& mesh_bulk,
                                                 stk::mesh::Selector selector,
                                                 const stk::mesh::FieldVector &fields,
                                                 const std::string& calling_file, int line_number, const std::string& basename = "internals");
void write_int_value(const std::string &tagString, int scalar, const std::string& calling_file, int numProcs, int localProc, int line_number);
void write_real_value(const std::string &tagString, double scalar, const std::string& calling_file, int numProcs, int localProc, int line_number);
void write_string_value(const std::string &tagString, const std::string &str, const std::string& callingFile, int numProcs, int localProc, int lineNumber);
void write_meta_data(const stk::mesh::BulkData& meshBulk,
                     const std::string& callingFile, int lineNumber, const std::string& basename);


}//namespace debug
}//namespace stk

#endif /* WRITESELECTEDENTITIES_HPP_ */

