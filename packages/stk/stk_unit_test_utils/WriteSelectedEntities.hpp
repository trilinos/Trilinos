// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
 // 
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
 // 
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef WRITESELECTEDENTITIES_HPP_
#define WRITESELECTEDENTITIES_HPP_

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <stk_mesh/base/Selector.hpp>  // for Selector
#include <stk_mesh/base/Types.hpp>     // for FieldVector
#include <string>                      // for string
namespace stk { namespace mesh { class BulkData; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

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

