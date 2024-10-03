// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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
//

#ifndef STK_PERFORMANCE_MULTI_BLOCK_HPP
#define STK_PERFORMANCE_MULTI_BLOCK_HPP

namespace stk {

namespace mesh { class MetaData; }
namespace mesh { class BulkData; }

namespace performance_tests {

void setup_multiple_blocks(stk::mesh::MetaData& meta, unsigned numBlocks);
void setup_elem_fields_on_blocks(stk::mesh::MetaData& meta, unsigned numFields, bool allFieldsSameSize = false);
void setup_sidesets_between_blocks(stk::mesh::MetaData& meta);
void setup_sidesets_for_blocks(stk::mesh::MetaData& meta);
void move_elements_to_other_blocks(stk::mesh::BulkData& bulk, unsigned numElemsInDimX);
void move_elements_to_other_contiguous_blocks(stk::mesh::BulkData& bulk, unsigned numElemsInDimX);
void fill_sidesets_between_blocks(stk::mesh::BulkData& bulk);

}}

#endif
