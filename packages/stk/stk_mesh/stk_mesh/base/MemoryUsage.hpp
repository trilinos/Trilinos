// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
// 

#ifndef stk_mesh_MemoryUsage_hpp
#define stk_mesh_MemoryUsage_hpp

//----------------------------------------------------------------------

#include <stddef.h>                     // for size_t
#include <iosfwd>                       // for ostream
#include <string>                       // for string
#include <vector>                       // for vector
namespace stk { namespace mesh { class BulkData; } }


namespace stk {
namespace mesh {

//----------------------------------------------------------------------
/** \addtogroup stk_mesh_module
 *  \{
 */

struct MemoryUsage {
  unsigned num_fields;
  unsigned field_bytes;
  unsigned num_parts;
  unsigned part_bytes;
  std::vector<std::string> entity_rank_names;
  std::vector<unsigned> entity_counts;
  unsigned bytes_per_entity;
  std::vector<unsigned> downward_relation_counts;
  std::vector<unsigned> upward_relation_counts;
  unsigned bytes_per_relation;
  std::vector<unsigned> bucket_counts;
  std::vector<unsigned> bucket_bytes;
  size_t total_bytes;
};

void compute_memory_usage(const BulkData& bulk, MemoryUsage& mem_usage);

void print_memory_usage(const MemoryUsage& mem_usage, std::ostream& os);

//----------------------------------------------------------------------


} // namespace mesh
} // namespace stk

#endif
