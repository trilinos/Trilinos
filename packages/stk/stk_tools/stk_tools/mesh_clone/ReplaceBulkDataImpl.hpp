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

#ifndef STK_STK_TOOLS_MESH_CLONE_REPLACEBULKDATAIMPL_HPP_
#define STK_STK_TOOLS_MESH_CLONE_REPLACEBULKDATAIMPL_HPP_

#include <stk_mesh/base/Types.hpp>

namespace stk {
namespace tools {

namespace impl {
  void clone_bulk_data_entities(const stk::mesh::BulkData & inMesh, stk::mesh::BulkData & outMesh);
  void copy_field_data(const stk::mesh::BulkData & inMesh, stk::mesh::BulkData & outMesh);
  stk::mesh::Selector translate_selector(const stk::mesh::Selector & in_selector, const stk::mesh::MetaData & out_meta);
  void translate_parts(const stk::mesh::PartVector & inParts, const stk::mesh::MetaData & outMeta, stk::mesh::PartVector & outParts);
  void clone_meta_data_parts_and_fields(const stk::mesh::MetaData & in_meta, stk::mesh::MetaData & out_meta);
  void copy_field_data(const stk::mesh::BulkData & inMesh, std::vector<stk::mesh::FieldBase*>& inFieldsToCopyFrom, stk::mesh::BulkData & outMesh, std::vector<stk::mesh::FieldBase*>& outFieldsToCopyTo);
}

} // namespace tools
} // namespace stk

#endif /* STK_STK_TOOLS_MESH_CLONE_REPLACEBULKDATAIMPL_HPP_ */
