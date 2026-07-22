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

#ifndef STK_STK_TOOLS_MESH_CLONE_REPLACEBULKDATA_HPP_
#define STK_STK_TOOLS_MESH_CLONE_REPLACEBULKDATA_HPP_

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_tools/mesh_clone/ReplaceBulkDataImpl.hpp>
#include <functional>

namespace stk {
namespace tools {

template<typename T>
void replace_bulk_data(const stk::mesh::BulkData & inMesh, T & outMesh, std::function<void(T& outMesh_)> op)
{
  size_t previousSyncCount = outMesh.synchronized_count();
  std::shared_ptr<stk::mesh::MetaData> outMeta = outMesh.mesh_meta_data_ptr();
  outMesh.~T();

  const stk::mesh::BulkData::AutomaticAuraOption aura_option =
    inMesh.is_automatic_aura_on() ? stk::mesh::BulkData::AUTO_AURA :  stk::mesh::BulkData::NO_AUTO_AURA;

  new (&outMesh) T(outMeta,
                   inMesh.parallel(),
                   aura_option
#ifdef SIERRA_MIGRATION
                   ,inMesh.add_fmwk_data()
#endif
    );

  op(outMesh);

  outMesh.modification_begin();
  impl::clone_bulk_data_entities(inMesh, outMesh);
  outMesh.modification_end();

  impl::copy_field_data(inMesh, outMesh);

  outMesh.m_meshModification.set_sync_count(previousSyncCount+1);
}

template<typename T>
void replace_bulk_data(const stk::mesh::BulkData & inMesh, T & outMesh)
{
  std::function<void(T& outMesh_)> op  = [](T& /*outMesh_*/) {};
  replace_bulk_data(inMesh, outMesh, op);
}

} // namespace tools
} // namespace stk

#endif /* STK_STK_TOOLS_MESH_CLONE_REPLACEBULKDATA_HPP_ */
