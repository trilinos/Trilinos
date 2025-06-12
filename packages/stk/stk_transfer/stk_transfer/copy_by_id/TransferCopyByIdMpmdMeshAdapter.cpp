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

#include "TransferCopyByIdMpmdMeshAdapter.hpp"
#include "TransferCopyByIdMeshAdapter.hpp"

namespace stk {
namespace transfer {

  TransferCopyByIdMpmdMeshAdapter::TransferCopyByIdMpmdMeshAdapter(stk::ParallelMachine pm, const int num_fields)
  : m_comm(pm),
    m_num_fields(num_fields)
  {
  }

  const void* TransferCopyByIdMpmdMeshAdapter::field_data(const Mesh_ID & /*id*/, const unsigned /*field_index*/) const
  {
    return nullptr;
  }
  void* TransferCopyByIdMpmdMeshAdapter::field_data(const Mesh_ID & /*id*/, const unsigned /*field_index*/)
  {
    return nullptr;
  }
  unsigned TransferCopyByIdMpmdMeshAdapter::field_data_size(const Mesh_ID & /*id*/, const unsigned /*field_index*/) const
  {
    return 0;
  }
  unsigned TransferCopyByIdMpmdMeshAdapter::num_fields() const
  {
    return m_num_fields;
  }

  ParallelMachine TransferCopyByIdMpmdMeshAdapter::comm() const
  {
    return m_comm;
  }

  const TransferCopyByIdMeshAdapter::MeshIDVector & TransferCopyByIdMpmdMeshAdapter::get_mesh_ids() const
  {
    return m_empty_ids;
  }
  bool TransferCopyByIdMpmdMeshAdapter::is_locally_owned(const Mesh_ID & /*id*/) const
  {
    return false;
  }
  void TransferCopyByIdMpmdMeshAdapter::centroid(const Mesh_ID & /*id*/, [[maybe_unused]] double coords[3]) const
  {
  }
  std::string TransferCopyByIdMpmdMeshAdapter::print_mesh_id(const Mesh_ID& /*id*/) const
  {
    return "";
  }
  std::string TransferCopyByIdMpmdMeshAdapter::field_name(const unsigned /*index*/) const
  {
    return "";
  }

  DataTypeKey::data_t TransferCopyByIdMpmdMeshAdapter::get_field_type(const unsigned /*fieldIndex*/) const
  { 
    return DataTypeKey::data_t::INVALID_TYPE;
  }
}
}
