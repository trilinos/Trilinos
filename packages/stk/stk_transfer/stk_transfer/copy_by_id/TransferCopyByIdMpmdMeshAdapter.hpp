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

#ifndef STK_STK_TRANSFER_STK_TRANSFER_COPY_BY_ID_TRANSFERCOPYBYIDMPMDMESHADAPTER_HPP_
#define STK_STK_TRANSFER_STK_TRANSFER_COPY_BY_ID_TRANSFERCOPYBYIDMPMDMESHADAPTER_HPP_

#include "TransferCopyByIdMeshAdapter.hpp"
#include <inttypes.h>
#include <stk_util/parallel/Parallel.hpp>
#include <string>
#include <vector>

namespace stk {
namespace transfer {

//For mpmd transfer cases, the stk mesh meta and bulk won't ven exist on some processors,
//depending on the communicator. In those cases, just serve up empty get_mesh_ids(),
//throw on everything else
class TransferCopyByIdMpmdMeshAdapter : public TransferCopyByIdMeshAdapter {
public:
  TransferCopyByIdMpmdMeshAdapter(stk::ParallelMachine pm, const int num_fields);
  const void* field_data(const Mesh_ID & id, const unsigned field_index) const override;
  void* field_data(const Mesh_ID & id, const unsigned field_index) override;
  std::string field_name(const unsigned field_index) const override;
  unsigned field_data_size(const Mesh_ID & id, const unsigned field_index) const override;
  unsigned num_fields() const override;
  ParallelMachine comm() const override;

  const MeshIDVector & get_mesh_ids() const override;
  bool is_locally_owned(const Mesh_ID & id) const override;
  void centroid(const Mesh_ID & id, double coords[3]) const override;
  std::string print_mesh_id(const Mesh_ID& id) const override;
  virtual ~TransferCopyByIdMpmdMeshAdapter() = default;
  DataTypeKey::data_t get_field_type(const unsigned fieldIndex) const override;
  
private:
  stk::ParallelMachine m_comm;
  MeshIDVector m_empty_ids;
  const int m_num_fields;
};

} } // namespace stk transfer


#endif /* STK_STK_TRANSFER_STK_TRANSFER_COPY_BY_ID_TRANSFERCOPYBYIDMPMDMESHADAPTER_HPP_ */
