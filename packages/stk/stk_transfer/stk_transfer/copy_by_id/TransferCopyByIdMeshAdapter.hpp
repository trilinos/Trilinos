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

#ifndef STK_TRANSFER_MESHBASE_HPP
#define STK_TRANSFER_MESHBASE_HPP

#include <inttypes.h>
#include <vector>
#include <string>
#include <stk_util/parallel/Parallel.hpp>
#include "TransferCopyTranslator.hpp"

namespace stk {
namespace transfer {

class TransferCopyByIdMeshAdapter {
public:
  typedef uint64_t Mesh_ID;
  typedef std::vector<Mesh_ID> MeshIDVector;

  virtual const void* field_data(const Mesh_ID & id, const unsigned field_index) const = 0;
  virtual       void* field_data(const Mesh_ID & id, const unsigned field_index)       = 0;
  virtual std::string field_name(const unsigned field_index) const = 0;
  virtual unsigned field_data_size(const Mesh_ID & id, const unsigned field_index) const = 0;
  virtual unsigned num_fields() const = 0;
  virtual ParallelMachine comm() const = 0;
  virtual void begin_transfer() const {}
  virtual void end_transfer() const {}

  virtual DataTypeKey::data_t get_field_type(const unsigned fieldIndex) const = 0;

  virtual const MeshIDVector & get_mesh_ids() const = 0;
  virtual bool is_locally_owned(const Mesh_ID & id) const = 0;
  virtual void centroid(const Mesh_ID & id, double coords[3]) const = 0;
  virtual std::string print_mesh_id(const Mesh_ID& id) const = 0;
  virtual ~TransferCopyByIdMeshAdapter() = default;
};

} } // namespace stk transfer



#endif // STK_TRANSFER_MESHBASE_HPP

