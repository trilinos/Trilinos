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
#ifndef STK_BALANCE_M2N_OUTPUTMESH_HPP
#define STK_BALANCE_M2N_OUTPUTMESH_HPP

#include <stk_balance/m2n/TransientFieldTransferById.hpp>
#include <stk_balance/m2n/OutputSerializerBulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <vector>

namespace stk { namespace balance { namespace m2n { class InputMesh; }}}

namespace stk {
namespace balance {
namespace m2n {

class OutputMesh
{
public:
  OutputMesh(const InputMesh& inputMesh,
             const std::vector<unsigned>& targetSubdomains);
  ~OutputMesh() = default;

  void transfer_and_write();

  OutputSerializerBulkData& get_bulk() { return m_bulk; }
  stk::mesh::MetaData& get_meta() { return m_meta; }

private:
  void clone_input_mesh();
  void move_subdomain_to_owning_processor();

  const InputMesh& m_inputMesh;
  const std::vector<unsigned>& m_targetSubdomains;
  stk::mesh::MetaData m_meta;
  OutputSerializerBulkData m_bulk;
};

}}}
#endif
