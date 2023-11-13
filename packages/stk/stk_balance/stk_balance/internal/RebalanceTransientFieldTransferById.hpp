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

#ifndef REBALANCETRANSIENTFIELDTRANSFERBYID_HPP
#define REBALANCETRANSIENTFIELDTRANSFERBYID_HPP

#include "stk_mesh/base/Types.hpp"
#include "stk_io/IossBridge.hpp"
#include <vector>
#include <map>

namespace stk { namespace io { class StkMeshIoBroker; } }
namespace stk { namespace transfer_utils { class TransientTransferByIdForRank; } }
namespace stk { namespace balance { class SubdomainWriterBase; } }
namespace stk { namespace balance { class OutputSerializerBulkData; } }

namespace stk {
namespace balance {

constexpr unsigned INVALID_SUBDOMAIN = std::numeric_limits<unsigned>::max();

inline bool is_valid_subdomain(unsigned subdomain) {
  return subdomain != INVALID_SUBDOMAIN;
}

class RebalanceTransientFieldTransferById
{
public:
  RebalanceTransientFieldTransferById() = delete;
  RebalanceTransientFieldTransferById(stk::io::StkMeshIoBroker& inputBroker, unsigned numSubDomain);

  ~RebalanceTransientFieldTransferById();

  void setup_subdomain(OutputSerializerBulkData& bulk, const std::string& filename,
                       unsigned subdomain, const stk::io::EntitySharingInfo& nodeSharingInfo,
                       int global_num_nodes, int global_num_elems);

  void transfer_transient_data(unsigned subdomain);
  void transfer_and_write_transient_data(unsigned subdomain);

  SubdomainWriterBase& get_subdomain_writer(unsigned subdomain);

private:
  void initialize_transfer(unsigned subDomainInfo);

  stk::io::StkMeshIoBroker& m_inputBroker;
  unsigned m_numSubDomain;
  std::map<unsigned, SubdomainWriterBase*> m_subdomainWriters;
  std::map<unsigned, std::vector<stk::transfer_utils::TransientTransferByIdForRank*>> m_subdomainTransfers;
  std::vector<stk::mesh::EntityRank> m_entityRanks;
};

}
}

#endif // REBALANCETRANSIENTFIELDTRANSFERBYID_HPP
