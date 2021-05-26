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
#ifndef M2NDECOMPOSER_HPP
#define M2NDECOMPOSER_HPP

#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace balance { class M2NBalanceSettings; } }

namespace stk {
namespace balance {
namespace internal {

class M2NDecomposer
{
public:
  M2NDecomposer(stk::mesh::BulkData & bulkData,
                const M2NBalanceSettings & balanceSettings);
  virtual ~M2NDecomposer() = default;

  virtual stk::mesh::EntityProcVec get_partition();

  virtual std::vector<unsigned> map_new_subdomains_to_original_processors();

  unsigned num_required_subdomains_for_each_proc();

protected:
  stk::mesh::BulkData & m_bulkData;
  const stk::balance::M2NBalanceSettings & m_balanceSettings;
};

class M2NDecomposerNested : public M2NDecomposer
{
public:
  M2NDecomposerNested(stk::mesh::BulkData & bulkData,
                      const stk::balance::M2NBalanceSettings & balanceSettings);
  virtual ~M2NDecomposerNested() override = default;

  virtual stk::mesh::EntityProcVec get_partition() override;

  virtual std::vector<unsigned> map_new_subdomains_to_original_processors() override;

private:
  stk::mesh::EntityProcVec get_partition_for_subdomain(int subdomainId);
  std::string get_initial_subdomain_part_name(int subdomainId);
  void declare_all_initial_subdomain_parts();
  void move_entities_into_initial_subdomain_part();

  int m_numFinalSubdomainsPerProc;
};

}
}
}
#endif // M2NDECOMPOSER_HPP
