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

#ifndef stk_mesh_impl_CommEntityMods_hpp
#define stk_mesh_impl_CommEntityMods_hpp

#include <stk_util/parallel/CommSparse.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/EntityCommListInfo.hpp>
#include <stk_mesh/base/EntityParallelState.hpp>
#include <vector>

namespace stk {
namespace mesh {

class BulkData;
class EntityCommDatabase;

namespace impl {

class CommEntityMods
{
public:
  enum PackOption { PACK_SHARED = 0, PACK_GHOSTED = 1, PACK_ALL = 2 };

  CommEntityMods(const BulkData& bulkData,
                 const EntityCommDatabase& commDB,
                 const EntityCommListInfoVector& commList);
  virtual ~CommEntityMods();

  void communicate(PackOption packOption);

  const std::vector<EntityParallelState>& get_shared_mods() const { return m_sharedMods; }
  const std::vector<EntityParallelState>& get_ghosted_mods() const { return m_ghostedMods; }

private:
  void pack();
  bool need_to_send();
  void pack_entity_mods();
  void unpack();

  const BulkData& m_bulkData;
  CommSparse m_commSparse;
  const EntityCommListInfoVector& m_commList;
  const EntityCommDatabase& m_commDB;
  PackOption m_packOption;
  std::vector<EntityParallelState> m_sharedMods; 
  std::vector<EntityParallelState> m_ghostedMods; 
};

void make_unique(std::vector<EntityParallelState>& pllStates);

}}} // end namepsace stk mesh impl


#endif
