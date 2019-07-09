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
#ifndef STK_TOOLS_CustomAura_hpp_
#define STK_TOOLS_CustomAura_hpp_

namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Ghosting; } }

#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Types.hpp"

namespace stk {
namespace tools {

  inline void fill_list_of_entities_to_send_for_aura_like_ghosting(stk::mesh::BulkData& bulkData,
                                                            stk::mesh::Selector const& selector,
                                                            stk::mesh::EntityProcVec &entitiesToGhost)
  {
      const stk::mesh::BucketVector& buckets = bulkData.get_buckets(stk::topology::NODE_RANK, selector);
      for(const stk::mesh::Bucket *bucket : buckets)
      {
          for(stk::mesh::Entity shared_node : *bucket)
          {
              unsigned num_elements = bulkData.num_elements(shared_node);
              const stk::mesh::Entity* elements = bulkData.begin_elements(shared_node);
              for(unsigned i=0;i<num_elements;++i)
              {
                  if(bulkData.bucket(elements[i]).owned())
                  {
                      std::vector<int> comm_shared_procs;
                      bulkData.comm_shared_procs(bulkData.entity_key(shared_node), comm_shared_procs);
                      for(int proc : comm_shared_procs )
                      {
                          entitiesToGhost.push_back(stk::mesh::EntityProc(elements[i], proc));
                      }
                  }
              }
          }
      }
  }

  inline stk::mesh::Ghosting * create_custom_aura(stk::mesh::BulkData & stkMeshBulkData,
                                                  stk::mesh::Selector const& selector,
                                                  const std::string & name)
  {
      stk::mesh::Ghosting * customAura = nullptr;
      if (!stkMeshBulkData.is_automatic_aura_on())
      {
          stkMeshBulkData.modification_begin();

          customAura = &stkMeshBulkData.create_ghosting(name);
          stk::mesh::EntityProcVec entitiesToGhost;
          fill_list_of_entities_to_send_for_aura_like_ghosting(stkMeshBulkData, selector, entitiesToGhost);
          stkMeshBulkData.change_ghosting(*customAura, entitiesToGhost);

          stkMeshBulkData.modification_end();
      }

      return customAura;
  }

  inline void destroy_custom_aura(stk::mesh::BulkData & stkMeshBulkData, stk::mesh::Ghosting * customAura)
  {
      if (nullptr != customAura)
      {
          stkMeshBulkData.modification_begin();
          stkMeshBulkData.destroy_ghosting(*customAura);
          stkMeshBulkData.modification_end();
      }
  }

}
}

#endif
