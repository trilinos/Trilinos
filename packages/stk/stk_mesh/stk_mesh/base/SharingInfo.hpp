// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
//
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
//
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

#ifndef STK_STK_MESH_STK_MESH_BASE_SHARINGINFO_HPP_
#define STK_STK_MESH_STK_MESH_BASE_SHARINGINFO_HPP_

#include <stddef.h>                     // for size_t

#include <algorithm>
#include <vector>

#include <stk_mesh/base/Entity.hpp>     // for Entity, etc
#include <stk_util/util/SortAndUnique.hpp>

namespace stk {
namespace mesh {

struct sharing_info
{
    stk::mesh::Entity m_entity;
    int m_sharing_proc;
    int m_owner;

    sharing_info(stk::mesh::Entity entity, int sharing_proc, int owner) :
        m_entity(entity), m_sharing_proc(sharing_proc), m_owner(owner) {}

    sharing_info(const sharing_info&) = default;

    sharing_info():
      m_entity(stk::mesh::Entity()), m_sharing_proc(-1), m_owner(-1) {}
};

struct SharingInfoLess
{
    bool operator()(const stk::mesh::sharing_info &a, const stk::mesh::sharing_info &b)
    {
        if(a.m_entity == b.m_entity)
            return a.m_owner < b.m_owner;
        else
            return a.m_entity < b.m_entity;
    }
};

struct SharingInfoLessByEntity
{
    bool operator()(const stk::mesh::sharing_info &a, const stk::mesh::sharing_info &b)
    {
      return a.m_entity < b.m_entity;
    }
};

struct SharingInfoEqual
{
    bool operator()(const stk::mesh::sharing_info &a, const stk::mesh::sharing_info &b)
    {
       return (a.m_entity == b.m_entity &&  a.m_owner == b.m_owner && a.m_sharing_proc == b.m_sharing_proc);
    }
};

inline void update_sharing_info_ownership(std::vector<stk::mesh::sharing_info>& sharedModified)
{
  stk::util::sort_and_unique(sharedModified, SharingInfoLess(), SharingInfoEqual());
  if(!sharedModified.empty()) {
    for(size_t i=1; i<sharedModified.size(); i++) {
      if(sharedModified[i].m_entity == sharedModified[i-1].m_entity) {
        sharedModified[i].m_owner = sharedModified[i-1].m_owner;
      }
    }
  }
}

} // namespace mesh
} // namespace stk

#endif /* STK_STK_MESH_STK_MESH_BASE_SHARINGINFO_HPP_ */
