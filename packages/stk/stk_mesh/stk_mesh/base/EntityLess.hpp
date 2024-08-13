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

#ifndef STK_ENTITYLESS_HPP
#define STK_ENTITYLESS_HPP

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/BulkData.hpp>

namespace stk {
namespace mesh {

class EntityLess {
public:
#ifdef SIERRA_MIGRATION
  EntityLess(const BulkData& mesh);
#else
  inline EntityLess(const BulkData& mesh);
#endif
  /** \brief  Comparison operator */
  inline bool operator()(const Entity lhs, const Entity rhs) const;
  inline bool operator()(const Entity lhs, const EntityKey & rhs) const;
  inline bool operator()( const EntityProc & lhs, const EntityProc & rhs) const;
  inline bool operator()( const EntityProc & lhs, const Entity rhs) const;
  inline bool operator()( const EntityProc & lhs, const EntityKey & rhs) const;
private:
  const BulkData* m_mesh;
#ifdef SIERRA_MIGRATION
  const bool m_shouldSortFacesByNodeIds;
  const EntityRank m_sideRank;
#endif
}; //struct EntityLess

#ifdef SIERRA_MIGRATION
inline bool 
EntityLess::operator()(const Entity lhs, const Entity rhs) const
{
  bool result = false;
  if (m_shouldSortFacesByNodeIds &&
      m_mesh->entity_rank(lhs) == m_sideRank &&
      m_mesh->entity_rank(rhs) == m_sideRank)
  {
    unsigned num_nodes_lhs = m_mesh->count_valid_connectivity(lhs, stk::topology::NODE_RANK);
    unsigned num_nodes_rhs = m_mesh->count_valid_connectivity(rhs, stk::topology::NODE_RANK);
    if (num_nodes_lhs != num_nodes_rhs)
    {    
      result = num_nodes_lhs < num_nodes_rhs;
    }    
    else if (num_nodes_lhs == 0) { 
      result = m_mesh->identifier(lhs) < m_mesh->identifier(rhs);
    }    
    else 
    {    
      const stk::mesh::Entity* nodes_lhs_ptr = m_mesh->begin_nodes(lhs);
      const stk::mesh::Entity* nodes_rhs_ptr = m_mesh->begin_nodes(rhs);
      unsigned i=0; 
      while(i<num_nodes_lhs &&
            (m_mesh->identifier(nodes_lhs_ptr[i]) == m_mesh->identifier(nodes_rhs_ptr[i])))
      {    
        ++i;
      }    
      result = (i<num_nodes_lhs) ?
            (m_mesh->identifier(nodes_lhs_ptr[i]) < m_mesh->identifier(nodes_rhs_ptr[i]))
          : false;
    }    
  }
  else 
  {
    const EntityKey lhs_key = m_mesh->entity_key(lhs);
    const EntityKey rhs_key = m_mesh->entity_key(rhs);
    result = lhs_key < rhs_key;
  }
  return result;
}

#else

inline
EntityLess::EntityLess(const BulkData& mesh) : m_mesh(&mesh) {}

inline bool 
EntityLess::operator()(const Entity lhs, const Entity rhs) const
{
  const EntityKey lhs_key = m_mesh->entity_key(lhs);
  const EntityKey rhs_key = m_mesh->entity_key(rhs);
  return (lhs_key < rhs_key);
}
#endif

inline bool 
EntityLess::operator()(const Entity lhs, const EntityKey & rhs) const
{
  const EntityKey lhs_key = m_mesh->entity_key(lhs);
  return lhs_key < rhs; 
}

inline bool 
EntityLess::operator()( const EntityProc & lhs, const EntityProc & rhs) const
{
  const EntityKey lhs_key = m_mesh->entity_key(lhs.first);
  const EntityKey rhs_key = m_mesh->entity_key(rhs.first);
  return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
}

inline bool 
EntityLess::operator()( const EntityProc & lhs, const Entity rhs) const
{
  const EntityKey lhs_key = m_mesh->entity_key(lhs.first);
  const EntityKey rhs_key = m_mesh->entity_key(rhs);
  return lhs_key < rhs_key;
}

inline bool 
EntityLess::operator()( const EntityProc & lhs, const EntityKey & rhs) const
{
  const EntityKey lhs_key = m_mesh->entity_key(lhs.first);
  return lhs_key < rhs ;
}

} // namespace mesh
} // namespace stk

#endif
