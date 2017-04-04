// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stk_mesh/base/Ghosting.hpp>
#include <ostream>                      // for operator<<, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Types.hpp"      // for EntityProc, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter

namespace stk {
namespace mesh {

void Ghosting::send_list( std::vector< EntityProc > & v ) const
{
  for ( const EntityCommListInfo& commListInfo : m_mesh.internal_comm_list() ) {
    if ( commListInfo.owner == m_mesh.parallel_rank() ) {
      for ( PairIterEntityComm ec = m_mesh.internal_entity_comm_map(commListInfo.key) ; ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == m_ordinal ) {
          v.push_back( EntityProc( commListInfo.entity , ec->proc ) );
        }
      }
    }
  }
}

void Ghosting::receive_list( std::vector<EntityKey> & v ) const
{
  for ( const EntityCommListInfo& commListInfo : m_mesh.internal_comm_list() ) {
    if ( commListInfo.owner != m_mesh.parallel_rank() ) {
      for ( PairIterEntityComm ec = m_mesh.internal_entity_comm_map(commListInfo.key) ; ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == m_ordinal ) {
          v.push_back(commListInfo.key);
        }
      }
    }
  }
}

std::ostream& Ghosting::operator<<(std::ostream& out) const
{
  out << "Ghosting object: name: " << name()
      << ", ordinal: " << ordinal() << "\n";

  out << "  Locally owned entities ghosted on other processors (send list):\n";

  for ( const EntityCommListInfo& commListInfo : m_mesh.internal_comm_list() ) {
    if ( commListInfo.owner == m_mesh.parallel_rank() ) {
      for ( PairIterEntityComm ec = m_mesh.internal_entity_comm_map(commListInfo.key) ; ! ec.empty() ; ++ec ) {
        if ( ec->ghost_id == m_ordinal ) {
          out << "    ";
          out << commListInfo.key.id();
          out << ", sending ghost to " << ec->proc << ", status is: "
              << m_mesh.state(commListInfo.entity) << "\n";
        }
      }
    }
  }

  out << "  Entities ghosted on this processor from the owner (recv list):\n";
  for ( const EntityCommListInfo& commListInfo : m_mesh.internal_comm_list() ) {
    if ( commListInfo.owner != m_mesh.parallel_rank() ) {
      for ( PairIterEntityComm ec = m_mesh.internal_entity_comm_map(commListInfo.key); !ec.empty(); ++ec ) {
        if ( ec->ghost_id == m_ordinal ) {
          out << "    ";
          out << commListInfo.key.id();
          out << ", owner of ghost is " << commListInfo.owner
              << ", status is: " << m_mesh.state(commListInfo.entity) << "\n";
        }
      }
    }
  }
  return out;
}

std::ostream& operator<<(std::ostream& out, const Ghosting& rhs)
{
  return rhs.operator<<(out);
}

}
}


