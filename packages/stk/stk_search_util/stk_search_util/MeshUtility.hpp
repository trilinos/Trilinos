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

#ifndef STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MESHUTILITY_HPP_
#define STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MESHUTILITY_HPP_

#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, operator<<
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for EntityRank, EntityId, PartVector
#include "stk_search/Box.hpp"           // for Box
#include "stk_search_util/spmd/EntityKeyPair.hpp"
#include "stk_topology/topology.hpp"    // for topology, topology::INVALID_RANK

#include <Kokkos_Core_fwd.hpp>

#include <ostream>                      // for operator<<, ostream, basic_os...
#include <string>                       // for string
#include <utility>                      // for pair
#include <vector>                       // for vector

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }

namespace stk {
namespace search {

inline int get_key(int ekey) { return ekey; }
inline int get_id(int ekey) { return ekey; }
inline int get_rank(int /*ekey*/) { return stk::topology::INVALID_RANK; }

inline stk::mesh::EntityKey get_key(stk::mesh::EntityKey ekey) { return ekey; }
inline stk::mesh::EntityId get_id(stk::mesh::EntityKey ekey) { return ekey.id(); }
inline stk::mesh::EntityRank get_rank(stk::mesh::EntityKey ekey) { return ekey.rank(); }

inline stk::mesh::EntityKey get_key(Kokkos::pair<stk::mesh::EntityKey,int> ekey) { return ekey.first; }
inline stk::mesh::EntityId get_id(Kokkos::pair<stk::mesh::EntityKey,int> ekey) { return ekey.first.id(); }
inline stk::mesh::EntityRank get_rank(Kokkos::pair<stk::mesh::EntityKey,int> ekey) { return ekey.first.rank(); }
inline std::ostream& operator<<(std::ostream& os, const Kokkos::pair<stk::mesh::EntityKey,int> ekey)
{
  os << "{" << ekey.first << ", INDEX: " << ekey.second << "}";
  return os;
}

inline stk::mesh::EntityKey get_key(std::pair<stk::mesh::EntityKey,int> ekey) { return ekey.first; }
inline stk::mesh::EntityId get_id(std::pair<stk::mesh::EntityKey,int> ekey) { return ekey.first.id(); }
inline stk::mesh::EntityRank get_rank(std::pair<stk::mesh::EntityKey,int> ekey) { return ekey.first.rank(); }
inline std::ostream& operator<<(std::ostream& os, std::pair<stk::mesh::EntityKey,int> ekey)
{
  os << "{" << ekey.first << ", INDEX: " << ekey.second << "}";
  return os;
}

inline stk::search::spmd::EntityKeyPair get_key(Kokkos::pair<stk::search::spmd::EntityKeyPair,int> ekey) { return ekey.first; }
inline stk::mesh::EntityId get_id(Kokkos::pair<stk::search::spmd::EntityKeyPair,int> ekey) { return ekey.first.id(); }
inline stk::mesh::EntityRank get_rank(Kokkos::pair<stk::search::spmd::EntityKeyPair,int> ekey) { return ekey.first.rank(); }
inline std::ostream& operator<<(std::ostream& os, const Kokkos::pair<stk::search::spmd::EntityKeyPair,int> ekey)
{
  os << "{" << ekey.first << ", INDEX: " << ekey.second << "}";
  return os;
}

inline stk::search::spmd::EntityKeyPair get_key(std::pair<stk::search::spmd::EntityKeyPair,int> ekey) { return ekey.first; }
inline stk::mesh::EntityId get_id(std::pair<stk::search::spmd::EntityKeyPair,int> ekey) { return ekey.first.id(); }
inline stk::mesh::EntityRank get_rank(std::pair<stk::search::spmd::EntityKeyPair,int> ekey) { return ekey.first.rank(); }
inline std::ostream& operator<<(std::ostream& os, std::pair<stk::search::spmd::EntityKeyPair,int> ekey)
{
  os << "{" << ekey.first << ", INDEX: " << ekey.second << "}";
  return os;
}

stk::search::Box<double> get_mesh_bounding_box(const stk::mesh::FieldBase* coords, const stk::mesh::BulkData& bulk);

stk::mesh::EntityRank get_objects_rank(const stk::mesh::EntityRank objectType, const stk::mesh::PartVector& parts);

unsigned get_number_of_parametric_coordinates(stk::topology topo);

bool is_valid_entity_rank(const stk::mesh::MetaData& meta, const stk::mesh::EntityRank rank);

bool part_has_proper_entity_rank(const stk::mesh::Part* part);

std::vector<std::string> get_part_membership(const stk::mesh::BulkData& bulk, const stk::mesh::Entity e, const stk::mesh::PartVector& parts);

std::vector<std::string> get_part_membership(const stk::mesh::BulkData& bulk, const stk::mesh::EntityKey k, const stk::mesh::PartVector& parts);

stk::mesh::Selector get_objects_selector(const stk::mesh::MetaData& stkMeta, const stk::mesh::PartVector& parts,
                                         const stk::mesh::Selector* activeSelector = nullptr, bool includeGhosts = false);

std::string get_time_stamp();

}
}

#endif /* STK_STK_SEARCH_UTIL_STK_SEARCH_UTIL_MESHUTILITY_HPP_ */
