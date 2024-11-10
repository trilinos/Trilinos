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

#ifndef STK_MESH_BUCKET_CONNECTIVITY_HPP
#define STK_MESH_BUCKET_CONNECTIVITY_HPP

#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_util/util/StridedArray.hpp>
#include "stk_util/util/ReportHandler.hpp"

namespace stk {
namespace mesh {

class BulkData;

using ConnectedNodes    = util::StridedArray<const stk::mesh::Entity>;
using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;
using ConnectedOrdinals = util::StridedArray<const stk::mesh::ConnectivityOrdinal>;
using Permutations      = util::StridedArray<const stk::mesh::Permutation>;

namespace impl {


template <typename Connectivity>
inline void check_bucket_ordinal(unsigned bucket_ordinal, Connectivity const* connectivity)
{
  STK_ThrowAssertMsg(bucket_ordinal < connectivity->size(),
                 "bucket_ordinal " << bucket_ordinal << " is out of range, bucket size is " << connectivity->size());
}

template<EntityRank TargetRank, ConnectivityType >
class BucketConnectivity;

template<typename VecType>
size_t capacity_in_bytes(const VecType& v)
{
  return sizeof(typename VecType::value_type)*v.capacity();
}

template<EntityRank TargetRank >
class BucketConnectivity<TargetRank, FIXED_CONNECTIVITY>
{
 public:
  typedef BucketConnectivity<TargetRank, FIXED_CONNECTIVITY> SelfType;

  static const EntityRank target_rank = TargetRank;
  static const ConnectivityType connectivity_type = FIXED_CONNECTIVITY;

  typedef std::vector<Entity> EntityVector;
  typedef std::vector<ConnectivityOrdinal> ConnectivityOrdinalVector;
  typedef std::vector<Permutation> PermutationVector;

  BucketConnectivity() //default constructed BucketConnectivity implies connectivity is not used
    : m_num_connectivity(0u)
    , m_targets()
    , m_ordinals()
    , m_permutations()
  {}

  BucketConnectivity(unsigned arg_num_connectivity)
    : m_num_connectivity(0)
    , m_targets()
    , m_ordinals()
    , m_permutations()
  {
    set_num_connectivity(arg_num_connectivity);
  }

  void set_num_connectivity(unsigned arg_num_connectivity)
  {
    STK_ThrowAssertMsg(m_num_connectivity == 0, "Cannot reset num_connectivity");
    STK_ThrowAssertMsg(arg_num_connectivity != 0, "Cannot set num connectivity to 0 for fixed connectivity");

    m_num_connectivity = arg_num_connectivity;

    // Ordinal is the same for all fixed, just counts up
    m_ordinals.resize(m_num_connectivity);
    for (ConnectivityOrdinal ord=static_cast<ConnectivityOrdinal>(0); ord < static_cast<ConnectivityOrdinal>(m_num_connectivity); ++ord) {
      m_ordinals[ord] = ord;
    }
  }

  // Entity iterator

  const ConnectedEntities get_connected_entities(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return ConnectedEntities(&m_targets[bucket_ordinal * m_num_connectivity], m_num_connectivity); }
  ConnectedEntities get_connected_entities(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return ConnectedEntities(&m_targets[bucket_ordinal * m_num_connectivity], m_num_connectivity); }

  Entity const* begin(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[bucket_ordinal * m_num_connectivity]; }

  Entity      * begin(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[bucket_ordinal * m_num_connectivity]; }

  Entity const* end(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[(bucket_ordinal + 1) * m_num_connectivity]; }

  Entity      * end(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return &m_targets[(bucket_ordinal + 1) * m_num_connectivity]; }

  // Ordinal iterator

  ConnectivityOrdinal const* begin_ordinals(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return m_ordinals.data(); }

  ConnectivityOrdinal      * begin_ordinals(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return m_ordinals.data(); }

  ConnectivityOrdinal const* end_ordinals(unsigned bucket_ordinal) const
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return m_ordinals.data() + m_num_connectivity; }

  ConnectivityOrdinal      * end_ordinals(unsigned bucket_ordinal)
  { impl::check_bucket_ordinal(bucket_ordinal, this);
    return m_ordinals.data() + m_num_connectivity; }

  // Permutation iterator

  Permutation const* begin_permutations(unsigned bucket_ordinal) const
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[bucket_ordinal * m_num_connectivity];
  }

  Permutation      * begin_permutations(unsigned bucket_ordinal)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[bucket_ordinal * m_num_connectivity];
  }

  Permutation const* end_permutations(unsigned bucket_ordinal) const
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[(bucket_ordinal + 1) * m_num_connectivity];
  }

  Permutation      * end_permutations(unsigned bucket_ordinal)
  {
    impl::check_bucket_ordinal(bucket_ordinal, this);
    if (!has_permutation()) return NULL;
    return &m_permutations[(bucket_ordinal + 1) * m_num_connectivity];
  }

  // Queries

  unsigned num_connectivity(unsigned /*bucket_ordinal*/) const
  { return m_num_connectivity; }

  // return number of entities
  unsigned size() const
  { return m_targets.size() / m_num_connectivity; }

  // Modification API

  bool add_connectivity(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal, Permutation permutation = INVALID_PERMUTATION)
  {
    STK_ThrowAssertMsg(ordinal < m_num_connectivity,
                   "Ordinal " <<  (uint32_t)ordinal << " exceeds topological limit: " << m_num_connectivity);
    impl::check_bucket_ordinal(bucket_ordinal, this);
#ifndef NDEBUG
    // TODO - Add topology invariant; target entity should match a sub topology
#endif

    unsigned index = m_num_connectivity*bucket_ordinal + ordinal;

    if (m_targets[index] == to) {
      STK_ThrowAssert(!has_permutation() || m_permutations[index] == permutation);
      // Already exists
      return false;
    }

    m_targets[index] = to;

    if (has_permutation()) {
      m_permutations[index] = permutation;
    }

    return true;
  }

  bool remove_connectivity(unsigned bucket_ordinal, Entity to, ConnectivityOrdinal ordinal)
  {
    STK_ThrowAssertMsg(ordinal < m_num_connectivity,
                   "Ordinal " <<  (uint32_t)ordinal << " exceeds topological limit: " << m_num_connectivity);
    impl::check_bucket_ordinal(bucket_ordinal, this);

    unsigned index = m_num_connectivity*bucket_ordinal + ordinal;
    if (m_targets[index] != to) {
      return false;
    }

    // Clear
    m_targets[index] = Entity();
    if (has_permutation()) {
      m_permutations[index] = INVALID_PERMUTATION;
    }

    return true;
  }

  bool replace_connectivity(unsigned bucket_ordinal, unsigned numConnectivity,
                            const Entity* connectivity,
                            const ConnectivityOrdinal* ordinals,
                            const Permutation* perms)
  {
    if (bucket_ordinal == size()) {
      add_entity();
    }

    impl::check_bucket_ordinal(bucket_ordinal, this);
    const unsigned index = m_num_connectivity*bucket_ordinal;
    for(unsigned i=0; i<numConnectivity; ++i) {
      m_targets[index+ordinals[i]] = connectivity[i];
      if (has_permutation()) {
        m_permutations[index+ordinals[i]] = perms[i];
      }
    }
    return true;
  }

  void begin_modification()
  {}

  template <typename BULKDATA> // hack to get around dependency
  void end_modification(BULKDATA* mesh = NULL);

  void add_entity()
  {
    const unsigned new_conn_size = m_targets.size() + m_num_connectivity;
    Entity invalid;
    m_targets.resize(new_conn_size, invalid); // Not a perf issue: vectors are smart when resizing
    if (has_permutation()) {
      m_permutations.resize(new_conn_size, INVALID_PERMUTATION);
    }
  }

  // Always removes last entity
  void remove_entity()
  {
    STK_ThrowAssertMsg(size() > 0, "Cannot remove, connectivity is already empty");

    const unsigned new_conn_size = m_targets.size() - m_num_connectivity;
    m_targets.resize(new_conn_size);
    if (has_permutation()) {
      m_permutations.resize(new_conn_size);
    }
  }

  void copy_entity(unsigned from_ordinal, SelfType& to, unsigned to_ordinal=-1u)
  {
    STK_ThrowAssertMsg(size() > 0, "Cannot move, connectivity is empty");

    if (to_ordinal == -1u) { // check if we should just append
      to_ordinal = to.size();
      to.add_entity(); // make room for new entity
    }

    // Copy connectivity to other BucketConnectivity
    copy_connectivity(from_ordinal, to, to_ordinal);
  }

  bool has_permutation() const
  {
    const static bool rv = TargetRank != stk::topology::NODE_RANK;
    return rv;
  }

  size_t heap_memory_in_bytes() const
  {
     return capacity_in_bytes(m_targets)
          + capacity_in_bytes(m_ordinals)
          + capacity_in_bytes(m_permutations);
  }

  void debug_dump(std::ostream& out) const
  {
    out << "For fixed connectivity to rank: " << TargetRank << "\n";
    out << "  size is: " << size() << "\n";
    for (int i = 0, ie = size(); i < ie; ++i) {
      out << "    At ordinal " << i << "\n";
      debug_dump(out, i);
    }
    out << std::endl;
  }

  void debug_dump(std::ostream& out, unsigned ordinal) const
  {
    for (int j = m_num_connectivity * ordinal, je = m_num_connectivity*(ordinal+1); j < je; ++j) {
      out << "        target:" << m_targets[j].local_offset() << "\n";
    }
  }

private:

  void copy_connectivity(unsigned from_ordinal, SelfType& to, unsigned to_ordinal)
  {
    unsigned to_offset   = to_ordinal * m_num_connectivity;
    unsigned from_offset = from_ordinal * m_num_connectivity;

    std::copy(m_targets.begin() + from_offset,
              m_targets.begin() + from_offset + m_num_connectivity,
              to.m_targets.begin() + to_offset);

    if (has_permutation()) {
      std::copy(m_permutations.begin() + from_offset,
                m_permutations.begin() + from_offset + m_num_connectivity,
                to.m_permutations.begin() + to_offset);
    }
  }

  // Illegal
  BucketConnectivity(const SelfType&);
  SelfType& operator=(const SelfType&);

  // MEMBERS

  unsigned m_num_connectivity;

  // connectivity data
  EntityVector              m_targets;
  ConnectivityOrdinalVector m_ordinals; // shared for all entities
  PermutationVector         m_permutations;
};

}

//
// BucketConnectivity
//

template <EntityRank TargetRank>
template <typename BULKDATA> // hack to get around dependency
inline
void impl::BucketConnectivity<TargetRank, FIXED_CONNECTIVITY>::end_modification(BULKDATA* mesh)
{
  //TODO: If bucket is blocked, no longer need to shrink to fit!

  if (m_targets.size() < m_targets.capacity()) {

    {
      EntityVector temp(m_targets.begin(), m_targets.end());
      m_targets.swap(temp);
    }

    {
      PermutationVector temp(m_permutations.begin(), m_permutations.end());
      m_permutations.swap(temp);
    }
  }

}

}} //namespace stk::mesh::impl

#endif
