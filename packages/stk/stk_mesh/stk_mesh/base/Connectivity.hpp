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

#ifndef STK_MESH_CONNECTIVITY_HPP
#define STK_MESH_CONNECTIVITY_HPP

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>      // for EntityRank
#include <stk_topology/topology.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <vector>

namespace stk {
namespace mesh {

class Connectivity {
public:
  Connectivity();

  static constexpr uint16_t MAX_NUM_RANKS = stk::topology::NUM_RANKS+1;

  Connectivity(const Connectivity& rhs)
   : m_connectivity(nullptr),
     m_offsetsAndOrdinals(nullptr)
  {
    *this = rhs;
  }

  Connectivity(Connectivity && rhs) noexcept
   : m_connectivity(nullptr),
     m_offsetsAndOrdinals(nullptr)
  {
    *this = std::move(rhs);
  }

  Connectivity& operator=(const Connectivity& rhs)
  {
    delete [] m_connectivity;
    delete [] m_offsetsAndOrdinals;

    unsigned capacity = rhs.get_capacity();
    if (capacity > 0) {
      m_connectivity = new Entity[capacity];
    }

    m_offsetsAndOrdinals = new ConnectivityOrdinal[m_firstOrdinalIndex + capacity];
    for(unsigned i=0; i<m_capacityIndex; ++i) {
      m_offsetsAndOrdinals[i] = rhs.m_offsetsAndOrdinals[i];
    }

    m_offsetsAndOrdinals[m_capacityIndex] = capacity;
    for(uint16_t i=0; i<m_offsetsAndOrdinals[MAX_NUM_RANKS]; ++i) {
      m_connectivity[i] = rhs.m_connectivity[i];
      m_offsetsAndOrdinals[m_capacityIndex+i] = rhs.m_offsetsAndOrdinals[m_capacityIndex+i];
    }

    return *this;
  }

  Connectivity& operator=(Connectivity&& rhs)
  {
    delete [] m_connectivity;
    delete [] m_offsetsAndOrdinals;

    m_connectivity = std::move(rhs.m_connectivity);
    m_offsetsAndOrdinals = std::move(rhs.m_offsetsAndOrdinals);
    rhs.m_connectivity = nullptr;
    rhs.m_offsetsAndOrdinals = nullptr;

    return *this;
  }

  ~Connectivity();

  inline
  unsigned num_connectivity(EntityRank connectedRank) const
  {
    ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
    return m_offsetsAndOrdinals[connectedRank+1]-m_offsetsAndOrdinals[connectedRank];
  }

  inline
  const Entity* begin_connectivity(EntityRank connectedRank) const
  {
    ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
    return m_connectivity+m_offsetsAndOrdinals[connectedRank];
  }

  inline
  const Entity* end_connectivity(EntityRank connectedRank) const
  {
    ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
    return m_connectivity+m_offsetsAndOrdinals[connectedRank+1];
  }

  inline
  PairIterEntity get_connectivity(EntityRank connectedRank) const
  {
    ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
    return PairIterEntity(m_connectivity+m_offsetsAndOrdinals[connectedRank], m_connectivity+m_offsetsAndOrdinals[connectedRank+1]);
  }

  inline
  PairIterEntity get_highest_rank_connectivity(EntityRank fromRank) const
  {
    ThrowAssertMsg(fromRank < MAX_NUM_RANKS, "fromRank="<<fromRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
    static constexpr int maxRank = MAX_NUM_RANKS-1;
    int highestRank = maxRank;
    const int self = fromRank;
    const ConnectivityOrdinal end = m_offsetsAndOrdinals[MAX_NUM_RANKS];
    while(highestRank > self && m_offsetsAndOrdinals[highestRank] == end) {
      --highestRank;
    }
    return PairIterEntity(m_connectivity+m_offsetsAndOrdinals[highestRank], m_connectivity+m_offsetsAndOrdinals[highestRank+1]);
  }

  inline
  PairIterEntity get_connectivity(EntityRank beginRank, EntityRank endRank) const
  {
    ThrowAssertMsg(endRank <= MAX_NUM_RANKS, "endRank="<<endRank<<" is required to be less-or-equal MAX_NUM_RANKS="<<MAX_NUM_RANKS);
    return PairIterEntity(m_connectivity+m_offsetsAndOrdinals[beginRank], m_connectivity+m_offsetsAndOrdinals[endRank]);
  }

  const ConnectivityOrdinal* begin_ordinals(EntityRank connectedRank) const
  {
    ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
    return m_offsetsAndOrdinals+m_firstOrdinalIndex+m_offsetsAndOrdinals[connectedRank];
  }

  const ConnectivityOrdinal* end_ordinals(EntityRank connectedRank) const
  {
    ThrowAssertMsg(connectedRank < MAX_NUM_RANKS, "connectedRank="<<connectedRank<<" is required to be less than MAX_NUM_RANKS="<<MAX_NUM_RANKS);
    return m_offsetsAndOrdinals+m_firstOrdinalIndex+m_offsetsAndOrdinals[connectedRank+1];
  }

  int add_connectivity(EntityRank connectedRank,
                       Entity connectedEntity,
                       ConnectivityOrdinal ordinal,
                       bool allowDuplicateOrdinals);

  int replace_or_add_connectivity(EntityRank connectedRank,
                                  Entity connectedEntity,
                                  ConnectivityOrdinal ordinal);

  int remove_connectivity(EntityRank connectedRank,
                          Entity connectedEntity,
                          ConnectivityOrdinal ordinal);

  void clear()
  {
    for(unsigned idx = 0; idx<=MAX_NUM_RANKS; ++idx) {
      m_offsetsAndOrdinals[idx] = 0;
    }
  }

protected:
  static constexpr unsigned m_capacityIndex = MAX_NUM_RANKS + 1;
  static constexpr unsigned m_firstOrdinalIndex = m_capacityIndex + 1;

  unsigned get_capacity() const
  { return m_offsetsAndOrdinals==nullptr ? 0 : m_offsetsAndOrdinals[m_capacityIndex]; }

  Entity* internal_begin_conn(EntityRank connectedRank)
  { return m_connectivity+m_offsetsAndOrdinals[connectedRank]; }

  inline
  Entity* internal_end_conn(EntityRank connectedRank)
  { return m_connectivity+m_offsetsAndOrdinals[connectedRank+1]; }

  ConnectivityOrdinal* internal_begin_ord(EntityRank connectedRank)
  { return m_offsetsAndOrdinals+m_firstOrdinalIndex+m_offsetsAndOrdinals[connectedRank]; }

  ConnectivityOrdinal* internal_end_ord(EntityRank connectedRank)
  { return m_offsetsAndOrdinals+m_firstOrdinalIndex+m_offsetsAndOrdinals[connectedRank+1]; }

  int find_insert_location(EntityRank connectedRank,
                           Entity connectedEntity,
                           ConnectivityOrdinal ordinal,
                           bool allowDuplicateOrdinals);

  int find_ordinal(EntityRank connectedRank,
                   ConnectivityOrdinal ordinal);

  void do_the_insert(unsigned insertOffset,
                     Entity newEntity,
                     ConnectivityOrdinal newOrdinal);

  void insert(EntityRank connectedRank,
              Entity connectedEntity,
              ConnectivityOrdinal ordinal,
              unsigned insertLocation);

  bool replace(EntityRank connectedRank,
               Entity connectedEntity,
               ConnectivityOrdinal ordinal,
               unsigned location);

  void remove(EntityRank connectedRank, unsigned location);

  Entity* m_connectivity;
  ConnectivityOrdinal* m_offsetsAndOrdinals;

  //
  //Ugly comment:
  //m_offsetsAndOrdinals is laid out like this:
  //  [0..MAX_NUM_RANKS] offsets into ordinal data and connectivity data for each entity-rank
  //  [m_capacityIndex] allocated capacity of m_connectivity and of ordinal data
  //  [m_firstOrdinalIndex...] ordinal data
  //
  //see begin_connectivity(rank) and begin_ordinals(rank) for indexing demonstration
  //
};

} // namespace mesh
} // namespace stk

#endif //STK_MESH_CONNECTIVITY_HPP

