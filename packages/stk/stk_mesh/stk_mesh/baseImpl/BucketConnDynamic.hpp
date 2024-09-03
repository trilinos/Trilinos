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

#ifndef stk_mesh_impl_BucketConnDynamic_hpp
#define stk_mesh_impl_BucketConnDynamic_hpp

//----------------------------------------------------------------------

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/util/StridedArray.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <algorithm>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {

using ConnectedEntities = util::StridedArray<const stk::mesh::Entity>;

namespace impl {

class BucketConnDynamic
{
public:
  BucketConnDynamic(unsigned bucketCapacity, bool hasPermutations = false)
  : m_bucketCapacity(bucketCapacity),
    m_hasPermutations(hasPermutations),
    m_offsets(),
    m_connectivity(),
    m_ordinals(),
    m_permutations(),
    m_numUnusedEntries(0),
    m_compressionThreshold(0.5)
  {
    STK_ThrowRequireMsg(bucketCapacity > 0, "BucketConnDynamic must have bucketCapacity strictly greater than 0");
  }

  ~BucketConnDynamic() { }

  bool has_permutation() const { return m_hasPermutations; }

  unsigned num_connectivity(unsigned bktOrdinal) const
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::num_connectivity: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");

    auto [first, second] = m_offsets[bktOrdinal];
    return second - first;
  }

  const ConnectedEntities get_connected_entities(unsigned bktOrdinal) const
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::get_connected_entities: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");
    auto [first, second] = m_connectivity.empty() ? IndexRange{0, 0} : m_offsets[bktOrdinal];
    const unsigned len = second - first;
    const Entity* ptr = len==0 ? nullptr : m_connectivity.data()+first;
    return ConnectedEntities(ptr, len);
  }

  const Entity* begin(unsigned bktOrdinal) const
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::begin: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");

    return m_connectivity.empty() ? nullptr : m_connectivity.data()+m_offsets[bktOrdinal].first;
  }
  Entity* begin(unsigned bktOrdinal)
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::begin: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");

    return m_connectivity.empty() ? nullptr : m_connectivity.data()+m_offsets[bktOrdinal].first;
  }

  const Entity* end(unsigned bktOrdinal) const
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::end: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");

    return m_connectivity.empty() ? nullptr : m_connectivity.data()+m_offsets[bktOrdinal].second;
  }
  Entity* end(unsigned bktOrdinal)
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::end: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");

    return m_connectivity.empty() ? nullptr : m_connectivity.data()+m_offsets[bktOrdinal].second;
  }

  const ConnectivityOrdinal* begin_ordinals(unsigned bktOrdinal) const
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::begin_ordinals: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");

    return m_ordinals.empty() ? nullptr : m_ordinals.data()+m_offsets[bktOrdinal].first;
  }
  const ConnectivityOrdinal* end_ordinals(unsigned bktOrdinal) const
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::end_ordinals: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");

    return m_ordinals.empty() ? nullptr : m_ordinals.data()+m_offsets[bktOrdinal].second;
  }

  const Permutation* begin_permutations(unsigned bktOrdinal) const
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::begin_permutations: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");
    return (has_permutation() && !m_permutations.empty()) ?
       m_permutations.data()+m_offsets[bktOrdinal].first : nullptr;
  }
  Permutation* begin_permutations(unsigned bktOrdinal)
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::begin_permutations: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");
    return (has_permutation() && !m_permutations.empty()) ?
       m_permutations.data()+m_offsets[bktOrdinal].first : nullptr;
  }

  const Permutation* end_permutations(unsigned bktOrdinal) const
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::end_permutations: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");
    return (has_permutation() && !m_permutations.empty()) ?
       m_permutations.data()+m_offsets[bktOrdinal].second : nullptr;
  }
  Permutation* end_permutations(unsigned bktOrdinal)
  {
    STK_ThrowAssertMsg(bktOrdinal < m_bucketCapacity,"BucketConnDynamic::end_permutations: bktOrdinal("<<bktOrdinal<<") must be less than bucketCapacity("<<m_bucketCapacity<<")");
    return (has_permutation() && !m_permutations.empty()) ?
       m_permutations.data()+m_offsets[bktOrdinal].second : nullptr;
  }

  bool add_connectivity(unsigned bktOrdinal,
                        Entity entity,
                        ConnectivityOrdinal ordinal,
                        Permutation perm = INVALID_PERMUTATION)
  {
    return insert_connectivity(bktOrdinal, entity, ordinal, perm);
  }

  bool remove_connectivity(unsigned bktOrdinal,
                           Entity entity,
                           ConnectivityOrdinal ordinal,
                           Permutation perm = INVALID_PERMUTATION)
  {
    IndexRange& indices = m_offsets[bktOrdinal];
    UpwardConnIndexType idx = indices.second;
    for(UpwardConnIndexType i=indices.first; i<indices.second; ++i) {
      if (m_connectivity[i] == entity && m_ordinals[i] == ordinal) {
        idx = i;
        break;
      }
    }

    if (idx < indices.second) {
      for(UpwardConnIndexType i=idx; i<indices.second - 1; ++i) {
        m_connectivity[i] = m_connectivity[i+1];
        m_ordinals[i] = m_ordinals[i+1];
        if (m_hasPermutations) {
          m_permutations[i] = m_permutations[i+1];
        }
      }

      --indices.second;
      m_ordinals[indices.second] = INVALID_CONNECTIVITY_ORDINAL;
      ++m_numUnusedEntries;
      return true;
    }

    return false;
  }

  bool remove_connectivity(unsigned bktOrdinal)
  {
    IndexRange& indices = m_offsets[bktOrdinal];
    for(UpwardConnIndexType i=indices.first; i<indices.second; ++i) {
      m_ordinals[i] = INVALID_CONNECTIVITY_ORDINAL;
    }

    UpwardConnIndexType numRemoved = indices.second - indices.first;
    indices.second = indices.first;
    m_numUnusedEntries += numRemoved;
    return true;
  }

  bool replace_connectivity(unsigned bktOrdinal,
                            unsigned numConnectivity,
                            const Entity* connectivity,
                            const ConnectivityOrdinal* ordinals,
                            const Permutation* permutations)
  {
    IndexRange& indices = m_offsets[bktOrdinal];
    const unsigned numExisting = num_connectivity(bktOrdinal);
    const unsigned numToReplace = std::min(numConnectivity, numExisting);
    for(unsigned i=0; i<numToReplace; ++i) {
      m_connectivity[i+indices.first] = connectivity[i];
      m_ordinals[i+indices.first] = ordinals[i];
      if (m_hasPermutations) {
        m_permutations[i+indices.first] = permutations[i];
      }
    }
    
    if (numToReplace < numExisting) {
      for(unsigned i=numToReplace; i<numExisting; ++i) {
        m_ordinals[i+indices.first] = INVALID_CONNECTIVITY_ORDINAL;
      }
      indices.second = indices.first + numToReplace;
      m_numUnusedEntries += numExisting - numToReplace;
      return true;
    }
    else {
      for(unsigned i=numExisting; i<numConnectivity; ++i) {
        Permutation perm = m_hasPermutations ? permutations[i] : INVALID_PERMUTATION;
        add_connectivity(bktOrdinal, connectivity[i], ordinals[i], perm);
      }
    }

    return true;
  }

  bool swap_connectivity(unsigned bktOrdinal1, unsigned bktOrdinal2)
  {
    IndexRange tmp = m_offsets[bktOrdinal1];
    m_offsets[bktOrdinal1] = m_offsets[bktOrdinal2];
    m_offsets[bktOrdinal2] = tmp;
    return true;
  }

  size_t bucket_size() const { return m_offsets.size(); }
  size_t bucket_capacity() const { return m_bucketCapacity; }
  size_t total_capacity() const { return m_connectivity.capacity(); }
  size_t total_num_connectivity() const { return m_connectivity.size() - m_numUnusedEntries; }
  size_t num_unused_entries() const { return m_numUnusedEntries; }

  void compress_connectivity(unsigned suggestedCapacity = 0)
  {
    if (m_numUnusedEntries == 0) {
      return;
    }

    std::vector<std::pair<IndexRange,unsigned>> sortedOffsets(m_offsets.size());
    for(unsigned i=0; i<m_offsets.size(); ++i) {
      sortedOffsets[i].first = m_offsets[i];
      sortedOffsets[i].second = i;
    }

    std::sort(sortedOffsets.begin(), sortedOffsets.end());

    if (sortedOffsets[0].first.first > 0) {
      IndexRange& sRange = sortedOffsets[0].first;
      const unsigned gap = sRange.first;
      slide_range_and_update(sRange, gap, sortedOffsets[0].second);
    }

    for(unsigned i=0; i<sortedOffsets.size()-1; ++i) {
      const unsigned thisRangeEnd = sortedOffsets[i].first.second;
      const unsigned nextRangeBegin = sortedOffsets[i+1].first.first;
      const unsigned gap = nextRangeBegin - thisRangeEnd;
      if (gap > 0) {
        slide_range_and_update(sortedOffsets[i+1].first, gap, sortedOffsets[i+1].second);
      }
    }

    const unsigned oldSize = m_connectivity.size();
    m_connectivity.resize(oldSize-m_numUnusedEntries);
    m_ordinals.resize(oldSize-m_numUnusedEntries);
    if (m_hasPermutations) {
      m_permutations.resize(oldSize-m_numUnusedEntries);
    }
    m_numUnusedEntries = 0;
    const unsigned lastIdx = sortedOffsets.size()-1;
    STK_ThrowRequireMsg(sortedOffsets[lastIdx].first.second == m_connectivity.size(),
                        "Internal BucketConnDynamic::compress_connectivity ERROR, indices out of sync with data.");
  }

  void grow_if_necessary(unsigned bktOrdinal)
  {
    m_bucketCapacity = std::max(bktOrdinal+1, m_bucketCapacity);
    if (bktOrdinal >= m_offsets.size()) {
      const unsigned candidate = m_offsets.empty() ? bktOrdinal+1 : 2*m_offsets.size();
      const unsigned newSize = std::min(m_bucketCapacity, candidate);
      m_offsets.resize(newSize, IndexRange(0u, 0u));

      if (m_offsets.capacity() > m_bucketCapacity) {
        std::vector<IndexRange>(m_offsets).swap(m_offsets);
      }
    }
  }

  void increase_bucket_capacity(unsigned newBucketCapacity)
  {
    STK_ThrowRequireMsg(newBucketCapacity >= m_bucketCapacity, "BucketDynamicConn::increase_bucket_capacity, old capacity="<<m_bucketCapacity<<" should be less than new capacity="<<newBucketCapacity);
    m_bucketCapacity = newBucketCapacity;
  }

  size_t heap_memory_in_bytes() const
  {
    return sizeof(IndexRange) * m_offsets.capacity()
         + sizeof(Entity) * m_connectivity.capacity()
         + sizeof(ConnectivityOrdinal) * m_ordinals.capacity()
         + sizeof(Permutation) * m_permutations.capacity();
  }

private:
  using IndexRange = std::pair<UpwardConnIndexType,UpwardConnIndexType>;

  bool is_valid(ConnectivityOrdinal ordinal)
  {
    return ordinal != INVALID_CONNECTIVITY_ORDINAL;
  }

  void slide_range_and_update(IndexRange& range, unsigned gap, unsigned rangeOrd)
  {
    for(unsigned idx=range.first; idx<range.second; ++idx) {
      m_connectivity[idx-gap] = m_connectivity[idx];
      m_ordinals[idx-gap] = m_ordinals[idx];
      if (m_hasPermutations) {
        m_permutations[idx-gap] = m_permutations[idx];
      }
    }
    range.first -= gap;
    range.second -= gap;
    m_offsets[rangeOrd].first -= gap;
    m_offsets[rangeOrd].second -= gap;
  }

  bool insert_connectivity(unsigned bktOrdinal,
                           Entity entity,
                           ConnectivityOrdinal ordinal,
                           Permutation perm = INVALID_PERMUTATION)
  {
    static constexpr unsigned minSizeHeuristic = 256;
    if (total_num_connectivity() > minSizeHeuristic && (static_cast<double>(m_numUnusedEntries)/total_num_connectivity()) > m_compressionThreshold)
    {
      compress_connectivity(total_num_connectivity()+m_numUnusedEntries/2);
    }

    grow_if_necessary(bktOrdinal);

    IndexRange& indices = m_offsets[bktOrdinal];

    if (indices.second >= m_connectivity.size()) {
      int insertIdx = find_sorted_insertion_index(indices, entity, ordinal);
      if (insertIdx < 0) {
        return false;
      }
      m_connectivity.emplace_back();
      m_ordinals.emplace_back();
      if (m_hasPermutations) {
        m_permutations.emplace_back();
      }

      STK_ThrowRequireMsg(m_connectivity.size() < INVALID_UPWARDCONN_INDEX,"Internal error, BucketConnDynamic size exceeds limitation of index type");

      return insert_connectivity_at_idx(indices, insertIdx, entity, ordinal, perm);
    }

    if (!is_valid(m_ordinals[indices.second])) {
      bool didInsert = insert_connectivity_into_sorted_range(indices, entity, ordinal, perm);
      if (didInsert) {
        --m_numUnusedEntries;
      }
      return didInsert;
    }
    else {
      const int insertIdx = find_sorted_insertion_index(indices, entity, ordinal);
      if (insertIdx < 0) {
        return false;
      }

      const UpwardConnIndexType distanceMoved = move_connectivity_to_end(bktOrdinal);
      m_connectivity.emplace_back();
      m_ordinals.emplace_back();
      if (m_hasPermutations) {
        m_permutations.emplace_back();
      }

      STK_ThrowRequireMsg(m_connectivity.size() < INVALID_UPWARDCONN_INDEX,"Internal error, BucketConnDynamic size exceeds limitation of index type");

      return insert_connectivity_at_idx(indices, (insertIdx+distanceMoved), entity, ordinal, perm);
    }

    return false;
  }

  int find_sorted_insertion_index(const IndexRange& indices,
                                  Entity entity,
                                  ConnectivityOrdinal ordinal)
  {
    const ConnectivityOrdinal* beg = m_ordinals.data()+indices.first;
    const ConnectivityOrdinal* end = m_ordinals.data()+indices.second;
    const ConnectivityOrdinal* it = std::lower_bound(beg, end, ordinal);
    if (it != end) {
      int idx = indices.first + (it - beg);
      if (*it==ordinal) {
        for(; static_cast<unsigned>(idx)<indices.second && m_ordinals[idx]==ordinal; ++idx) {
          if (m_connectivity[idx] > entity) {
            return idx;
          }
          if (m_connectivity[idx] == entity) {
            return -1;
          }
        }
      }
      return idx;
    }
    return indices.second;
  }

  bool insert_connectivity_at_idx(IndexRange& indices,
                                  int insertIdx,
                                  Entity entity,
                                  ConnectivityOrdinal ordinal,
                                  Permutation perm)
  {
    if (insertIdx < 0) {
      return false;
    }
    unsigned uinsertIdx = static_cast<unsigned>(insertIdx);
    for(UpwardConnIndexType i=indices.second; i>uinsertIdx; --i) {
      m_connectivity[i] = m_connectivity[i-1];
      m_ordinals[i] = m_ordinals[i-1];
      if (m_hasPermutations) {
        m_permutations[i] = m_permutations[i-1];
      }
    }

    m_connectivity[insertIdx] = entity;
    m_ordinals[insertIdx] = ordinal;
    if (m_hasPermutations) {
      m_permutations[insertIdx] = perm;
    }

    ++indices.second;
    return true;
  }
  
  bool insert_connectivity_into_sorted_range(IndexRange& indices,
                                             Entity entity,
                                             ConnectivityOrdinal ordinal,
                                             Permutation perm)
  {
    return insert_connectivity_at_idx(indices,
                                      find_sorted_insertion_index(indices, entity, ordinal),
                                      entity, ordinal, perm);
  }
  
  UpwardConnIndexType move_connectivity_to_end(unsigned bktOrdinal)
  {
    IndexRange& indices = m_offsets[bktOrdinal];
    const UpwardConnIndexType newStartIdx = m_connectivity.size();

    for(UpwardConnIndexType i=indices.first; i<indices.second; ++i) {
      m_connectivity.push_back(m_connectivity[i]);
      m_ordinals.push_back(m_ordinals[i]);
      m_ordinals[i] = INVALID_CONNECTIVITY_ORDINAL;
      if (m_hasPermutations) {
        m_permutations.push_back(m_permutations[i]);
      }
    }

    STK_ThrowRequireMsg(m_connectivity.size() < INVALID_UPWARDCONN_INDEX,"Internal error, BucketConnDynamic size exceeds limitation of index type");

    m_numUnusedEntries += indices.second - indices.first;
    const UpwardConnIndexType distanceMoved = newStartIdx - indices.first;
    m_offsets[bktOrdinal] = IndexRange(newStartIdx, m_connectivity.size());
    return distanceMoved;
  }

  unsigned m_bucketCapacity;
  bool m_hasPermutations;

  std::vector<IndexRange> m_offsets;
  std::vector<Entity> m_connectivity;
  std::vector<ConnectivityOrdinal> m_ordinals;
  std::vector<Permutation> m_permutations;
  unsigned m_numUnusedEntries;
  double m_compressionThreshold;
};

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

