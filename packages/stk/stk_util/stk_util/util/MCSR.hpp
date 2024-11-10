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

#ifndef stk_util_util_MCSR_hpp
#define stk_util_util_MCSR_hpp

//----------------------------------------------------------------------

#include <stk_util/stk_config.h>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/PairIter.hpp>
#include <algorithm>

//----------------------------------------------------------------------

namespace stk {
namespace util {

using IndexRange = std::pair<unsigned,unsigned>;

template<typename T>
int find_sorted_insertion_index(const std::vector<T>& items, const IndexRange& indices, const T& item)
{
  const T* begItems = items.data()+indices.first;
  const T* endItems = items.data()+indices.second;
  const T* it = std::lower_bound(begItems, endItems, item);
  if (it != endItems) {
    if (*it==item) {
      return -1;
    }
    return indices.first + (it - begItems);
  }
  return indices.second;
}

template<typename T>
class MCSR
{
public:
  MCSR(unsigned numRows, const T& invalidItem)
  : m_offsets(numRows, IndexRange(0u,0u)),
    m_items(),
    m_numUnusedEntries(0),
    m_compressionThreshold(0.5),
    m_invalidItem(invalidItem)
  {
  }

  ~MCSR() { }

  unsigned num_rows() const { return m_offsets.size(); }

  unsigned add_row()
  {
    unsigned newRow = m_offsets.size();
    m_offsets.push_back(IndexRange(0, 0));
    return newRow;
  }

  unsigned size(unsigned row) const
  {
    STK_ThrowAssertMsg(row < m_offsets.size(),"MCSR::size: row("<<row<<") must be less than num_rows()("<<m_offsets.size()<<")");

    const IndexRange& range = m_offsets[row];
    return range.second - range.first;
  }

  PairIter<const T*> items(unsigned row) const
  {
    STK_ThrowAssertMsg(row < m_offsets.size(),"MCSR::items: row("<<row<<") must be less than num_rows()("<<m_offsets.size()<<")");
    const IndexRange& indices = m_offsets[row];
    return PairIter<const T*>(&m_items[indices.first], &m_items[indices.second]);
  }

  const T* begin(unsigned row) const
  {
    STK_ThrowAssertMsg(row < m_offsets.size(),"MCSR::begin: row("<<row<<") must be less than num_rows()("<<m_offsets.size()<<")");

    return &m_items[m_offsets[row].first];
  }

  const T* end(unsigned row) const
  {
    STK_ThrowAssertMsg(row < m_offsets.size(),"MCSR::end: row("<<row<<") must be less than num_rows()("<<m_offsets.size()<<")");

    return &m_items[m_offsets[row].second];
  }

  bool add_item(unsigned row, const T& item)
  {
    return insert_item(row, item);
  }

  bool remove_item(unsigned row, const T& item)
  {
    IndexRange& indices = m_offsets[row];
    unsigned idx = indices.second;
    for(unsigned i=indices.first; i<indices.second; ++i) {
      if (m_items[i] == item) {
        idx = i;
        break;
      }
    }

    if (idx < indices.second) {
      for(unsigned i=idx; i<indices.second - 1; ++i) {
        m_items[i] = m_items[i+1];
      }

      --indices.second;
      m_items[indices.second] = m_invalidItem;
      ++m_numUnusedEntries;
      return true;
    }

    return false;
  }

  bool remove_items(unsigned row)
  {
    IndexRange& indices = m_offsets[row];
    for(unsigned i=indices.first; i<indices.second; ++i) {
      m_items[i] = m_invalidItem;
    }

    unsigned numRemoved = indices.second - indices.first;
    indices.second = indices.first;
    m_numUnusedEntries += numRemoved;
    return numRemoved > 0;
  }

  template<class Matcher>
  bool remove_items_if(unsigned row, const Matcher& matcher)
  {
    IndexRange& indices = m_offsets[row];
    unsigned numRemoved = 0;
    for(unsigned i=indices.first; i<indices.second; ++i) {
      if (matcher(m_items[i])) {
        m_items[i] = m_invalidItem;
        ++numRemoved;
      }
    }

    if (numRemoved > 0) {
      if (numRemoved < (indices.second-indices.first)) {
        unsigned idx = indices.first;
        for(unsigned i=indices.first; i<indices.second; ++i) {
          if (m_items[i] != m_invalidItem) {
            if (i > idx) {
              m_items[idx] = m_items[i];
            }
            ++idx;
          }
        }
      }
      indices.second -= numRemoved;
      m_numUnusedEntries += numRemoved;
      return true;
    }
    return false;
  }

  bool remove_items_greater_equal(unsigned row, T item)
  {
    IndexRange& indices = m_offsets[row];
    unsigned numRemoved = 0;
    for(unsigned i=indices.first; i<indices.second; ++i) {
      if (m_items[i] >= item) {
        m_items[i] = m_invalidItem;
        ++numRemoved;
      }
    }

    indices.second -= numRemoved;
    m_numUnusedEntries += numRemoved;
    return true;
  }

  size_t total_capacity() const { return m_items.capacity(); }
  size_t total_num_items() const { return m_items.size() - m_numUnusedEntries; }
  size_t num_unused_entries() const { return m_numUnusedEntries; }

  void compress(unsigned suggestedCapacity = 0)
  {
    if (m_numUnusedEntries == 0) {
      return;
    }

    size_t initialCapacity = suggestedCapacity == 0 ? total_num_items() : suggestedCapacity;
    std::vector<T> m_newItems;
    m_newItems.reserve(initialCapacity);

    for(IndexRange& indices : m_offsets) {
      unsigned newStartIdx = m_newItems.size();
      for(unsigned i=indices.first; i<indices.second; ++i) {
        m_newItems.push_back(m_items[i]);
      }
      indices.first = newStartIdx;
      indices.second = m_newItems.size();
    }

    m_items.swap(m_newItems);

    m_numUnusedEntries = 0;
  }

  void clear(unsigned initialCapacity, std::vector<T>& itemsToSwap)
  {
    std::vector<IndexRange> m_newOffsets;
    m_newOffsets.reserve(initialCapacity);
    itemsToSwap.clear();
    itemsToSwap.reserve(initialCapacity);
    m_offsets.swap(m_newOffsets);
    m_items.swap(itemsToSwap);
    m_numUnusedEntries = 0;
  }

private:
  bool is_valid(const T& item)
  {
    return item != m_invalidItem;
  }

  bool insert_item(unsigned row, const T& item)
  {
    size_t numItems = m_items.size();
    if (numItems > 256u && (static_cast<double>(m_numUnusedEntries)/numItems) > m_compressionThreshold)
    {
      compress(numItems-m_numUnusedEntries/2);
      numItems = m_items.size();
    }

    IndexRange& indices = m_offsets[row];

    if (indices.second >= numItems) {
      int insertIdx = find_sorted_insertion_index(m_items, indices, item);
      if (insertIdx < 0) {
        return false;
      }
      m_items.emplace_back();
      return insert_item_at_index(indices, insertIdx, item);
    }

    if (!is_valid(m_items[indices.second])) {
      bool didInsert = insert_item_into_sorted_range(indices, item);
      if (didInsert) {
        --m_numUnusedEntries;
      }
      return didInsert;
    }
    else {
      const int insertIdx = find_sorted_insertion_index(m_items, indices, item);
      if (insertIdx < 0) {
        return false;
      }

      const unsigned distanceMoved = move_items_to_end(row);
      m_items.emplace_back();

      return insert_item_at_index(indices, (insertIdx+distanceMoved), item);
    }

    return false;
  }

  bool insert_item_at_index(IndexRange& indices, int insertIdx, const T& item)
  {
    if (insertIdx < 0) {
      return false;
    }
    unsigned uinsertIdx = static_cast<unsigned>(insertIdx);
    for(unsigned i=indices.second; i>uinsertIdx; --i) {
      m_items[i] = m_items[i-1];
    }

    m_items[insertIdx] = item;

    ++indices.second;
    return true;
  }
  
  bool insert_item_into_sorted_range(IndexRange& indices, const T& item)
  {
    return insert_item_at_index(indices, find_sorted_insertion_index(m_items, indices, item), item);
  }
  
  unsigned move_items_to_end(unsigned row)
  {
    IndexRange& indices = m_offsets[row];
    const unsigned newStartIdx = m_items.size();

    for(unsigned i=indices.first; i<indices.second; ++i) {
      m_items.push_back(m_items[i]);
      m_items[i] = m_invalidItem;
    }

    m_numUnusedEntries += indices.second - indices.first;
    const unsigned distanceMoved = newStartIdx - indices.first;
    m_offsets[row] = IndexRange(newStartIdx, m_items.size());
    return distanceMoved;
  }

  std::vector<IndexRange> m_offsets;
  std::vector<T> m_items;
  unsigned m_numUnusedEntries;
  double m_compressionThreshold;
  T m_invalidItem;
};

} // namespace util
} // namespace stk

#endif

