/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_UNORDERED_MAP2_IMPL_HPP
#define KOKKOS_UNORDERED_MAP2_IMPL_HPP

#include <Kokkos_Macros.hpp>
#include <stdint.h>

#include <cstdio>
#include <climits>
#include <iostream>
#include <iomanip>

namespace Kokkos { namespace Impl {

uint32_t find_hash_size( uint32_t size );

KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_forward(uint32_t i)
{
#if defined( __CUDA_ARCH__ )
  return __ffs(i) - 1;
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_ffs(i) - 1;
#elif defined( __INTEL_COMPILER )
  return _bit_scan_forward(i);
#else

  uint32_t t = 1;
  int r = 0;
  while (i && (i & t == 0))
  {
    t = t << 1;
    ++r;
  }
  return r;
#endif
}

KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_forward(uint64_t i)
{
#if defined( __CUDA_ARCH__ )
  return __ffsll(i) - 1;
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_ffsll(i) - 1;
#else
  uint64_t t = 1;
  int r = 0;
  while (i && (i & t == 0))
  {
    t = t << 1;
    ++r;
  }
  return r;
#endif
}

KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_reverse(uint32_t i)
{
#if defined( __CUDA_ARCH__ )
  return 31 - __clz(i);
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return 31 - __builtin_clz(i);
#elif defined( __INTEL_COMPILER )
  return _bit_scan_reverse(i);
#else
  uint32_t t = 1 << 31;
  int r = 0;
  while (i && (i & t == 0))
  {
    t = t >> 1;
    ++r;
  }
  return r;
#endif
}

KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_reverse(uint64_t i)
{
#if defined( __CUDA_ARCH__ )
  return 63 - __clzll(i);
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return 63 - __builtin_clzll(i);
#else
  uint64_t t = 1 << 63;
  int r = 0;
  while (i && (i & t == 0))
  {
    t = t >> 1;
    ++r;
  }
  return r;
#endif
}

// count the bits set
// http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
KOKKOS_FORCEINLINE_FUNCTION
int popcount(uint32_t i)
{
#if defined( __CUDA_ARCH__ )
  return __popc(i);
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_popcount(i);
#elif defined ( __INTEL_COMPILER )
  return _popcnt32(i);
#else
  i = i - ((i >> 1) & (uint32_t)~(uint32_t)0/3);                                     // temp
  i = (i & (uint32_t)~(uint32_t)0/15*3) + ((i >> 2) & (uint32_t)~(uint32_t)0/15*3);  // temp
  i = (i + (i >> 4)) & (uint32_t)~(uint32_t)0/255*15;                                // temp
  return (int)((uint32_t)(i * ((uint32_t)~(uint32_t)0/255)) >> (sizeof(uint32_t) - 1) * CHAR_BIT); // count
#endif
}

KOKKOS_FORCEINLINE_FUNCTION
int popcount(uint64_t i)
{
#if defined( __CUDA_ARCH__ )
  return __popcll(i);
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_popcountll(i);
#elif defined ( __INTEL_COMPILER )
  return _popcnt64(i);
#else
  return popcount( static_cast<uint32_t>(i >> 32) ) + popcount( static_cast<uint32_t>(i & ( ~0u )) );
#endif
}

struct UnorderedMapScalars
{
  bool modified;
  bool erasable;
  bool has_failed_inserts;
  uint32_t size;
  uint32_t failed_inserts;
};



template <typename Map>
struct UnorderedMapRehash
{
  typedef Map map_type;
  typedef typename map_type::const_map_type const_map_type;
  typedef typename map_type::device_type device_type;
  typedef typename map_type::size_type size_type;

  map_type       m_dst;
  const_map_type m_src;

  UnorderedMapRehash( map_type const& dst, const_map_type const& src)
    : m_dst(dst), m_src(src)
  {}

  void apply() const
  {
    parallel_for(m_src.capacity(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const
  {
    if ( m_src.valid_at(i) )
      m_dst.insert(m_src.key_at(i), m_src.value_at(i));
  }

};

template <typename UMap>
struct UnorderedMapSize
{
  typedef UMap map_type;
  typedef typename map_type::device_type device_type;
  typedef typename map_type::size_type size_type;
  typedef uint32_t value_type;

  map_type m_map;

  UnorderedMapSize( map_type const& map)
    : m_map(map)
  {}

  void apply() const
  {
    parallel_reduce(m_map.m_available_indexes.size(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & size)
  {
    size = 0;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & size, const volatile size_type & incr )
  {
    size += incr;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & size) const
  {
    size += popcount(~m_map.m_available_indexes[i]);
  }

  KOKKOS_INLINE_FUNCTION
  void final( value_type & size ) const
  {
    m_map.m_scalars().size = size;
  }
};

template <typename UMap>
struct UnorderedMapErase
{
  typedef UMap map_type;
  typedef typename map_type::device_type device_type;
  typedef typename map_type::size_type size_type;
  typedef typename map_type::key_type key_type;
  typedef typename map_type::impl_value_type value_type;

  map_type m_map;

  UnorderedMapErase( map_type const& map)
    : m_map(map)
  {}

  void apply() const
  {
    parallel_for(m_map.m_hash_lists.size(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    const size_type invalid_index = map_type::invalid_index;

    size_type curr = m_map.m_hash_lists(i);
    size_type next = invalid_index;

    // remove erased head of the linked-list
    while (curr != invalid_index && !m_map.valid_at(curr)) {
      next = m_map.m_next_index[curr];
      m_map.m_next_index[curr] = invalid_index;
      m_map.m_keys[curr] = key_type();
      if (m_map.is_set) m_map.m_values[curr] = value_type();
      curr = next;
      m_map.m_hash_lists(i) = next;
    }

    // if the list is non-empty and the head is valid
    if (curr != invalid_index && m_map.valid_at(curr) ) {
      size_type prev = curr;
      curr = m_map.m_next_index[prev];

      while (curr != invalid_index) {
        next = m_map.m_next_index[curr];
        if (m_map.valid_at(curr)) {
          prev = curr;
        }
        else {
          // remove curr from list
          m_map.m_next_index[prev] = next;
          m_map.m_next_index[curr] = invalid_index;
          m_map.m_keys[curr] = key_type();
          if (map_type::is_set) m_map.m_values[curr] = value_type();
        }
        curr = next;
      }
    }
  }
};

template <typename UMap>
struct UnorderedMapHistogram
{
  typedef UMap map_type;
  typedef typename map_type::device_type device_type;
  typedef typename map_type::size_type size_type;

  typedef View<int[100], device_type> histogram_view;
  typedef typename histogram_view::HostMirror host_histogram_view;

  map_type m_map;
  histogram_view m_length;
  histogram_view m_distance;
  histogram_view m_block_distance;

  UnorderedMapHistogram( map_type const& map)
    : m_map(map)
    , m_length("UnorderedMap Histogram")
    , m_distance("UnorderedMap Histogram")
    , m_block_distance("UnorderedMap Histogram")
  {}

  void calculate()
  {
    parallel_for(m_map.m_hash_lists.size(), *this);
  }

  void clear()
  {
    Kokkos::deep_copy(m_length, 0);
    Kokkos::deep_copy(m_distance, 0);
    Kokkos::deep_copy(m_block_distance, 0);
  }

  void print_length(std::ostream &out)
  {
    host_histogram_view host_copy = create_mirror_view(m_length);
    Kokkos::deep_copy(host_copy, m_length);

    for (int i=0, size = host_copy.size(); i<size; ++i)
    {
      out << host_copy[i] << " , ";
    }
    out << "\b\b\b   " << std::endl;
  }

  void print_distance(std::ostream &out)
  {
    host_histogram_view host_copy = create_mirror_view(m_distance);
    Kokkos::deep_copy(host_copy, m_distance);

    for (int i=0, size = host_copy.size(); i<size; ++i)
    {
      out << host_copy[i] << " , ";
    }
    out << "\b\b\b   " << std::endl;
  }

  void print_block_distance(std::ostream &out)
  {
    host_histogram_view host_copy = create_mirror_view(m_block_distance);
    Kokkos::deep_copy(host_copy, m_block_distance);

    for (int i=0, size = host_copy.size(); i<size; ++i)
    {
      out << host_copy[i] << " , ";
    }
    out << "\b\b\b   " << std::endl;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i ) const
  {
    const size_type invalid_index = map_type::invalid_index;

    uint32_t length = 0;
    size_type min_index = ~0u, max_index = 0;
    for (size_type curr = m_map.m_hash_lists(i); curr != invalid_index; curr = m_map.m_next_index[curr]) {
      ++length;
      min_index = (curr < min_index) ? curr : min_index;
      max_index = (max_index < curr) ? curr : max_index;
    }

    size_type distance = (0u < length) ? max_index - min_index : 0u;
    size_type blocks = (0u < length) ? max_index/map_type::block_size - min_index/map_type::block_size : 0u;

    // normalize data
    length   = length   < 100u ? length   : 99u;
    distance = distance < 100u ? distance : 99u;
    blocks   = blocks   < 100u ? blocks   : 99u;

    if (0u < length)
    {
      atomic_fetch_add( &m_length(length), 1);
      atomic_fetch_add( &m_distance(distance), 1);
      atomic_fetch_add( &m_block_distance(blocks), 1);
    }
  }
};

template <typename DKey, typename DValue, typename SKey, typename SValue>
struct UnorderedMapCanAssign : public false_ {};

template <typename Key, typename Value>
struct UnorderedMapCanAssign<Key,Value,Key,Value> : public true_ {};

template <typename Key, typename Value>
struct UnorderedMapCanAssign<const Key,Value,Key,Value> : public true_ {};

template <typename Key, typename Value>
struct UnorderedMapCanAssign<const Key,const Value,Key,Value> : public true_ {};

template <typename Key, typename Value>
struct UnorderedMapCanAssign<const Key,const Value,const Key,Value> : public true_ {};


}} //Kokkos::Impl

#endif // KOKKOS_UNORDERED_MAP2_IMPL_HPP
