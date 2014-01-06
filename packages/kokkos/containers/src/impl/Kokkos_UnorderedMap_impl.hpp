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

namespace Kokkos { namespace Impl {

uint32_t find_hash_size( uint32_t size );

// find the first set bit of the integer
KOKKOS_FORCEINLINE_FUNCTION
int find_first_set(uint32_t i)
{
#if defined( __CUDA_ARCH__ )
  return __ffs(i);
#elif defined( __INTEL_COMPILER )
  return i ? _bit_scan_forward(i) + 1 : 0;
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_ffs(i);
#else
  if(!i)
    return 0;
  uint32_t t = 1;
  int r = 1;
  while (i & t == 0)
  {
    t <<= 1;
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
  i = i - ((i >> 1) & (uint32_t)~(uint32_t)0/3);                                     // temp
  i = (i & (uint32_t)~(uint32_t)0/15*3) + ((i >> 2) & (uint32_t)~(uint32_t)0/15*3);  // temp
  i = (i + (i >> 4)) & (uint32_t)~(uint32_t)0/255*15;                                // temp
  return (int)((uint32_t)(i * ((uint32_t)~(uint32_t)0/255)) >> (sizeof(uint32_t) - 1) * CHAR_BIT); // count
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
    m_map.m_scalars().modified = false;
    m_map.m_scalars().size = size;
  }
};

template <typename UMap>
struct UnorderedMapPrint
{
  typedef UMap map_type;
  typedef typename map_type::device_type device_type;
  typedef typename map_type::size_type size_type;

  map_type m_map;

  UnorderedMapPrint( map_type const& map)
    : m_map(map)
  {}

  void apply() const
  {
    parallel_for(m_map.m_hash_lists.size(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const
  {
    const size_type invalid_index = map_type::invalid_index;
    size_type curr = m_map.m_hash_lists(i);
    while (curr != invalid_index) {
      printf("%d %d\n", m_map.m_keys[curr], i);
      curr = m_map.m_next_index[curr];
    }
  }
};

template <typename UMap>
struct UnorderedMapCountFailedInserts
{
  typedef UMap map_type;
  typedef typename map_type::device_type device_type;
  typedef typename map_type::size_type size_type;
  typedef uint32_t value_type;

  map_type m_map;

  UnorderedMapCountFailedInserts( map_type const& map)
    : m_map(map)
  {}

  void apply() const
  {
    parallel_reduce(m_map.m_failed_insert_scratch.size(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & failed_inserts)
  {
    failed_inserts = 0;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & failed_inserts, const volatile size_type & incr )
  {
    failed_inserts += incr;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & failed_inserts) const
  {
    failed_inserts += m_map.m_failed_insert_scratch[i];
    m_map.m_failed_insert_scratch[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void final( value_type & failed_inserts ) const
  {
    m_map.m_scalars().failed_inserts = failed_inserts;
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
