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

/// \file Kokkos_UnorderedMap.hpp
/// \brief Declaration and definition of Kokkos::UnorderedMap.
///
/// This header file declares and defines Kokkos::UnorderedMap and its
/// related nonmember functions.

#ifndef KOKKOS_UNORDERED_MAP_HPP
#define KOKKOS_UNORDERED_MAP_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Functional.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Atomic.hpp>
#include <Kokkos_HostSpace.hpp>

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_UnorderedMap_impl.hpp>

#include <iostream>

#include <stdint.h>

#include <stdexcept>

#if (defined( __GNUC__ ) || defined( __GNUG__ )) && not defined( __CUDACC__ )

#define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) __builtin_prefetch(addr,0,0)
#define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) __builtin_prefetch(addr,1,0)

#else

#define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) ((void)0)
#define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) ((void)0)

#endif

namespace Kokkos {

enum { UnorderedMapInvalidIndex = ~0u };

/// \brief First element of the return value of UnorderedMap::insert().
///
/// Inserting an element into an UnorderedMap is not guaranteed to
/// succeed.  There are three possible conditions:
/// <ol>
/// <li> <tt>INSERT_FAILED</tt>: The insert failed.  This usually
///      means that the UnorderedMap ran out of space. </li>
/// <li> <tt>INSERT_SUCCESS</tt>: The insert succeeded, and the key
///      did <i>not</i> exist in the table before. </li>
/// <li> <tt>INSERT_EXISTING</tt>: The insert succeeded, and the key
///      <i>did</i> exist in the table before.  The new value was
///      ignored and the old value was left in place. </li>
/// </ol>

class UnorderedMapInsertResult
{
private:
  enum Status{
     SUCCESS = 1u << 31
   , EXISTING = 1u << 30
   , FREED_EXISTING = 1u << 29
   , LOOP_COUNT_MASK = ~(SUCCESS | EXISTING | FREED_EXISTING)
  };

public:
  /// Did the map successful insert the key/value pair
  KOKKOS_FORCEINLINE_FUNCTION
  bool success() const { return (m_status & SUCCESS); }

  /// Was the key already present in the map
  KOKKOS_FORCEINLINE_FUNCTION
  bool existing() const { return (m_status & EXISTING); }

  /// Did the map fail to insert the key due to insufficent capacity
  KOKKOS_FORCEINLINE_FUNCTION
  bool failed() const { return m_index == UnorderedMapInvalidIndex; }

  /// Did the map lose a race condition to insert a dupulicate key/value pair
  /// where an index was claimed that needed to be released
  KOKKOS_FORCEINLINE_FUNCTION
  bool freed_existing() const { return (m_status & FREED_EXISTING); }

  /// How many iterations through the insert loop did it take before the
  /// map returned
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t loop_iterations() const { return (m_status & LOOP_COUNT_MASK); }

  /// Index where the key can be found as long as the insert did not fail
  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t index() const { return m_index; }

  KOKKOS_FORCEINLINE_FUNCTION
  UnorderedMapInsertResult()
    : m_index(UnorderedMapInvalidIndex)
    , m_status(0)
  {}

  KOKKOS_FORCEINLINE_FUNCTION
  void increment_loop_counter()
  {
    m_status += (loop_iterations() < LOOP_COUNT_MASK) ? 1u : 0u;
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void set_existing(uint32_t i, bool arg_freed_existing)
  {
    m_index = i;
    m_status = EXISTING | (arg_freed_existing ? FREED_EXISTING : 0u) | loop_iterations();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  void set_success(uint32_t i)
  {
    m_index = i;
    m_status = SUCCESS | loop_iterations();
  }

private:
  uint32_t m_index;
  uint32_t m_status;
};

/// \class UnorderedMap
/// \brief Thread-safe, performance-portable lookup table.
///
/// This class provides a lookup table.  In terms of functionality,
/// this class compares to std::unordered_map (new in C++11).
/// "Unordered" means that keys are not stored in any particular
/// order, unlike (for example) std::map.  "Thread-safe" means that
/// lookups, insertion, and deletion are safe to call by multiple
/// threads in parallel.  "Performance-portable" means that parallel
/// performance of these operations is reasonable, on multiple
/// hardware platforms.  Platforms on which performance has been
/// tested include conventional Intel x86 multicore processors, Intel
/// Xeon Phi ("MIC"), and NVIDIA GPUs.
///
/// Parallel performance portability entails design decisions that
/// might differ from one's expectation for a sequential interface.
/// This particularly affects insertion of single elements.  In an
/// interface intended for sequential use, insertion might reallocate
/// memory if the original allocation did not suffice to hold the new
/// element.  In this class, insertion does <i>not</i> reallocate
/// memory.  This means that it might fail.  insert() returns an enum
/// which indicates whether the insert failed.  There are three
/// possible conditions:
/// <ol>
/// <li> <tt>INSERT_FAILED</tt>: The insert failed.  This usually
///      means that the UnorderedMap ran out of space. </li>
/// <li> <tt>INSERT_SUCCESS</tt>: The insert succeeded, and the key
///      did <i>not</i> exist in the table before. </li>
/// <li> <tt>INSERT_EXISTING</tt>: The insert succeeded, and the key
///      <i>did</i> exist in the table before.  The new value was
///      ignored and the old value was left in place. </li>
/// </ol>
///
/// \tparam Key Type of keys of the lookup table.  If \c const, users
///   are not allowed to add or remove keys, though they are allowed
///   to change values.  In that case, the implementation may make
///   optimizations specific to the <tt>Device</tt>.  For example, if
///   <tt>Device</tt> is \c Cuda, it may use texture fetches to access
///   keys.
///
/// \tparam Value Type of values stored in the lookup table.  You may use
///   \c void here, in which case the table will be a set of keys.  If
///   \c const, users are not allowed to change entries.
///   In that case, the implementation may make
///   optimizations specific to the \c Device, such as using texture
///   fetches to access values.
///
/// \tparam Device The Kokkos Device type.
///
/// \tparam Hasher Definition of the hash function for instances of
///   <tt>Key</tt>.  The default will calculate a bitwise hash.
///
/// \tparam EqualTo Definition of the equality function for instances of
///   <tt>Key</tt>.  The default will do a bitwise equality comparison.
///
template <   typename Key
           , typename Value
           , typename Device
           , typename Hasher = pod_hash<typename Impl::remove_const<Key>::type>
           , typename EqualTo = pod_equal_to<typename Impl::remove_const<Key>::type>
        >
class UnorderedMap
{
public:
  //! \name Public types and constants
  //@{

  //key_types
  typedef Key declared_key_type;
  typedef typename Impl::remove_const<declared_key_type>::type key_type;
  typedef typename Impl::add_const<key_type>::type const_key_type;

  //value_types
  typedef Value declared_value_type;
  typedef typename Impl::remove_const<declared_value_type>::type value_type;
  typedef typename Impl::add_const<value_type>::type const_value_type;

  typedef Device device_type;
  typedef Hasher hasher_type;
  typedef EqualTo  equal_to_type;
  typedef uint32_t size_type;

  //map_types
  typedef UnorderedMap<declared_key_type,declared_value_type,device_type,hasher_type,equal_to_type> declared_map_type;
  typedef UnorderedMap<key_type,value_type,device_type,hasher_type,equal_to_type>                   insertable_map_type;
  typedef UnorderedMap<const_key_type,value_type,device_type,hasher_type,equal_to_type>             modifiable_map_type;
  typedef UnorderedMap<const_key_type,const_value_type,device_type,hasher_type,equal_to_type>       const_map_type;

  static const bool is_set = Impl::is_same<void,value_type>::value;
  static const bool has_const_key = Impl::is_same<const_key_type,declared_key_type>::value;
  static const bool has_const_value = is_set || Impl::is_same<const_value_type,declared_value_type>::value;

  static const bool is_insertable_map = !has_const_key && (is_set || !has_const_value);
  static const bool is_modifiable_map = has_const_key && !has_const_value;
  static const bool is_const_map = has_const_key && has_const_value;


  typedef UnorderedMapInsertResult insert_result;

  typedef UnorderedMap<Key,Value,typename Device::host_mirror_device_type,Hasher,EqualTo> HostMirror;

  typedef Impl::UnorderedMapHistogram<const_map_type> histogram_type;

  //@}

private:
  typedef uint32_t block_type;
  enum { invalid_index = UnorderedMapInvalidIndex   };
  enum { block_size = sizeof(block_type)*CHAR_BIT   };
  enum { block_shift = Impl::power_of_two<block_size>::value };
  enum { intial_block = ~static_cast<block_type>(0) };


  typedef typename Impl::if_c< is_set, int, declared_value_type>::type impl_value_type;

  typedef typename Impl::if_c<   is_insertable_map
                               , View< key_type *, device_type>
                               , View< const key_type *, device_type, MemoryTraits<RandomAccess> >
                             >::type key_type_view;

  typedef typename Impl::if_c<   is_insertable_map || is_modifiable_map
                               , View< impl_value_type *, device_type>
                               , View< const impl_value_type *, device_type, MemoryTraits<RandomAccess> >
                             >::type value_type_view;

  typedef typename Impl::if_c<   is_insertable_map
                               , View< size_type *, device_type>
                               , View< const size_type *, device_type, MemoryTraits<RandomAccess> >
                             >::type size_type_view;

  typedef typename Impl::if_c<   is_insertable_map
                               , View< block_type *, device_type>
                               , View< const block_type *, device_type, MemoryTraits<RandomAccess> >
                             >::type block_type_view;

  typedef View< Impl::UnorderedMapScalars, device_type> scalars_view;

  typedef Kokkos::Impl::DeepCopy< Kokkos::HostSpace, typename device_type::memory_space > raw_deep_copy;

public:
  //! \name Public member functions
  //@{

  UnorderedMap()
    : m_hasher()
    , m_equal_to()
    , m_capacity()
    , m_available_indexes()
    , m_hash_lists()
    , m_next_index()
    , m_keys()
    , m_values()
    , m_scalars()
  {}

  /// \brief Constructor
  ///
  /// \param requested_capacity [in] Initial requested maximum number of
  ///   entries in the hash table.
  /// \param hash [in] Hasher function for \c Key instances.  The
  ///   default value usually suffices.
  UnorderedMap(  size_type requested_capacity , hasher_type hasher = hasher_type(), equal_to_type equal_to = equal_to_type() )
    : m_hasher(hasher)
    , m_equal_to(equal_to)
    , m_capacity( (((requested_capacity + block_size - 1) >> block_shift) << block_shift) )
    , m_available_indexes(AllocateWithoutInitializing(), "UnorderedMap available indexes", (m_capacity >> block_shift))
    , m_hash_lists(AllocateWithoutInitializing(), "UnorderedMap hash list", Impl::find_hash_size(m_capacity))
    , m_next_index(AllocateWithoutInitializing(), "UnorderedMap next index", m_capacity+1) // +1 so that the *_at functions can always return a valid reference
    , m_keys("UnorderedMap keys",m_capacity+1)
    , m_values("UnorderedMap values",(is_set? 1 : m_capacity+1))
    , m_scalars("UnorderedMap scalars")
  {
    if (!is_insertable_map) {
      throw std::runtime_error("Cannot construct a non-insertable (i.e. const key_type) unordered_map");
    }

    Kokkos::deep_copy(m_available_indexes, intial_block);
    Kokkos::deep_copy(m_hash_lists,invalid_index);
    Kokkos::deep_copy(m_next_index,invalid_index);
  }

  histogram_type get_histogram()
  {
    return histogram_type(*this);
  }

  //! Clear all entries in the table.
  void clear()
  {
    if (m_capacity == 0) return;

    Kokkos::deep_copy(m_available_indexes, intial_block);
    Kokkos::deep_copy(m_hash_lists,invalid_index);
    Kokkos::deep_copy(m_next_index,invalid_index);
    {
      const key_type tmp = key_type();
      Kokkos::deep_copy(m_keys,tmp);
    }
    if (is_set){
      const impl_value_type tmp = impl_value_type();
      Kokkos::deep_copy(m_values,tmp);
    }
    {
      const Impl::UnorderedMapScalars tmp = Impl::UnorderedMapScalars();
      Kokkos::deep_copy(m_scalars,tmp);
    }

  }

  /// \brief Change the capacity of the the map
  ///
  /// If there are no failed inserts the current size of the map will
  /// be used as a lower bound for the input capacity.
  /// If the map is not empty and does not have failed inserts
  /// and the capacity changes then the current data is copied
  /// into the resized / rehashed map.
  ///
  /// This is <i>not</i> a device function; it may <i>not</i> be
  /// called in a parallel kernel.
  bool rehash(size_type new_capacity = 0)
  {
    if(!is_insertable_map) return false;

    if ( new_capacity != m_capacity ) {

      const size_type curr_size = size();
      const bool copy_data = (curr_size > 0u) && !has_failed_inserts();
      new_capacity = (copy_data && (new_capacity < curr_size)) ? curr_size : new_capacity;

      declared_map_type tmp(new_capacity, m_hasher);

      if (copy_data ) {
        Impl::UnorderedMapRehash<declared_map_type> f(tmp,*this);
        f.apply();
      }
      *this = tmp;
    }
    else if ( has_failed_inserts() ) {
      clear();
    }

    return true;
  }

  bool rehash(size_type new_capacity, hasher_type hasher)
  {
    if(!is_insertable_map) return false;

    const size_type curr_size = size();
    const bool copy_data = (curr_size > 0u) && !has_failed_inserts();
    new_capacity = (copy_data && (new_capacity < curr_size)) ? curr_size : new_capacity;

    declared_map_type tmp(new_capacity, hasher);

    if (copy_data ) {
      Impl::UnorderedMapRehash<declared_map_type> f(tmp,*this);
      f.apply();
    }

    *this = tmp;

    return true;
  }

  /// \brief The number of entries in the table.
  ///
  /// This method has undefined behavior when erasable() is true.
  ///
  /// Note that this is not a device function; it cannot be called in
  /// a parallel kernel.  The value is not stored as a variable; it
  /// must be computed.
  size_type size() const
  {
    if( m_capacity == 0u ) return 0u;
    sync_scalars();
    size_type result;
    raw_deep_copy(&result,&m_scalars.ptr_on_device()->size, sizeof(size_type));
    return result;
  }

  /// \brief The current number of failed insert() calls.
  ///
  /// This is <i>not</i> a device function; it may <i>not</i> be
  /// called in a parallel kernel.  The value is not stored as a
  /// variable; it must be computed.
  bool has_failed_inserts() const
  {
    if( m_capacity == 0u ) return false;
    bool result;
    raw_deep_copy(&result,&m_scalars.ptr_on_device()->has_failed_inserts, sizeof(size_type));
    return result;
  }

  bool erasable() const
  {
    if( m_capacity == 0u ) return false;
    bool result = false;
    if (is_insertable_map){
      raw_deep_copy(&result,&m_scalars.ptr_on_device()->erasable, sizeof(bool));
    }
    return result;
  }

  bool begin_erase()
  {
    bool result = !erasable();
    if (is_insertable_map && result) {
      device_type::fence();
      const bool true_ = true;
      typedef Kokkos::Impl::DeepCopy< typename device_type::memory_space, Kokkos::HostSpace > copy_to_device;
      copy_to_device(&m_scalars.ptr_on_device()->erasable, &true_, sizeof(bool) );
      device_type::fence();
    }
    return result;
  }

  bool end_erase()
  {
    bool result = erasable();
    if (is_insertable_map && result) {
      device_type::fence();
      Impl::UnorderedMapErase<declared_map_type> f(*this);
      f.apply();
      const bool false_ = false;
      typedef Kokkos::Impl::DeepCopy< typename device_type::memory_space, Kokkos::HostSpace > copy_to_device;
      copy_to_device(&m_scalars.ptr_on_device()->erasable, &false_, sizeof(bool) );
      sync_scalars(true);
    }
    return result;
  }

  /// \brief The maximum number of entries that the table can hold.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_INLINE_FUNCTION
  size_type capacity() const
  { return m_capacity; }

  /// \brief The number of hash table "buckets."
  ///
  /// This is different than the number of entries that the table can
  /// hold.  Each key hashes to an index in [0, hash_capacity() - 1].
  /// That index can hold zero or more entries.  This class decides
  /// what hash_capacity() should be, given the user's upper bound on
  /// the number of entries the table must be able to hold.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_INLINE_FUNCTION
  size_type hash_capacity() const
  { return m_hash_lists.size(); }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.  As discussed in the class documentation, it need not
  /// succeed.  The return value tells you if it did.
  ///
  /// \param k [in] The key to attempt to insert.
  /// \param v [in] The corresponding value to attempt to insert.  If
  ///   using this class as a set (with Value = void), then you need not
  ///   provide this value.
  KOKKOS_INLINE_FUNCTION
  insert_result insert(key_type const& k, impl_value_type const&v = impl_value_type()) const
  {
    insert_result result;

    if ( is_insertable_map && 0u < m_capacity && ! m_scalars().erasable ) {

      bool volatile & has_failed_inserts_ref = m_scalars().has_failed_inserts ;

      const size_type hash_value = m_hasher(k);
      const size_type hash_list = hash_value % m_hash_lists.size();
      const size_type num_blocks = m_available_indexes.size();

      // Force integer multiply to long
      const size_type start_block_index =
        static_cast<size_type>( ( static_cast<uint64_t>(hash_list) * num_blocks ) / m_hash_lists.size() );

      size_type * curr_ptr   = & m_hash_lists[ hash_list ];
      size_type new_index    = invalid_index ;

      size_type curr_block_index     = start_block_index ;
      block_type * curr_block_ptr    = & m_available_indexes[ curr_block_index ];
      block_type new_index_mask  = 0 ;

      bool not_done = true ;

#if defined( __MIC__ )
      #pragma noprefetch
#endif
      while ( not_done && !has_failed_inserts_ref ) {
        result.increment_loop_counter();

        // Continue searching the unordered list for this key,
        // list will only be appended during insert phase.
        // Need volatile_load as other threads may be appending.
        size_type curr = volatile_load(curr_ptr);

        KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
#if defined( __MIC__ )
        #pragma noprefetch
#endif
        while ( curr != invalid_index && ! m_equal_to( volatile_load(&m_keys[curr]), k) ) {
          curr_ptr = &m_next_index[curr];
          curr = volatile_load(curr_ptr);
          KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
        }

        //------------------------------------------------------------
        // If key already present then return that index.
        if ( curr != invalid_index ) {

          const bool free_existing = new_index != invalid_index;
          if ( free_existing ) {
            // Previously claimed an unused entry that was not inserted.
            // Release this unused entry immediately.
            atomic_fetch_or( curr_block_ptr , new_index_mask );
          }

          result.set_existing(curr, free_existing);
          not_done = false ;
        }
        //------------------------------------------------------------
        // Key is not currently in the map.
        // If the thread has claimed an entry try to insert now.
        else if ( new_index != invalid_index ) {

          // Attempt to append claimed entry into the list.
          // Another thread may also be trying to append the same list so protect with atomic.
          curr = atomic_compare_exchange(curr_ptr, static_cast<size_type>(invalid_index), new_index);

          // Succeeded in appending
          if ( curr == invalid_index ) {
            if ( ! m_scalars().modified ) { m_scalars().modified = true ; }
            result.set_success(new_index);
            not_done = false ;
          }
        }
        //------------------------------------------------------------
        // If have not already claimed an unused entry then do so now.
        else {

          const block_type curr_block = volatile_load(curr_block_ptr) ;

          // Block has an unused entry, try to claim that entry, a race condition managed via atomics
          if ( curr_block ) {

            // Find an unused offset within the current block
            const int offset = get_offset(hash_list, curr_block);

            // Try to unset that bit:
            new_index_mask = static_cast<block_type>(1) << offset ;

            // Try to claim the offset
            if ( atomic_fetch_and(curr_block_ptr, ~new_index_mask) & new_index_mask) {
              new_index = ( curr_block_index << block_shift ) + static_cast<size_type>(offset);

              // Set key and value
              KOKKOS_NONTEMPORAL_PREFETCH_STORE(&m_keys[new_index]);
              m_keys[new_index] = k ;

              if (!is_set) {
                KOKKOS_NONTEMPORAL_PREFETCH_STORE(&m_values[new_index]);
                m_values[new_index] = v ;
              }

              // Do not proceed until key and value are updated in global memory
              memory_fence();
            }
          }
          // curr_block is empty, advanced to the next nonempty block
          else {
            do {
              curr_block_index = (curr_block_index + 1u) % num_blocks;
              curr_block_ptr = & m_available_indexes[ curr_block_index ];
              // Wrapped completely around to start, map is (nearly) full, set failure flag
              if ( start_block_index == curr_block_index ) {
                has_failed_inserts_ref = true;
                not_done = false ;
              }
            } while ( !volatile_load(curr_block_ptr) && !has_failed_inserts_ref );
          }
        }
      } // while ( not_done )
    }

    return result ;
  }

  KOKKOS_INLINE_FUNCTION
  bool erase(key_type const& k) const
  {
    bool result = false;

    if(is_insertable_map && 0u < m_capacity && m_scalars().erasable) {
      size_type index = find(k);
      if (valid_at(index)) {
        free_index(index);
        result = true;
      }
    }

    return result;
  }

  /// \brief Find the given key \c k, if it exists in the table.
  ///
  /// \return If the key exists in the table, the index of the
  ///   value corresponding to that key; otherwise, an invalid index.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_INLINE_FUNCTION
  size_type find( const key_type & k) const
  {
    size_type curr = 0u < m_capacity ? m_hash_lists( m_hasher(k) % m_hash_lists.size() ) : invalid_index ;

    KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
    while (curr != invalid_index && !m_equal_to( m_keys[curr], k) ) {
      KOKKOS_NONTEMPORAL_PREFETCH_LOAD(&m_keys[curr != invalid_index ? curr : 0]);
      curr = m_next_index[curr];
    }

    return curr;
  }

  /// \brief Does the key exist in the map
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_INLINE_FUNCTION
  bool exists( const key_type & k) const
  {
    return find(k) != invalid_index;
  }


  /// \brief Get the value with \c i as its direct index.
  ///
  /// \param i [in] Index directly into the array of entries.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  ///
  /// 'const value_type' via Cuda texture fetch must return by value.
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::if_c< (is_set || has_const_value), impl_value_type, impl_value_type &>::type
  value_at(size_type i) const
  {
    return m_values[ is_set ? 0 : (i < m_capacity ? i : m_capacity) ];
  }

  /// \brief Get the key with \c i as its direct index.
  ///
  /// \param i [in] Index directly into the array of entries.
  ///
  /// This <i>is</i> a device function; it may be called in a parallel
  /// kernel.
  KOKKOS_FORCEINLINE_FUNCTION
  key_type key_at(size_type i) const
  {
    return m_keys[ i < m_capacity ? i : m_capacity ];
  }

  KOKKOS_FORCEINLINE_FUNCTION
  bool valid_at(size_type i) const
  {
    if (i >= m_capacity) return false;
    const block_type block = m_available_indexes[i >> block_shift];
    const int offset = static_cast<int>(i & (block_size-1u));

    return !(block & ( static_cast<block_type>(1u) << offset));
  }

  template <typename SKey, typename SValue>
  UnorderedMap( UnorderedMap<SKey,SValue,Device,Hasher,EqualTo> const& src,
                typename Impl::enable_if< Impl::UnorderedMapCanAssign<declared_key_type,declared_value_type,SKey,SValue>::value,int>::type = 0
              )
    : m_hasher(src.m_hasher)
    , m_equal_to(src.m_equal_to)
    , m_capacity(src.m_capacity)
    , m_available_indexes(src.m_available_indexes)
    , m_hash_lists(src.m_hash_lists)
    , m_next_index(src.m_next_index)
    , m_keys(src.m_keys)
    , m_values(src.m_values)
    , m_scalars(src.m_scalars)
  {}


  template <typename SKey, typename SValue>
  typename Impl::enable_if< Impl::UnorderedMapCanAssign<declared_key_type,declared_value_type,SKey,SValue>::value
                           ,declared_map_type & >::type
  operator=( UnorderedMap<SKey,SValue,Device,Hasher,EqualTo> const& src)
  {
    m_hasher = src.m_hasher;
    m_equal_to = src.m_equal_to;
    m_capacity = src.m_capacity;
    m_available_indexes = src.m_available_indexes;
    m_hash_lists = src.m_hash_lists;
    m_next_index = src.m_next_index;
    m_keys = src.m_keys;
    m_values = src.m_values;
    m_scalars = src.m_scalars;
    return *this;
  }

  template <typename SKey, typename SValue, typename SDevice>
  typename Impl::enable_if< Impl::is_same< typename Impl::remove_const<SKey>::type, key_type>::value &&
                            Impl::is_same< typename Impl::remove_const<SValue>::type, value_type>::value
                          >::type
  create_copy_view( UnorderedMap<SKey, SValue, SDevice, Hasher,EqualTo> const& src)
  {
    if (m_available_indexes.ptr_on_device() != src.m_available_indexes.ptr_on_device()) {
      typedef Kokkos::Impl::DeepCopy< typename device_type::memory_space, typename SDevice::memory_space > deep_copy_pointer;

      src.sync_scalars();
      *this = insertable_map_type(src.capacity(), src.m_hasher);

      deep_copy_pointer(m_available_indexes.ptr_on_device(), src.m_available_indexes.ptr_on_device(), sizeof(block_type) * src.m_available_indexes.size());
      deep_copy_pointer(m_hash_lists.ptr_on_device(), src.m_hash_lists.ptr_on_device(), sizeof(size_type) * src.m_hash_lists.size());
      deep_copy_pointer(m_next_index.ptr_on_device(), src.m_next_index.ptr_on_device(), sizeof(size_type) * src.m_next_index.size());
      deep_copy_pointer(m_keys.ptr_on_device(), src.m_keys.ptr_on_device(), sizeof(key_type) * src.m_keys.size());
      if (!is_set) deep_copy_pointer(m_values.ptr_on_device(), src.m_values.ptr_on_device(), sizeof(value_type) * src.m_values.size());
      deep_copy_pointer(m_scalars.ptr_on_device(), src.m_scalars.ptr_on_device(), sizeof(Impl::UnorderedMapScalars));
    }
  }

  //@}
private: // private member functions

  KOKKOS_FORCEINLINE_FUNCTION
  int get_offset(size_type hash_list, block_type block) const
  {
    const int rotate = (hash_list & 31u);
    block = rotate ? ((block >> rotate) | (block << static_cast<int>(block_size-rotate))) : block;
    return (Impl::bit_scan_forward(block) + rotate) & (block_size-1);
  }

  KOKKOS_INLINE_FUNCTION
  bool free_index(size_type i) const
  {
    if (!is_insertable_map) return false;

    block_type * available_ptr = &m_available_indexes[ i >> block_shift ];
    const block_type bit_mask = static_cast<block_type>(1u) << static_cast<int>(i & (block_size-1));

    atomic_or(available_ptr, bit_mask);

    return true;
  }

  void sync_scalars(bool force_sync = false) const
  {
    if( m_capacity == 0u ) return;
    bool modified = false;
    raw_deep_copy(&modified, &m_scalars.ptr_on_device()->modified, sizeof(bool) );
    if (force_sync || modified)
    {
      device_type::fence();
      {
        Impl::UnorderedMapSize<const_map_type> f(*this);
        f.apply();
      }

      // make sure the results are stored before continuing
      device_type::fence();
    }
  }

private: // private members
  hasher_type      m_hasher;
  equal_to_type    m_equal_to;
  size_type        m_capacity;
  block_type_view  m_available_indexes;
  size_type_view   m_hash_lists;
  size_type_view   m_next_index;
  key_type_view    m_keys;
  value_type_view  m_values;
  scalars_view     m_scalars;

  template <typename KKey, typename VValue, typename DDevice, typename HHash, typename EEqualTo>
  friend class UnorderedMap;

  template <typename UMap>
  friend struct Impl::UnorderedMapSize;

  template <typename UMap>
  friend struct Impl::UnorderedMapErase;

  template <typename UMap>
  friend struct Impl::UnorderedMapHistogram;
};

// Specialization of deep_copy for two UnorderedMap objects.
template <  typename DKey, typename DT, typename DDevice
          , typename SKey, typename ST, typename SDevice
          , typename Hasher, typename EqualTo >
inline void deep_copy(         UnorderedMap<DKey, DT, DDevice, Hasher, EqualTo> & dst
                       , const UnorderedMap<SKey, ST, SDevice, Hasher, EqualTo> & src )
{
  dst.create_copy_view(src);
}


} // namespace Kokkos

#endif //KOKKOS_UNORDERED_MAP_HPP
