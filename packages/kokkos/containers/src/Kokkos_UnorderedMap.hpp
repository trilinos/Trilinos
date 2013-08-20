#ifndef KOKKOS_UNORDERED_MAP_HPP
#define KOKKOS_UNORDERED_MAP_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Functional.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Atomic.hpp>

#include <impl/Kokkos_ArrayTraits.hpp>
#include <impl/Kokkos_UnorderedMap_impl.hpp>

#include <iostream>

#include <stdint.h>

namespace Kokkos {
enum unordered_map_insert_state
{
    INSERT_FAILED
  , INSERT_SUCCESS
  , INSERT_EXISTING
};

template <   typename Key
           , typename T
           , typename Device
           , typename Compare = less<typename Impl::remove_const<Key>::type>
           , typename Hash = hash<typename Impl::remove_const<Key>::type>
        >
class unordered_map
{
public: // public types and constants
  typedef Impl::UnorderedMap::map_data<Key,T,Device,Compare,Hash> map_data;
  typedef Impl::UnorderedMap::node_atomic node_atomic;

  typedef typename map_data::device_type device_type;
  typedef typename map_data::compare_type compare_type;
  typedef typename map_data::hash_type hash_type;
  typedef typename map_data::key_type key_type;
  typedef typename map_data::mapped_type mapped_type;
  typedef typename map_data::value_type value_type;
  typedef typename map_data::pointer pointer;
  typedef typename map_data::const_pointer const_pointer;
  typedef typename map_data::node_type node_type;
  typedef typename map_data::node_block_type node_block_type;
  typedef typename map_data::size_type size_type;

  typedef pair<unordered_map_insert_state, pointer> insert_result;

private:

  typedef typename Impl::if_c<  map_data::has_void_mapped_type
                              , int
                              , mapped_type
                             >::type insert_mapped_type;


public: //public member functions

  unordered_map(  uint32_t arg_num_nodes
                , compare_type compare = compare_type()
                , hash_type hash = hash_type()
               )
    : m_data(  arg_num_nodes
             , compare
             , hash
            )
  {}

  void check_sanity() const
  { m_data.check_sanity(); }

  uint32_t size() const
  {  return m_data.size(); }

  uint32_t unused() const
  {  return m_data.unused(); }

  uint32_t pending_delete() const
  {  return m_data.pending_delete(); }

  uint32_t failed_inserts() const
  { return m_data.failed_inserts(); }

  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return m_data.capacity(); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_capacity() const
  { return m_data.hash_capacity(); }

  void remove_pending_delete() const
  {  return m_data.remove_pending_delete_keys(); }

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


  KOKKOS_INLINE_FUNCTION
  insert_result insert(const key_type & k, const insert_mapped_type & v = insert_mapped_type()) const
  {
    m_data.set_modified();

    insert_result result(INSERT_FAILED,NULL);

    const uint32_t hash_value = m_data.key_hash(k);
    const uint32_t hash_index = hash_value % m_data.hashes.size();

    uint32_t node_index = node_atomic::invalid_next;

    bool curr_equal = false;
    uint32_t curr_index = node_atomic::invalid_next;
    volatile uint64_t * prev_atomic = & m_data.hashes[hash_index].value;
    uint64_t prev = 0u;

    find_previous(k,prev_atomic,prev,curr_equal,curr_index);

    do {
      if (curr_equal) {
        if (node_index != node_atomic::invalid_next) {
          // release any node that was claimed by this thread
          m_data.get_node(node_index).atomic = node_atomic::make_atomic(node_atomic::invalid_next, Impl::UnorderedMap::UNUSED);
#if defined( __CUDA_ARCH__ )
          __threadfence();
#endif
          volatile int * used_count = &m_data.node_blocks[node_index>>node_block_type::shift].used_count;
          atomic_fetch_add(used_count, -1);
        }
        // Node already exist
        result = insert_result(INSERT_EXISTING, &m_data.get_node(curr_index).value);
        break;
      }
      else {
        // try to insert here
        if (node_index == node_atomic::invalid_next) {
          node_index = find_unused_node(hash_value);
          if (node_index == node_atomic::invalid_next) {
            // unable to obtain an unused node
            break;
          }
        }
        // this thread has unique control of the node
        // so can construct the value and set up the state and next index
        node_type & n = m_data.get_node(node_index);
        n.destruct_value();
        n.construct_value(value_type(k,v));
        n.atomic = node_atomic::make_atomic( curr_index, Impl::UnorderedMap::USED);

        uint64_t new_atomic = node_atomic::make_atomic( node_index, node_atomic::state(prev));

#if defined( __CUDA_ARCH__ )
        __threadfence();
#endif
        const bool ok = atomic_compare_exchange_strong( prev_atomic, prev, new_atomic);
        if ( ok ) {
          // successfully inserted the node
          result = insert_result(INSERT_SUCCESS, &n.value);
          break;
        }
      }
      // insert failed -- find correct insertion point again
      find_previous(k,prev_atomic,prev,curr_equal,curr_index);
    } while (true);
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  void mark_pending_delete(const key_type & k) const
  {
    m_data.set_modified();

    const uint32_t hash_value = m_data.key_hash(k);
    const uint32_t hash_index = hash_value % m_data.hashes.size();

    uint32_t node_index = node_atomic::invalid_next;

    bool curr_equal = false;
    uint32_t curr_index = node_atomic::invalid_next;
    volatile uint64_t * prev_atomic = & m_data.hashes[hash_index].value;
    uint64_t prev = 0u;

    find_previous(k,prev_atomic,prev,curr_equal,curr_index);

    do {
      if (curr_equal) {
        if (node_index != node_atomic::invalid_next) {
          // release any node that was claimed by this thread
          m_data.get_node(node_index).atomic = node_atomic::make_atomic(node_atomic::invalid_next, Impl::UnorderedMap::UNUSED);
#if defined( __CUDA_ARCH__ )
          __threadfence();
#endif
          volatile int * used_count = &m_data.node_blocks[node_index>>node_block_type::shift].used_count;
          atomic_fetch_add(used_count, -1);
        }
        // mark the current node as deleted
        volatile uint64_t * curr_atomic_ptr = &m_data.get_node(curr_index).atomic.value;
        uint64_t curr_atomic = *curr_atomic_ptr;
        while ( node_atomic::state(curr_atomic) == Impl::UnorderedMap::USED) {
          uint64_t new_atomic = node_atomic::make_atomic( node_atomic::next(curr_atomic), Impl::UnorderedMap::PENDING_DELETE);
          curr_atomic = atomic_compare_exchange(curr_atomic_ptr,curr_atomic,new_atomic);
        }
        return;
      }
      else {
        // key does not exist
        // insert a node with the given key marked as deleted
        if (node_index == node_atomic::invalid_next) {
          node_index = find_unused_node(hash_value);
          if (node_index == node_atomic::invalid_next) {
            return;
          }
        }

        // this thread has unique control of the node
        // so can construct the value and set up the state and next index
        node_type & n = m_data.get_node(node_index);
        n.destruct_value();
        n.construct_value(value_type(k,insert_mapped_type()));
        n.atomic = node_atomic::make_atomic( curr_index, Impl::UnorderedMap::PENDING_DELETE);

        uint64_t new_atomic = node_atomic::make_atomic( node_index, node_atomic::state(prev));

#if defined( __CUDA_ARCH__ )
        __threadfence();
#endif

        const bool ok = atomic_compare_exchange_strong( prev_atomic, prev, new_atomic);
        if ( ok ) {
          return;
        }
      }
      // insert failed -- find correct insertion point again
      find_previous(k,prev_atomic,prev,curr_equal,curr_index);
    } while (true);
  }

  // TODO protect with enable_if
  KOKKOS_INLINE_FUNCTION
  void mark_pending_delete( const_pointer p ) const
  {
    if (p) mark_pending_delete(p->first);
  }

  KOKKOS_INLINE_FUNCTION
  pointer find( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    return (node_index != node_atomic::invalid_next) ? &m_data.get_node(node_index).value : NULL;
  }


  KOKKOS_INLINE_FUNCTION
  pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i < m_data.capacity();
    const bool used_node  = node_atomic::state(m_data.get_node(i).atomic) == Impl::UnorderedMap::USED;

    return valid_range && used_node ? &m_data.get_node(i).value : NULL;
  }

private: // private member functions

  KOKKOS_INLINE_FUNCTION
  uint32_t find_unused_node(uint32_t hash_value) const
  {
    const uint32_t num_blocks = m_data.node_blocks.size();
    const uint32_t start_block = hash_value % num_blocks;
    const uint32_t end_block = start_block + num_blocks;

    if (m_data.no_failed_inserts()) {

      for (uint32_t i = start_block; i < end_block; ++i) {
        if (!m_data.no_failed_inserts()) break;

        const uint32_t block = i % num_blocks;
        volatile int * used_count = &m_data.node_blocks[block].used_count;
        int count = * used_count;
        if (static_cast<unsigned>(count) < node_block_type::size) {
          //stores the old value into count
          const int old_count = atomic_fetch_add(used_count, 1);
          if (static_cast<unsigned>(old_count) < node_block_type::size) {
            //claimed a node in this block keep looping block utill successful at claming a node
            for (uint32_t start_node = (hash_value & node_block_type::mask); true; ++start_node) {
              if (!m_data.no_failed_inserts()) break;
              const uint32_t n = (block*node_block_type::size) + (start_node & node_block_type::mask);
              volatile uint64_t * atomic = &m_data.get_node(n).atomic.value;
              uint64_t value = *atomic;
              if (    (node_atomic::state(value) == Impl::UnorderedMap::UNUSED)
                  && atomic_compare_exchange_strong(atomic, value, node_atomic::make_atomic(node_atomic::invalid_next,Impl::UnorderedMap::PENDING_INSERT)) )
              {
                return n;
              }
            }
          }
          else {
            //unable to claim a node from this block
            atomic_fetch_add(used_count, -1);
          }
        }
      }
      // unable to get a free node -- insert failed
      m_data.set_failed_insert();
    }
    // count the failed insert
    volatile int * failed_inserts = &m_data.node_blocks[start_block].failed_inserts;
    atomic_fetch_add(failed_inserts, 1);
    return node_atomic::invalid_next;
  }

  KOKKOS_INLINE_FUNCTION
  void find_previous(const key_type & k, volatile uint64_t *& prev_atomic, uint64_t & prev, bool &curr_equal, uint32_t & curr_index) const
  {
    curr_equal = false;
    do {
      prev = *prev_atomic;
      curr_index = node_atomic::next(prev);
      const bool curr_invalid = curr_index == node_atomic::invalid_next;

      if (curr_invalid) break;

       // global read of the key
      volatile const key_type * const key_ptr = &m_data.get_node(curr_index).value.first;
      const key_type curr_key = *key_ptr;

      const bool curr_less = m_data.key_compare( curr_key, k);
      const bool curr_greater = m_data.key_compare( k, curr_key);
      curr_equal = !curr_less && !curr_greater;

      if (!curr_less) break;

      prev_atomic = & m_data.get_node(node_atomic::next(prev)).atomic.value;
    } while (true);
  }

//private: // private members
public:
  map_data m_data;

};

template <   typename Key
           , typename T
           , typename Device
           , typename Compare
           , typename Hash
        >
class unordered_map< const Key, T, Device, Compare, Hash>
{
public: // public types and constants
  typedef Impl::UnorderedMap::map_data<const Key, T,Device,Compare,Hash> map_data;
  typedef Impl::UnorderedMap::node_atomic node_atomic;

  typedef typename map_data::device_type device_type;
  typedef typename map_data::compare_type compare_type;
  typedef typename map_data::hash_type hash_type;
  typedef typename map_data::key_type key_type;
  typedef typename map_data::mapped_type mapped_type;
  typedef typename map_data::value_type value_type;
  typedef typename map_data::pointer pointer;
  typedef typename map_data::const_pointer const_pointer;
  typedef typename map_data::node_type node_type;
  typedef typename map_data::size_type size_type;

public: //public member functions

  template <typename UMap>
  unordered_map(  const UMap & umap )
    : m_data( umap.m_data )
  {}

  void check_sanity() const
  { m_data.check_sanity(); }

  uint32_t size() const
  {  return m_data.size(); }

  uint32_t unused() const
  {  return m_data.unused(); }

  uint32_t pending_delete() const
  {  return m_data.pending_delete(); }

  uint32_t failed_inserts() const
  { return m_data.failed_inserts(); }

  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return m_data.capacity(); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_capacity() const
  { return m_data.hash_capacity(); }

  void remove_pending_delete() const
  {  return m_data.remove_pending_delete_keys(); }


  KOKKOS_INLINE_FUNCTION
  pointer find( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    return (node_index != node_atomic::invalid_next) ? &m_data.get_node(node_index).value : NULL;
  }


  KOKKOS_INLINE_FUNCTION
  pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i < m_data.capacity();
    const bool used_node  = node_atomic::state(m_data.get_node(i).atomic) == Impl::UnorderedMap::USED;

    return valid_range && used_node ? &m_data.get_node(i).value : NULL;
  }


  KOKKOS_INLINE_FUNCTION
  void mark_pending_delete(const key_type & k) const
  {
    m_data.set_modified();

    const uint32_t hash_value = m_data.key_hash(k);
    const uint32_t hash_index = hash_value % m_data.hashes.size();

    bool curr_equal = false;
    uint32_t curr_index = node_atomic::invalid_next;
    const volatile uint64_t * prev_atomic = & m_data.hashes[hash_index].value;
    uint64_t prev = 0u;

    find_previous(k,prev_atomic,prev,curr_equal,curr_index);

    do {
      if (curr_equal) {
        // mark the current node as deleted
        volatile uint64_t * curr_atomic_ptr = &m_data.get_node(curr_index).atomic.value;
        uint64_t curr_atomic = *curr_atomic_ptr;
        while ( node_atomic::state(curr_atomic) == Impl::UnorderedMap::USED) {
          uint64_t new_atomic = node_atomic::make_atomic( node_atomic::next(curr_atomic), Impl::UnorderedMap::PENDING_DELETE);
          curr_atomic = atomic_compare_exchange(curr_atomic_ptr,curr_atomic,new_atomic);
        }
        return;
      }
    } while (true);
  }


private:
  KOKKOS_INLINE_FUNCTION
  void find_previous(const key_type & k, const volatile uint64_t *& prev_atomic, uint64_t & prev, bool &curr_equal, uint32_t & curr_index) const
  {
    curr_equal = false;
    do {
      prev = *prev_atomic;
      curr_index = node_atomic::next(prev);
      const bool curr_invalid = curr_index == node_atomic::invalid_next;

      if (curr_invalid) break;

       // global read of the key
      volatile const key_type * const key_ptr = &m_data.get_node(curr_index).value.first;
      const key_type curr_key = *key_ptr;

      const bool curr_less = m_data.key_compare( curr_key, k);
      const bool curr_greater = m_data.key_compare( k, curr_key);
      curr_equal = !curr_less && !curr_greater;

      if (!curr_less) break;

      prev_atomic = & m_data.get_node(node_atomic::next(prev)).atomic.value;
    } while (true);
  }

//private: // private members
public:
  map_data m_data;
};


template <   typename Key
           , typename T
           , typename Device
           , typename Compare
           , typename Hash
        >
class unordered_map< const Key, const T, Device, Compare, Hash>
{
public: // public types and constants
  typedef Impl::UnorderedMap::map_data<const Key, const T,Device,Compare,Hash> map_data;
  typedef Impl::UnorderedMap::node_atomic node_atomic;

  typedef typename map_data::device_type device_type;
  typedef typename map_data::compare_type compare_type;
  typedef typename map_data::hash_type hash_type;
  typedef typename map_data::key_type key_type;
  typedef typename map_data::mapped_type mapped_type;
  typedef typename map_data::value_type value_type;
  typedef typename map_data::pointer pointer;
  typedef typename map_data::const_pointer const_pointer;
  typedef typename map_data::node_type node_type;
  typedef typename map_data::size_type size_type;

public: //public member functions

  template <typename UMap>
  unordered_map(  const UMap & umap )
    : m_data( umap.m_data )
  {}

  void check_sanity() const
  { m_data.check_sanity(); }

  uint32_t size() const
  {  return m_data.size(); }

  uint32_t unused() const
  {  return m_data.unused(); }

  uint32_t pending_delete() const
  {  return m_data.pending_delete(); }

  uint32_t failed_inserts() const
  { return m_data.failed_inserts(); }

  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return m_data.capacity(); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_capacity() const
  { return m_data.hash_capacity(); }

  void remove_pending_delete() const
  {  return m_data.remove_pending_delete_keys(); }

  KOKKOS_INLINE_FUNCTION
  const_pointer find( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    return (node_index != node_atomic::invalid_next) ? &m_data.get_node(node_index).value : NULL;
  }

  KOKKOS_INLINE_FUNCTION
  const_pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i < m_data.capacity();
    const bool used_node  = node_atomic::state(m_data.get_node(i).atomic) == Impl::UnorderedMap::USED;

    return valid_range && used_node ? &m_data.get_node(i).value : NULL;
  }

//private: // private members
public:
  map_data m_data;

};

template <   typename Key
           , typename Device
           , typename Compare
           , typename Hash
        >
class unordered_map< const Key, void, Device, Compare, Hash>
{
public: // public types and constants
  typedef Impl::UnorderedMap::map_data<const Key, void,Device,Compare,Hash> map_data;
  typedef Impl::UnorderedMap::node_atomic node_atomic;

  typedef typename map_data::device_type device_type;
  typedef typename map_data::compare_type compare_type;
  typedef typename map_data::hash_type hash_type;
  typedef typename map_data::key_type key_type;
  typedef typename map_data::mapped_type mapped_type;
  typedef typename map_data::value_type value_type;
  typedef typename map_data::pointer pointer;
  typedef typename map_data::const_pointer const_pointer;
  typedef typename map_data::node_type node_type;
  typedef typename map_data::size_type size_type;

public: //public member functions

  template <typename UMap>
  unordered_map(  const UMap & umap )
    : m_data( umap.m_data )
  {}

  void check_sanity() const
  { m_data.check_sanity(); }

  uint32_t size() const
  {  return m_data.size(); }

  uint32_t unused() const
  {  return m_data.unused(); }

  uint32_t pending_delete() const
  {  return m_data.pending_delete(); }

  uint32_t failed_inserts() const
  { return m_data.failed_inserts(); }

  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return m_data.capacity(); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_capacity() const
  { return m_data.hash_capacity(); }

  void remove_pending_delete() const
  {  return m_data.remove_pending_delete_keys(); }

  KOKKOS_INLINE_FUNCTION
  const_pointer find( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    return (node_index != node_atomic::invalid_next) ? &m_data.get_node(node_index).value : NULL;
  }

  KOKKOS_INLINE_FUNCTION
  const_pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i < m_data.capacity();
    const bool used_node  = node_atomic::state(m_data.get_node(i).atomic) == Impl::UnorderedMap::USED;

    return valid_range && used_node ? &m_data.get_node(i).value : NULL;
  }

//private: // private members
public:
  map_data m_data;
};

} // namespace Kokkos

#endif //KOKKOS_UNORDERED_MAP_HPP
