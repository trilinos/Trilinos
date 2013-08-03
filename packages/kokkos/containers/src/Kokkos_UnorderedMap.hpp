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
  typedef Impl::unordered_map_data<Key,T,Device,Compare,Hash> impl_data_type;

  typedef typename impl_data_type::device_type device_type;
  typedef typename impl_data_type::compare_type compare_type;
  typedef typename impl_data_type::hash_type hash_type;
  typedef typename impl_data_type::key_type key_type;
  typedef typename impl_data_type::mapped_type mapped_type;
  typedef typename impl_data_type::value_type value_type;
  typedef typename impl_data_type::pointer pointer;
  typedef typename impl_data_type::const_pointer const_pointer;
  typedef typename impl_data_type::node_type node_type;
  typedef typename impl_data_type::insert_result insert_result;
  typedef typename impl_data_type::size_type size_type;

  typedef typename impl_data_type::uint64_t_view uint64_t_view;
  typedef typename impl_data_type::scalar_view scalar_view;
  typedef typename impl_data_type::host_scalar_view host_scalar_view;

  typedef Impl::unordered_map_node_state node_state;

private:

  typedef typename impl_data_type::insert_mapped_type insert_mapped_type;

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
  {  return m_data.get_node_states().value[ node_state::USED ]; }

  uint32_t unused() const
  {  return m_data.get_node_states().value[ node_state::UNUSED ]; }

  uint32_t marked_deleted() const
  {  return m_data.get_node_states().value[ node_state::MARKED_DELETED ]; }

  uint32_t num_failed_inserts() const
  { return m_data.get_num_failed_inserts(); }


  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return static_cast<uint32_t>(m_data.nodes.size()-1); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_length() const
  { return static_cast<uint32_t>(m_data.hashes.size()); }


  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  void remove_keys_marked_deleted() const
  {  return m_data.remove_keys_marked_deleted(); }


#if defined( __CUDACC__ )
    __device__ inline
#else
    KOKKOS_INLINE_FUNCTION
#endif
  insert_result insert(const key_type & k, const insert_mapped_type & v = insert_mapped_type()) const
  {
    insert_result result(INSERT_FAILED,NULL);

    const uint32_t hash_value = m_data.key_hash(k);
    const uint32_t hash_index = hash_value % m_data.hashes.size();

    volatile uint64_t * prev_atomic = & m_data.hashes[hash_index];

    uint32_t node_index = 0u;


    do {
      uint64_t prev = *prev_atomic;

      uint32_t curr_index = node_type::next(prev);

      const bool curr_invalid = curr_index == 0u;
      bool curr_greater = false;
      bool curr_less = false;
      bool curr_equal = false;

      if (!curr_invalid) {
        // global read of the key
        //volatile node_type * curr_node = &m_data.nodes[curr_index];
        volatile const key_type * const key_ptr = &m_data.nodes[curr_index].value.first;
        const key_type curr_key = *key_ptr;
        //const key_type curr_key = m_data.nodes[curr_index].value.first;

        curr_greater = m_data.key_compare( k, curr_key);
        curr_less = m_data.key_compare( curr_key, k);
        curr_equal = !curr_less && !curr_greater;
      }

      const bool insert_here = curr_invalid || curr_greater;

      // Can insert node
      if (insert_here) {
        if (node_index == 0u) {
          node_index = get_node_index(hash_value);
          if (node_index == 0u) {
            // unable to obtain a free node
            break;
          }
        }

        // this thread has unique control of the node
        // so can construct the value and set up the state and next index
        node_type & n = m_data.nodes[node_index];
        n.destruct_value();
        n.construct_value(value_type(k,v));
        n.atomic = node_type::make_atomic( curr_index, node_state::USED);

        uint64_t new_atomic = node_type::make_atomic( node_index, node_type::state(prev));

#if defined( __CUDACC__ )
        __threadfence();
#endif

        const bool ok = atomic_compare_exchange_strong( prev_atomic, prev, new_atomic);
        if ( ok ) {
          // successfully inserted the node
          result = insert_result(INSERT_SUCCESS, &n.value);
        }

        // atomic update failed
        // set the node to UNUSED
        // try again from prev node
      }
      else if (curr_equal) {
        if (node_index != 0u) {
          // release any node that was claimed by this thread
          m_data.nodes[node_index].atomic = node_type::make_atomic(0, node_state::UNUSED);
        }
        // Node already exist
        // DO NOT resurrect marked_deleted node
        result = insert_result(INSERT_EXISTING, &m_data.nodes[curr_index].value);
      }
      // Current is less -- advance to next node
      else {
        prev_atomic = & m_data.nodes[node_type::next(prev)].atomic;
      }
    } while (result.second == NULL);
    return result;
  }

  // TODO protect with enable_if
  KOKKOS_INLINE_FUNCTION
  void mark_for_deletion( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    if (node_index != 0u) {
      volatile node_type * const nodes = &m_data.nodes[0];
      volatile uint64_t * atomic = &((nodes+node_index)->atomic);
      uint64_t value = *atomic;

      while ( node_type::state(value) != node_state::MARKED_DELETED) {
        uint64_t new_value = node_type::make_atomic( node_type::next(value), node_state::MARKED_DELETED);
        atomic_compare_exchange_strong( atomic, value, new_value);
        value = *atomic;
      }
    }
  }

  // TODO protect with enable_if
  KOKKOS_INLINE_FUNCTION
  void mark_for_deletion( const_pointer p ) const
  {
    if (p) mark_for_deletion(p->first);
  }

  KOKKOS_INLINE_FUNCTION
  pointer find( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    return (node_index != 0u) ? &m_data.nodes[node_index].value : NULL;
  }


  KOKKOS_INLINE_FUNCTION
  pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i+1 < m_data.nodes.size();
    const bool used_node  = node_type::state(m_data.nodes[i+1].atomic) == node_state::USED;

    return valid_range && used_node ? &m_data.nodes[i+1].value : NULL;
  }

private: // private member functions

#if defined( __CUDACC__ )
    __device__ inline
#else
    KOKKOS_INLINE_FUNCTION
#endif
  uint32_t get_node_index(uint32_t hash_value) const
  {
    if (*m_data.num_failed_inserts == 0u) {
      const uint32_t num_nodes = m_data.nodes.size();
      const uint32_t start_node = hash_value % num_nodes;
      const uint32_t end_node = start_node + num_nodes;

      for (uint32_t i = start_node; i < end_node; ++i) {
        const uint32_t n = i % num_nodes;
        volatile uint64_t * atomic = &m_data.nodes[n].atomic;
        uint64_t value = *atomic;
        if (    (node_type::state(value) == node_state::UNUSED)
            && atomic_compare_exchange_strong(atomic, value, node_type::make_atomic(0,node_state::PENDING_INSERT)) )
        {
          return n;
        }
      }
    }
    // unable to get a free node -- insert failed
    volatile uint32_t * addr = &*m_data.num_failed_inserts;
    atomic_fetch_add( addr, 1u);
    return 0;
  }


//private: // private members
public:
  impl_data_type m_data;

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
  typedef Impl::unordered_map_data<const Key, T,Device,Compare,Hash> impl_data_type;

  typedef typename impl_data_type::device_type device_type;
  typedef typename impl_data_type::compare_type compare_type;
  typedef typename impl_data_type::hash_type hash_type;
  typedef typename impl_data_type::key_type key_type;
  typedef typename impl_data_type::mapped_type mapped_type;
  typedef typename impl_data_type::value_type value_type;
  typedef typename impl_data_type::pointer pointer;
  typedef typename impl_data_type::const_pointer const_pointer;
  typedef typename impl_data_type::node_type node_type;
  typedef typename impl_data_type::insert_result insert_result;
  typedef typename impl_data_type::size_type size_type;

  typedef typename impl_data_type::uint64_t_view uint64_t_view;
  typedef typename impl_data_type::scalar_view scalar_view;
  typedef typename impl_data_type::host_scalar_view host_scalar_view;

  typedef Impl::unordered_map_node_state node_state;

private:

  typedef typename impl_data_type::insert_mapped_type insert_mapped_type;

public: //public member functions

  template <typename UMap>
  unordered_map(  const UMap & umap )
    : m_data( umap.m_data )
  {}

  void check_sanity() const
  { m_data.check_sanity(); }


  uint32_t size() const
  {  return m_data.get_node_states().value[ node_state::USED ]; }

  uint32_t unused() const
  {  return m_data.get_node_states().value[ node_state::UNUSED ]; }

  uint32_t marked_deleted() const
  {  return m_data.get_node_states().value[ node_state::MARKED_DELETED ]; }

  uint32_t num_failed_inserts() const
  { return m_data.get_num_failed_inserts(); }


  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return static_cast<uint32_t>(m_data.nodes.size()-1); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_length() const
  { return static_cast<uint32_t>(m_data.hashes.size()); }


  KOKKOS_INLINE_FUNCTION
  pointer find( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    return (node_index != 0u) ? &m_data.nodes[node_index].value : NULL;
  }


  KOKKOS_INLINE_FUNCTION
  pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i+1 < m_data.nodes.size();
    const bool used_node  = node_type::state(m_data.nodes[i+1].atomic) == node_state::USED;

    return valid_range && used_node ? &m_data.nodes[i+1].value : NULL;
  }

//private: // private members
public:
  impl_data_type m_data;

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
  typedef Impl::unordered_map_data<const Key, const T,Device,Compare,Hash> impl_data_type;

  typedef typename impl_data_type::device_type device_type;
  typedef typename impl_data_type::compare_type compare_type;
  typedef typename impl_data_type::hash_type hash_type;
  typedef typename impl_data_type::key_type key_type;
  typedef typename impl_data_type::mapped_type mapped_type;
  typedef typename impl_data_type::value_type value_type;
  typedef typename impl_data_type::pointer pointer;
  typedef typename impl_data_type::const_pointer const_pointer;
  typedef typename impl_data_type::node_type node_type;
  typedef typename impl_data_type::insert_result insert_result;
  typedef typename impl_data_type::size_type size_type;

  typedef typename impl_data_type::uint64_t_view uint64_t_view;
  typedef typename impl_data_type::scalar_view scalar_view;
  typedef typename impl_data_type::host_scalar_view host_scalar_view;

  typedef Impl::unordered_map_node_state node_state;

private:

  typedef typename impl_data_type::insert_mapped_type insert_mapped_type;

public: //public member functions

  template <typename UMap>
  unordered_map(  const UMap & umap )
    : m_data( umap.m_data )
  {}

  void check_sanity() const
  { m_data.check_sanity(); }


  uint32_t size() const
  {  return m_data.get_node_states().value[ node_state::USED ]; }

  uint32_t unused() const
  {  return m_data.get_node_states().value[ node_state::UNUSED ]; }

  uint32_t marked_deleted() const
  {  return m_data.get_node_states().value[ node_state::MARKED_DELETED ]; }

  uint32_t num_failed_inserts() const
  { return m_data.get_num_failed_inserts(); }


  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return static_cast<uint32_t>(m_data.nodes.size()-1); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_length() const
  { return static_cast<uint32_t>(m_data.hashes.size()); }


  KOKKOS_INLINE_FUNCTION
  const_pointer find( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    return (node_index != 0u) ? &m_data.nodes[node_index].value : NULL;
  }


  KOKKOS_INLINE_FUNCTION
  const_pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i+1 < m_data.nodes.size();
    const bool used_node  = node_type::state(m_data.nodes[i+1].atomic) == node_state::USED;

    return valid_range && used_node ? &m_data.nodes[i+1].value : NULL;
  }

//private: // private members
public:
  impl_data_type m_data;

};

template <   typename Key
           , typename Device
           , typename Compare
           , typename Hash
        >
class unordered_map< const Key, void, Device, Compare, Hash>
{
public: // public types and constants
  typedef Impl::unordered_map_data<const Key, void,Device,Compare,Hash> impl_data_type;

  typedef typename impl_data_type::device_type device_type;
  typedef typename impl_data_type::compare_type compare_type;
  typedef typename impl_data_type::hash_type hash_type;
  typedef typename impl_data_type::key_type key_type;
  typedef typename impl_data_type::mapped_type mapped_type;
  typedef typename impl_data_type::value_type value_type;
  typedef typename impl_data_type::pointer pointer;
  typedef typename impl_data_type::const_pointer const_pointer;
  typedef typename impl_data_type::node_type node_type;
  typedef typename impl_data_type::insert_result insert_result;
  typedef typename impl_data_type::size_type size_type;

  typedef typename impl_data_type::uint64_t_view uint64_t_view;
  typedef typename impl_data_type::scalar_view scalar_view;
  typedef typename impl_data_type::host_scalar_view host_scalar_view;

  typedef Impl::unordered_map_node_state node_state;

private:

  typedef typename impl_data_type::insert_mapped_type insert_mapped_type;

public: //public member functions

  template <typename UMap>
  unordered_map(  const UMap & umap )
    : m_data( umap.m_data )
  {}

  void check_sanity() const
  { m_data.check_sanity(); }


  uint32_t size() const
  {  return m_data.get_node_states().value[ node_state::USED ]; }

  uint32_t unused() const
  {  return m_data.get_node_states().value[ node_state::UNUSED ]; }

  uint32_t marked_deleted() const
  {  return m_data.get_node_states().value[ node_state::MARKED_DELETED ]; }

  uint32_t num_failed_inserts() const
  { return m_data.get_num_failed_inserts(); }


  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return static_cast<uint32_t>(m_data.nodes.size()-1); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_length() const
  { return static_cast<uint32_t>(m_data.hashes.size()); }


  KOKKOS_INLINE_FUNCTION
  const_pointer find( const key_type & k) const
  {
    const uint32_t node_index = m_data.find_node_index(k);
    return (node_index != 0u) ? &m_data.nodes[node_index].value : NULL;
  }


  KOKKOS_INLINE_FUNCTION
  const_pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i+1 < m_data.nodes.size();
    const bool used_node  = node_type::state(m_data.nodes[i+1].atomic) == node_state::USED;

    return valid_range && used_node ? &m_data.nodes[i+1].value : NULL;
  }

//private: // private members
public:
  impl_data_type m_data;

};

} // namespace Kokkos

#endif //KOKKOS_UNORDERED_MAP_HPP
