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
  static const bool has_const_key_type    = Impl::is_const<Key>::value;
  static const bool has_void_mapped_type  = Impl::is_same<T,void>::value;
  static const bool has_const_mapped_type = has_void_mapped_type || Impl::is_const<T>::value;
  static const bool is_const_map          = has_const_key_type && has_const_mapped_type;


  typedef unordered_map<Key,T,Device,Compare,Hash> self_type;

  typedef Key    key_type;
  typedef T      mapped_type;
  typedef Device device_type;

  typedef pair<const key_type, typename Impl::remove_const<mapped_type>::type> value_type;

  typedef       value_type *       pointer;
  typedef const value_type * const_pointer;

  // insert is not allow on an const_map
  typedef pair<unordered_map_insert_state, pointer> insert_result;

  typedef uint32_t size_type;

  typedef Compare compare_type;
  typedef Hash    hash_type;

  typedef Impl::unordered_map_node<value_type> node_type;

  typedef View< Impl::unordered_map_atomic*, device_type > hash_view_type;
  typedef View< uint32_t*, device_type > free_count_view_type;
  typedef View< node_type*, device_type > node_view_type;

private:

  typedef typename Impl::if_c<  has_void_mapped_type
                              , int
                              , mapped_type
                             >::type insert_mapped_type;

  typedef View<uint32_t, device_type> ScalarView;
  typedef typename ScalarView::HostMirror HostScalarView;


public: //public member functions

  unordered_map(  uint32_t arg_num_nodes
                , compare_type compare = compare_type()
                , hash_type hash = hash_type()
                , uint32_t free_node_compression = 32u
               )
    : m_hashes("map_hash", Impl::find_hash_size(arg_num_nodes+1))
    , m_free_counts("map_free_counts", Impl::find_hash_size( arg_num_nodes/free_node_compression))
    , m_nodes("map_nodes", arg_num_nodes + 1)
    , m_num_failed_inserts("map_num_failed_inserts")
    , m_free_block_size( static_cast<uint32_t>((m_nodes.size() + m_free_counts.size() -1u) / m_free_counts.size()) )
    , m_key_compare(compare)
    , m_key_hash(hash)
  {
    Impl::init_unordered_map<self_type> init( m_nodes, m_free_counts, m_free_block_size );
  }

  uint32_t size() const
  {
    uint32_t value = 0;
    Impl::unordered_map_size_functor<self_type> s(m_nodes, value);
    return value;
  }

  uint32_t count_duplicate_keys() const
  {
    uint32_t value = 0;
    Impl::unordered_map_count_duplicate_keys_functor<self_type> s(m_nodes, m_hashes, value);
    return value;
  }

  uint32_t unused() const
  {
    uint32_t value = 0;
    Impl::unordered_map_count_unused_functor<self_type> s(m_free_counts, value);
    return value;
  }

  uint32_t delete_marked_keys() const
  {
    uint32_t value;
    Impl::unordered_map_delete_marked_keys_functor<self_type> s(m_nodes, m_free_counts, m_hashes, m_free_block_size, value);
    return value;
  }

  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return static_cast<uint32_t>(m_nodes.size()-1); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_length() const
  { return static_cast<uint32_t>(m_hashes.size()); }

  KOKKOS_INLINE_FUNCTION
  uint32_t free_length() const
  { return static_cast<uint32_t>(m_free_counts.size()); }

  uint32_t num_failed_inserts() const
  {
    HostScalarView tmp = create_mirror_view(m_num_failed_inserts);
    deep_copy(tmp, m_num_failed_inserts);
    return *tmp;
  }

  // TODO protect with enable_if
  KOKKOS_INLINE_FUNCTION
  insert_result insert(const key_type & k, const insert_mapped_type & v = insert_mapped_type()) const
  {
    const uint32_t hash_value = m_key_hash(k);
    const uint32_t hash_index = hash_value % m_hashes.size();
    Impl::unordered_map_atomic * prev_atomic_ptr = & m_hashes[hash_index];

    insert_result result;
    do {
      uint64_t prev_node_atomic = prev_atomic_ptr->m_atomic;

      node_type & curr_node = m_nodes[Impl::unordered_map_atomic(prev_node_atomic).next()];

      uint64_t  curr_node_atomic = curr_node.atomic.m_atomic;

      const bool curr_invalid = Impl::unordered_map_atomic(curr_node_atomic).invalid();
      const bool curr_greater = !curr_invalid && m_key_compare( k, curr_node.value.first);
      const bool curr_less = !curr_invalid && m_key_compare( curr_node.value.first, k);
      const bool curr_equal = !curr_less && !curr_greater;
      const bool insert_here = curr_greater || curr_invalid;

      // Can insert node
      if (insert_here) {

        uint32_t node_index = get_free_node_index(hash_value);
        if (node_index == 0u) {
          result = insert_result(INSERT_FAILED, NULL);
          break;
        }
        // this thread has unique control of the node
        // so can construct the value and set up the state and next index
        node_type & n = m_nodes[node_index];
        n.destruct_value();
        n.construct_value(value_type(k,v));
        n.atomic.m_value.next = Impl::unordered_map_atomic(prev_node_atomic).next();
        n.atomic.m_value.state = Impl::unordered_map_atomic::USED;

        if ( prev_atomic_ptr->atomic_set_next( Impl::unordered_map_atomic(prev_node_atomic), node_index) ) {
          // successfully inserted the node
          result = insert_result(INSERT_SUCCESS, &n.value);
          break;
        }
        else {
          // atomic update failed -- return node to free list and try again
          return_pending_create_node_to_free_list(node_index);
          continue;
        }
      }
      // Node already exist
      else if (curr_equal) {
        // resurrect marked_deleted node
        while ( Impl::unordered_map_atomic(curr_node_atomic).marked_deleted() ) {
          curr_node.atomic.atomic_set_state(Impl::unordered_map_atomic(curr_node_atomic), Impl::unordered_map_atomic::USED);
          curr_node_atomic = curr_node.atomic.m_atomic;
        }
        result = insert_result(INSERT_EXISTING, &curr_node.value);
        break;
      }
      // Current is less -- advance to next node
      else {
        prev_atomic_ptr = & m_nodes[Impl::unordered_map_atomic(prev_node_atomic).next()].atomic;
      }
    } while (true);

    return result;
  }

  // TODO protect with enable_if
  KOKKOS_INLINE_FUNCTION
  void mark_for_deletion( const key_type & k) const
  {
    node_type & n = find_node(k);

    if (&n != &m_nodes[0]) {
      Impl::unordered_map_atomic atomic = n.atomic;
      while ( !atomic.marked_deleted() ) {
        n.atomic.atomic_set_state( atomic, Impl::unordered_map_atomic::MARKED_DELETED );
        atomic = n.atomic;
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
    node_type & n = find_node(k);
    return (&n != &m_nodes[0]) ? &n.value : NULL;
  }


  KOKKOS_INLINE_FUNCTION
  pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i+1 < m_nodes.size();
    const bool used_node  = m_nodes[i+1].atomic.used();

    return valid_range && used_node ? &m_nodes[i+1].value : NULL;
  }

private: // private member functions

  KOKKOS_INLINE_FUNCTION
  uint32_t get_free_node_index(uint32_t hash_value) const
  {
    const uint32_t free_size = static_cast<uint32_t>(m_free_counts.size());
    uint32_t free_start = hash_value % free_size;
    const uint32_t free_end = free_size + free_start;

    if (*m_num_failed_inserts == 0u) {
      for (; free_start < free_end; ++free_start) {
        const uint32_t free_block = free_start % free_size;
        uint32_t count = m_free_counts[free_block];
        while ( count > 0u) {
          // decrement node count
          if ( atomic_compare_exchange_strong( & m_free_counts[free_block], count, count-1) ) {
            const uint32_t node_start = free_block * m_free_block_size;
            const uint32_t node_end = node_start + m_free_block_size <= m_nodes.size() ? node_start + m_free_block_size : m_nodes.size();
            const uint32_t block_size = node_end - node_start;

            uint32_t i = hash_value % block_size;
            const uint32_t end = i+block_size;

            // must have at least one free node in range
            for (; i < end; ++i) {
              uint32_t node_index = node_start + (i%block_size);
              uint64_t atomic = m_nodes[node_index].atomic.m_atomic;
              while ( Impl::unordered_map_atomic(atomic).unused() ) {
              // try to claim the unused node
                if (m_nodes[node_index].atomic.atomic_set_state(Impl::unordered_map_atomic(atomic), Impl::unordered_map_atomic::PENDING_INSERT)) {
                  return node_index;
                }
                atomic = m_nodes[node_index].atomic.m_atomic;
              }
            }
          }
          count = m_free_counts[free_block];
        }
      }
    }
    volatile uint32_t * addr = &*m_num_failed_inserts;
    atomic_fetch_add( addr, 1u);
    return 0;
  }

  KOKKOS_INLINE_FUNCTION
  void return_pending_create_node_to_free_list(uint64_t node_index) const
  {
    node_type & n = m_nodes[node_index];
    n.atomic.m_value.next = 0;
    n.atomic.m_value.state = Impl::unordered_map_atomic::UNUSED;

    const uint64_t free_block = node_index / m_free_block_size;


    volatile uint32_t * addr = &m_free_counts[free_block];
    atomic_fetch_add( addr, 1u);
  }

  // TODO protect with enable_if
  KOKKOS_INLINE_FUNCTION
  node_type & find_node( const key_type & k) const
  {
    const uint32_t hash_value = m_key_hash(k);
    const uint32_t hash_index = hash_value % m_hashes.size();
    Impl::unordered_map_atomic * prev_atomic_ptr = & m_hashes[hash_index];

    uint32_t index = ~0u;
    do {
      uint64_t prev_node_atomic = prev_atomic_ptr->m_atomic;

      node_type & curr_node = m_nodes[Impl::unordered_map_atomic(prev_node_atomic).next()];

      uint64_t  curr_node_atomic = curr_node.atomic;

      const bool curr_invalid = Impl::unordered_map_atomic(curr_node_atomic).invalid();
      const bool curr_greater = m_key_compare( k, curr_node.value.first);
      const bool curr_less =  m_key_compare( curr_node.value.first, k);
      const bool curr_equal = !(curr_less || curr_greater);
      const bool not_found = curr_invalid || curr_greater;


      if (not_found) {
        index = 0u;
      } else if (curr_equal) {
        // return existing node
        index = Impl::unordered_map_atomic(prev_node_atomic).next();
      }
      else {
        // Current is less -- advance to next node
        prev_atomic_ptr = & m_nodes[Impl::unordered_map_atomic(prev_node_atomic).next()].atomic;
      }
    } while (index == ~0u);

    return m_nodes[index];
  }


//private: // private members
public:

  hash_view_type       m_hashes;
  free_count_view_type m_free_counts;
  node_view_type       m_nodes;

  ScalarView m_num_failed_inserts;

  uint32_t     m_free_block_size;
  compare_type m_key_compare;
  hash_type    m_key_hash;


};


template <   typename Key
           , typename T
           , typename Device
           , typename Compare
           , typename Hash
        >
class unordered_map<const Key, T, Device, Compare, Hash>
{
public: // public types and constants
  static const bool has_const_key_type    = true;
  static const bool has_void_mapped_type  = Impl::is_same<T,void>::value;
  static const bool has_const_mapped_type = has_void_mapped_type || Impl::is_const<T>::value;
  static const bool is_const_map          = has_const_key_type && has_const_mapped_type;


  typedef unordered_map<Key,T,Device,Compare,Hash> self_type;

  typedef Key    key_type;
  typedef T      mapped_type;
  typedef Device device_type;

  typedef pair<const key_type, typename Impl::remove_const<mapped_type>::type> value_type;


  typedef typename Impl::if_c<  is_const_map
                              , const value_type *
                              , value_type *
                            >::type pointer;

  typedef const value_type * const_pointer;

  typedef uint32_t size_type;

  typedef Compare compare_type;
  typedef Hash    hash_type;


  typedef Impl::unordered_map_node<value_type> node_type;

  typedef  View< const Impl::unordered_map_atomic*, device_type, MemoryTraits<RandomRead> > hash_view_type;
  typedef  View< const uint32_t*, device_type, MemoryTraits<RandomRead> > free_count_view_type;

  typedef typename Impl::if_c<   is_const_map
                              , View< const node_type*, device_type, MemoryTraits<RandomRead> >
                              , View< node_type*, device_type >
                            >::type node_view_type;


public: //public member functions

  template < typename Map>
  KOKKOS_INLINE_FUNCTION
  unordered_map( const Map & m)
    : m_hashes(m.m_hashes)
    , m_free_counts(m.m_free_counts)
    , m_nodes(m.m_nodes)
    , m_free_block_size(m.m_free_block_size)
    , m_key_compare(m.m_key_compare)
    , m_key_hash(m.m_key_hash)
  {}

  uint32_t size() const
  {
    uint32_t value = 0;
    Impl::unordered_map_size_functor<self_type> s(m_nodes, value);
    return value;
  }

  uint32_t count_duplicate_keys() const
  {
    uint32_t value = 0;
    Impl::unordered_map_count_duplicate_keys_functor<self_type> s(m_nodes, m_hashes, value);
    return value;
  }

  uint32_t unused() const
  {
    uint32_t value = 0;
    Impl::unordered_map_count_unused_functor<self_type> s(m_free_counts, value);
    return value;
  }

  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  { return static_cast<uint32_t>(m_nodes.size()-1); }

  KOKKOS_INLINE_FUNCTION
  uint32_t hash_length() const
  { return static_cast<uint32_t>(m_hashes.size()); }

  KOKKOS_INLINE_FUNCTION
  uint32_t free_length() const
  { return static_cast<uint32_t>(m_free_counts.size()); }

  KOKKOS_INLINE_FUNCTION
  pointer find( const key_type & k) const
  {
    const node_type & n = find_node(k);
    return (&n != &m_nodes[0]) ? &n.value : NULL;
  }

  KOKKOS_INLINE_FUNCTION
  pointer get_value(uint64_t i) const
  {
    // add one to pass 0th node
    const bool valid_range = i+1 < m_nodes.size();
    const bool used_node  = m_nodes[i+1].atomic.used();

    return valid_range && used_node ? &m_nodes[i+1].value : NULL;
  }

private: // private member functions

  // TODO protect with enable_if
  KOKKOS_INLINE_FUNCTION
  const node_type & find_node( const key_type & k) const
  {
    const uint32_t hash_value = m_key_hash(k);
    const uint32_t hash_index = hash_value % m_hashes.size();
    const Impl::unordered_map_atomic * prev_atomic_ptr = & m_hashes[hash_index];

    uint32_t index = ~0u;
    do {
      uint64_t prev_node_atomic = prev_atomic_ptr->m_atomic;

      const node_type & curr_node = m_nodes[Impl::unordered_map_atomic(prev_node_atomic).next()];

      uint64_t  curr_node_atomic = curr_node.atomic.m_atomic;

      const bool curr_invalid = Impl::unordered_map_atomic(curr_node_atomic).invalid();
      const bool curr_greater = m_key_compare( k, curr_node.value.first);
      const bool curr_less =  m_key_compare( curr_node.value.first, k);
      const bool curr_equal = !(curr_less || curr_greater);
      const bool not_found = curr_invalid || curr_greater;


      if (not_found) {
        index = 0u;
      } else if (curr_equal) {
        // return existing node
        index = Impl::unordered_map_atomic(prev_node_atomic).next();
      }
      else {
        // Current is less -- advance to next node
        prev_atomic_ptr = & m_nodes[Impl::unordered_map_atomic(prev_node_atomic).next()].atomic;
      }
    } while (index == ~0u);

    return m_nodes[index];
  }

//private: // private members
public:

  hash_view_type       m_hashes;
  free_count_view_type m_free_counts;
  node_view_type       m_nodes;

  uint32_t     m_free_block_size;
  compare_type m_key_compare;
  hash_type    m_key_hash;
};


} // namespace Kokkos

#endif //KOKKOS_UNORDERED_MAP_HPP
