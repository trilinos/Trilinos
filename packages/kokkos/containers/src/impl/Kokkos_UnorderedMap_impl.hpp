#ifndef KOKKOS_UNORDERED_MAP_IMPL_HPP
#define KOKKOS_UNORDERED_MAP_IMPL_HPP

#include <Kokkos_Functional.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_View.hpp>

#include <stdexcept>
#include <string>
#include <stdint.h>
#include <iostream>
#include <sstream>
#include <cstdio>

namespace Kokkos { namespace Impl {

inline uint32_t find_hash_size(uint32_t size)
{
  // these primes try to preserve randomness of hash
  static const uint32_t primes [] = {
        3, 7, 13, 23, 53, 97, 193, 389, 769, 1543
      , 2237, 2423, 2617, 2797, 2999, 3167, 3359, 3539
      , 3727, 3911, 4441 , 4787 , 5119 , 5471 , 5801 , 6143 , 6521 , 6827
      , 7177 , 7517 , 7853 , 8887 , 9587 , 10243 , 10937 , 11617 , 12289
      , 12967 , 13649 , 14341 , 15013 , 15727
      , 17749 , 19121 , 20479 , 21859 , 23209 , 24593 , 25939 , 27329
      , 28669 , 30047 , 31469 , 35507 , 38231 , 40961 , 43711 , 46439
      , 49157 , 51893 , 54617 , 57347 , 60077 , 62801 , 70583 , 75619
      , 80669 , 85703 , 90749 , 95783 , 100823 , 105871 , 110909 , 115963
      , 120997 , 126031 , 141157 , 151237 , 161323 , 171401 , 181499 , 191579
      , 201653 , 211741 , 221813 , 231893 , 241979 , 252079
      , 282311 , 302483 , 322649 , 342803 , 362969 , 383143 , 403301 , 423457
      , 443629 , 463787 , 483953 , 504121 , 564617 , 604949 , 645313 , 685609
      , 725939 , 766273 , 806609 , 846931 , 887261 , 927587 , 967919 , 1008239
      , 1123477 , 1198397 , 1273289 , 1348177 , 1423067 , 1497983 , 1572869
      , 1647761 , 1722667 , 1797581 , 1872461 , 1947359 , 2022253
      , 2246953 , 2396759 , 2546543 , 2696363 , 2846161 , 2995973 , 3145739
      , 3295541 , 3445357 , 3595117 , 3744941 , 3894707 , 4044503
      , 4493921 , 4793501 , 5093089 , 5392679 , 5692279 , 5991883 , 6291469
      , 6591059 , 6890641 , 7190243 , 7489829 , 7789447 , 8089033
      , 8987807 , 9586981 , 10186177 , 10785371 , 11384539 , 11983729
      , 12582917 , 13182109 , 13781291 , 14380469 , 14979667 , 15578861
      , 16178053 , 17895707 , 19014187 , 20132683 , 21251141 , 22369661
      , 23488103 , 24606583 , 25725083 , 26843549 , 27962027 , 29080529
      , 30198989 , 31317469 , 32435981 , 35791397 , 38028379 , 40265327
      , 42502283 , 44739259 , 46976221 , 49213237 , 51450131 , 53687099
      , 55924061 , 58161041 , 60397993 , 62634959 , 64871921
      , 71582857 , 76056727 , 80530643 , 85004567 , 89478503 , 93952427
      , 98426347 , 102900263 , 107374217 , 111848111 , 116322053 , 120795971
      , 125269877 , 129743807 , 143165587 , 152113427 , 161061283 , 170009141
      , 178956983 , 187904819 , 196852693 , 205800547 , 214748383 , 223696237
      , 232644089 , 241591943 , 250539763 , 259487603 , 268435399
  };


  const size_t num_primes = sizeof(primes)/sizeof(uint32_t);

  uint32_t hsize = primes[num_primes-1] ;
  for (size_t i = 0; i < num_primes; ++i) {
    if (size <= primes[i]) {
      hsize = primes[i];
      break;
    }
  }
  return hsize;
}

struct unordered_map_node_state
{
  enum type
  {
      UNUSED = 0          // not used in a list
    , USED = 1            // used in a list
    , PENDING_INSERT = 2  // not used in a list, but reserved by a thread for inserting
    , PENDING_DELETE = 3  // node in the list is marked deleted
    , INVALID = 4         // the 0th node in the node view is set to invalid
    , NUM_STATES = 5
  };
};

struct unordered_map_node_state_counts
{
  KOKKOS_INLINE_FUNCTION
  unordered_map_node_state_counts() : value() {}

  uint32_t value[unordered_map_node_state::NUM_STATES];
};

struct unordered_map_hash_list_sanity_type
{
  KOKKOS_INLINE_FUNCTION
  unordered_map_hash_list_sanity_type()
    : duplicate_keys_errors(0)
    , unordered_list_errors(0)
    , incorrect_hash_index_errors(0)
  {}

  uint32_t duplicate_keys_errors;
  uint32_t unordered_list_errors;
  uint32_t incorrect_hash_index_errors;
};

struct unordered_map_sanity_type
{
  unordered_map_sanity_type()
    : state_count()
    , hash_list()
    , free_node_count(0)
    , total_errors(0)
  {}

  unordered_map_node_state_counts state_count;
  unordered_map_hash_list_sanity_type hash_list;
  uint32_t free_node_count;
  uint32_t total_errors;
};


template <typename ValueType>
struct unordered_map_node
{
  typedef ValueType value_type;
  typedef unordered_map_node_state::type node_state;

  enum { word_mask = 0xFFFFFFFFu };
  enum { word_shift = 32u };

  KOKKOS_INLINE_FUNCTION
  static uint32_t next(uint64_t v)
  { return static_cast<uint32_t>(v & word_mask); }

  KOKKOS_INLINE_FUNCTION
  static node_state state(uint64_t v)
  { return static_cast<node_state>((v >> word_shift)); }

  KOKKOS_INLINE_FUNCTION
  static uint64_t make_atomic( uint32_t n, node_state s)
  { return (static_cast<uint64_t>(s) << word_shift) | static_cast<uint64_t>(n); }

  // contruct a new value at the current node
  KOKKOS_INLINE_FUNCTION
  void construct_value( const value_type & v )
  { new (&value) value_type(v); }

  // destruct the value at the current node
  KOKKOS_INLINE_FUNCTION
  void destruct_value()
  { value.~value_type(); }

  uint64_t   atomic;
  value_type value;
};


template <class MapData>
struct unordered_map_count_node_states_functor
{
  typedef typename MapData::device_type device_type;
  typedef typename device_type::size_type size_type;
  typedef typename MapData::node_type node_type;
  typedef unordered_map_node_state_counts value_type;

  MapData  map;

  unordered_map_count_node_states_functor(MapData arg_map, value_type & value)
    : map(arg_map)
  {
    parallel_reduce( map.nodes.size(), *this, value);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & dst)
  {
    dst = value_type();
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & dst, const volatile value_type & src)
  {
    for (int i=0; i<unordered_map_node_state::NUM_STATES; ++i) {
      dst.value[i] += src.value[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & dst) const
  {
    unordered_map_node_state::type index = node_type::state(map.nodes[i].atomic);
    ++dst.value[ index < unordered_map_node_state::NUM_STATES ? index : unordered_map_node_state::INVALID];
  }
};


template <class MapData>
struct unordered_map_check_hash_list_functor
{
  typedef typename MapData::device_type device_type;
  typedef typename device_type::size_type size_type;
  typedef typename MapData::node_type node_type;
  typedef unordered_map_hash_list_sanity_type value_type;

  MapData map;

  unordered_map_check_hash_list_functor(MapData arg_map, value_type & value)
    : map(arg_map)
  {
    parallel_reduce( map.hashes.size(), *this, value);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & dst)
  {
    dst.duplicate_keys_errors       = 0;
    dst.unordered_list_errors       = 0;
    dst.incorrect_hash_index_errors = 0;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & dst, const volatile value_type & src)
  {
    dst.duplicate_keys_errors       += src.duplicate_keys_errors;
    dst.unordered_list_errors       += src.unordered_list_errors;
    dst.incorrect_hash_index_errors += src.incorrect_hash_index_errors;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & errors) const
  {
    const uint64_t * prev_atomic = &map.hashes[i];
    const uint64_t * curr_atomic = &map.nodes[node_type::next(*prev_atomic)].atomic;
    const uint64_t * next_atomic = &map.nodes[node_type::next(*curr_atomic)].atomic;

    uint32_t incorrect_hash_index_errors = 0;
    uint32_t duplicate_keys_errors = 0;
    uint32_t unordered_list_errors = 0;

    //traverse the list
    while ( node_type::state(*curr_atomic) != unordered_map_node_state::INVALID) {
      const uint32_t curr_index = node_type::next(*prev_atomic);
      const uint32_t next_index = node_type::next(*curr_atomic);

      //check that the key hashes to this index
      const uint32_t hash_value = map.key_hash(map.nodes[curr_index].value.first);
      const uint32_t hash_index = hash_value%map.hashes.size();

      if ( static_cast<uint32_t>(i) != hash_index) {
        ++incorrect_hash_index_errors;
      }

      if (next_index != 0u) {
        //check that the list is ordered and has no duplicates
        const bool key_less = map.key_compare( map.nodes[curr_index].value.first, map.nodes[next_index].value.first );
        const bool key_greater = map.key_compare( map.nodes[next_index].value.first, map.nodes[curr_index].value.first );
        const bool key_equal = !key_less && !key_greater;

        if (key_equal) {
          ++duplicate_keys_errors;
        }
        else if (key_greater) {
          ++unordered_list_errors;
        }
      }

      prev_atomic = curr_atomic;
      curr_atomic = next_atomic;
      next_atomic = &map.nodes[node_type::next(*curr_atomic)].atomic;
    }

    errors.incorrect_hash_index_errors += incorrect_hash_index_errors;
    errors.duplicate_keys_errors += duplicate_keys_errors;
    errors.unordered_list_errors += unordered_list_errors;
  }
};

template <class MapData>
struct unordered_map_remove_pending_delete_keys_functor
{
  typedef typename MapData::device_type device_type;
  typedef typename device_type::size_type size_type;
  typedef typename MapData::node_type node_type;

  node_type * nodes;
  uint64_t  * hashes;

  unordered_map_remove_pending_delete_keys_functor( MapData arg_map )
    : nodes( const_cast<node_type *>(arg_map.nodes.ptr_on_device()) )
    , hashes( const_cast<uint64_t *>(arg_map.hashes.ptr_on_device()) )
  {
    parallel_for( arg_map.hashes.size(), *this);
    device_type::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i) const
  {
    volatile uint64_t * prev_atomic = &hashes[i];
    volatile uint64_t * curr_atomic = &nodes[ node_type::next(*prev_atomic)].atomic;

    while (node_type::state(*curr_atomic) != unordered_map_node_state::INVALID) {
      uint64_t prev = *prev_atomic;
      uint64_t curr = *curr_atomic;
      if (node_type::state(curr) == unordered_map_node_state::PENDING_DELETE) {
        //remove the node
        *prev_atomic = node_type::make_atomic( node_type::next(curr), node_type::state(prev) );
        *curr_atomic = node_type::make_atomic( 0, unordered_map_node_state::UNUSED );

        curr_atomic = &nodes[ node_type::next(*prev_atomic)].atomic;
      }
      else {
        prev_atomic = curr_atomic;
        curr_atomic = &nodes[ node_type::next(*prev_atomic)].atomic;
      }
    }
  }
};

template <typename Key, typename T, typename Device, typename Compare, typename Hash>
struct unordered_map_data
{
  typedef unordered_map_data<Key,T,Device,Compare,Hash> self_type;

  typedef typename remove_const<Key>::type key_type;
  typedef typename add_const<Key>::type const_key_type;

  typedef typename remove_const<T>::type mapped_type;
  typedef typename add_const<T>::type const_mapped_type;

  typedef Device device_type;
  typedef Compare compare_type;
  typedef Hash hash_type;

  typedef unordered_map_data< key_type, mapped_type, Device, Compare, Hash>              insertable_map_type;
  typedef unordered_map_data< const_key_type, mapped_type, Device, Compare, Hash>        modifiable_map_type;
  typedef unordered_map_data< const_key_type, const_mapped_type, Device, Compare, Hash>  const_map_type;

  static const bool has_const_key_type = is_const<Key>::value;
  static const bool has_void_mapped_type = is_same<T,void>::value;
  static const bool has_const_mapped_type = has_void_mapped_type || is_const<T>::value;
  static const bool is_const_map = has_const_key_type && has_const_mapped_type;


  typedef pair<const_key_type, mapped_type> value_type;

  typedef typename if_c< is_const_map, value_type const *, value_type *>::type pointer;
  typedef value_type const * const_pointer;

  typedef unordered_map_node<value_type> node_type;

  typedef uint32_t size_type;

  typedef typename if_c<   has_const_key_type
                         , View< const uint64_t *, device_type, MemoryTraits<RandomRead> >
                         , View< uint64_t *, device_type >
                       >::type uint64_t_view;

  typedef typename if_c<   is_const_map
                         , View< const node_type *, device_type, MemoryTraits<RandomRead> >
                         , View< node_type *, device_type >
                       >::type node_view;


  typedef View< uint32_t, device_type > scalar_view;
  typedef typename scalar_view::HostMirror host_scalar_view;

  unordered_map_data(  uint32_t num_nodes
                     , compare_type compare
                     , hash_type hash
                    )
    : nodes("unordered_map_nodes", num_nodes+1)
    , hashes("unordered_map_hashes", find_hash_size(nodes.size()) )
    , num_failed_inserts("unordered_map_num_failed_inserts")
    , key_compare(compare)
    , key_hash(hash)
  {
    uint64_t * first_atomic = &nodes.ptr_on_device()->atomic;
    uint64_t atomic = node_type::make_atomic(0, unordered_map_node_state::INVALID);

    typedef Kokkos::DeepCopy< typename device_type::memory_space, Kokkos::HostSpace > deep_copy;
    deep_copy(first_atomic, &atomic, sizeof(uint64_t));
  }

  template <typename MMapType>
  KOKKOS_INLINE_FUNCTION
  unordered_map_data( const MMapType & m)
    : nodes(m.nodes)
    , hashes(m.hashes)
    , num_failed_inserts(m.num_failed_inserts)
    , key_compare(m.key_compare)
    , key_hash(m.key_hash)
  {}

  template <typename MMapType>
  KOKKOS_INLINE_FUNCTION
  unordered_map_data & operator=( const MMapType & m)
  {
    nodes = m.nodes;
    hashes = m.hashes;
    num_failed_inserts = m.num_failed_inserts;
    key_compare = m.key_compare;
    key_hash = m.key_hash;

    return *this;
  }

  unordered_map_node_state_counts get_node_states() const
  {
    unordered_map_node_state_counts result;
    unordered_map_count_node_states_functor<const_map_type>(*this, result);
    device_type::fence();
    return result;
  }

  unordered_map_hash_list_sanity_type check_hash_sanity() const
  {
    unordered_map_hash_list_sanity_type result;
    unordered_map_check_hash_list_functor<const_map_type>(*this, result);
    device_type::fence();
    return result;
  }

  void check_sanity() const
  {
    unordered_map_sanity_type result;

    unordered_map_count_node_states_functor<const_map_type>(*this, result.state_count);
    unordered_map_check_hash_list_functor<const_map_type>(*this, result.hash_list);

    device_type::fence();

    std::ostringstream out;

    const uint32_t failed = failed_inserts();

    if (failed > 0u) {
      out << "Error: " << failed << " failed insertions\n";
      result.total_errors+=failed;
    }


    if (result.hash_list.duplicate_keys_errors > 0u) {
      out << "Error: found " << result.hash_list.duplicate_keys_errors << " duplicate keys found in lists\n";
      result.total_errors+=result.hash_list.duplicate_keys_errors;
    }

    if (result.hash_list.unordered_list_errors > 0u) {
      out << "Error: found " << result.hash_list.unordered_list_errors << " unsorted lists\n";
      result.total_errors+=result.hash_list.unordered_list_errors;
    }

    if (result.hash_list.incorrect_hash_index_errors > 0u) {
      out << "Error: found " << result.hash_list.incorrect_hash_index_errors << " keys incorrectly hashed\n";
      result.total_errors+=result.hash_list.incorrect_hash_index_errors;
    }

    // state_count[ INVALID ] == 1
    if (result.state_count.value[ unordered_map_node_state::INVALID ] != 1u) {
      out << "Error: found " << result.state_count.value[ unordered_map_node_state::INVALID] << " invalid nodes (should always be 1)\n";
      ++result.total_errors;
    }

    // state_count[ PENDING_INSERT ] == 0
    if (result.state_count.value[ unordered_map_node_state::PENDING_INSERT ] > 0u) {
      out << "Error: found " << result.state_count.value[ unordered_map_node_state::PENDING_INSERT] << " pending insert nodes (should always be 0)\n";
      ++result.total_errors;
    }

    if (result.total_errors > 0u) {
      out << "  Details:\n";
      out << "    Node States: ";
      out << " {UNUSED: " << result.state_count.value[ unordered_map_node_state::UNUSED ] << "}";
      out << " {USED: " << result.state_count.value[ unordered_map_node_state::USED ] << "}";
      out << " {PENDING_DELETE: " << result.state_count.value[ unordered_map_node_state::PENDING_DELETE ] << "}";
      out << " {PENDING_INSERT: " << result.state_count.value[ unordered_map_node_state::PENDING_INSERT ] << "}";
      out << " {INVALID: " << result.state_count.value[ unordered_map_node_state::INVALID ] << "}\n";
      out << "   Errors: ";
      out << " {duplicate_keys: " << result.hash_list.duplicate_keys_errors << "}";
      out << " {unordered_lists: " << result.hash_list.unordered_list_errors << "}";
      out << " {incorrect hashes: " << result.hash_list.incorrect_hash_index_errors << "}";
      out << " {failed inserts: " << failed << "}\n";
      out << "   TOTAL Errors: " << result.total_errors;
      out << std::endl;

      throw std::runtime_error( out.str() );
    }
  }


  void remove_pending_delete_keys() const
  {
    unordered_map_remove_pending_delete_keys_functor<self_type> remove_keys(*this);
  }

  uint32_t failed_inserts() const
  {
    host_scalar_view tmp = create_mirror_view(num_failed_inserts);
    deep_copy(tmp, num_failed_inserts);
    return *tmp;
  }

  KOKKOS_INLINE_FUNCTION
  uint32_t find_node_index( const key_type & k) const
  {
    const uint32_t hash_value = key_hash(k);
    const uint32_t hash_index = hash_value % hashes.size();

    uint64_t prev = hashes[hash_index];

    uint32_t index = 0;
    do {
      const uint32_t curr_index = node_type::next(prev);

      if ( curr_index != 0u ) {
        const node_type & curr_node = nodes[curr_index];
        const uint64_t curr = nodes[curr_index].atomic;

        const bool curr_greater = key_compare( k, curr_node.value.first);
        const bool curr_less =  key_compare( curr_node.value.first, k);
        const bool curr_equal = !curr_less && !curr_greater;

        if (curr_greater) {
          index = 0u;
          break;
        } else if (curr_equal) {
          // return existing node
          index = curr_index;
          break;
        }
        else {
          // Current is less -- advance to next node
          prev = curr;
        }
      }
      else {
        break;
      }
    } while (true);

    return index;
  }

  // Data members
  node_view     nodes;
  uint64_t_view hashes;
  scalar_view   num_failed_inserts;
  compare_type  key_compare;
  hash_type     key_hash;
};












}} // namespace Kokkos::Impl

#endif //KOKKOS_UNORDERED_MAP_IMPL_HPP

