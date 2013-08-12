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

namespace Kokkos { namespace Impl { namespace UnorderedMap {

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

enum node_state
{
      UNUSED = 0          // not used in a list
    , USED = 1            // used in a list
    , PENDING_INSERT = 2  // not used in a list, but reserved by a thread for inserting
    , PENDING_DELETE = 3  // node in the list is marked deleted
    , INVALID = 4         // the 0th node in the node view is set to invalid
    , NUM_STATES = 5
};

struct node_atomic
{
  static const uint64_t word_mask = 0x00000000FFFFFFFFu;
  static const uint64_t word_shift = 32u;
  static const uint32_t invalid_next = 0xFFFFFFFFu;

  KOKKOS_FORCEINLINE_FUNCTION
  static uint32_t next(uint64_t v)
  { return static_cast<uint32_t>(v & word_mask); }

  KOKKOS_FORCEINLINE_FUNCTION
  static node_state state(uint64_t v)
  { return static_cast<node_state>((v >> word_shift)); }

  KOKKOS_FORCEINLINE_FUNCTION
  static uint64_t make_atomic( uint32_t n, node_state s)
  { return (static_cast<uint64_t>(s) << word_shift) | static_cast<uint64_t>(n); }

  KOKKOS_FORCEINLINE_FUNCTION
  node_atomic(uint64_t v = make_atomic(invalid_next, UNUSED) )
    : value(v)
  {}

  KOKKOS_FORCEINLINE_FUNCTION
  operator uint64_t() const
  { return value; }

  uint64_t value;
};

template <size_t Size>
struct Align16
{
  static const size_t value = (Size & 15ull);
};

template <typename ValueType, size_t AlignPad = Align16<sizeof(ValueType) + sizeof(uint64_t) >::value >
struct node
{
  typedef ValueType value_type;

  // contruct a new value at the current node
  KOKKOS_FORCEINLINE_FUNCTION
  void construct_value( const value_type & v )
  { new (&value) value_type(v); }

  // destruct the value at the current node
  KOKKOS_FORCEINLINE_FUNCTION
  void destruct_value()
  { value.~value_type(); }

  value_type value;
  uint8_t pad[AlignPad];
  node_atomic atomic;
};

template <typename Node>
struct node_block
{
  typedef Node node_type;

  static const uint32_t size = 32;
  static const uint32_t mask = 31;
  static const uint32_t shift = 5;

  typedef typename StaticAssert< (sizeof(node_type) % 16u == 0u)>::type node_okay;

  KOKKOS_FORCEINLINE_FUNCTION
  node_block()
    : used_count(0)
    , failed_insert_count(0)
    , nodes()
  {}

  int used_count;
  int failed_insert_count;
  node_type nodes[size];
};

template <typename ValueType>
struct node<ValueType, 0u>
{
  typedef ValueType value_type;

  // contruct a new value at the current node
  KOKKOS_FORCEINLINE_FUNCTION
  void construct_value( const value_type & v )
  { new (&value) value_type(v); }

  // destruct the value at the current node
  KOKKOS_FORCEINLINE_FUNCTION
  void destruct_value()
  { value.~value_type(); }

  value_type value;
  node_atomic atomic;
};

struct node_state_counts
{
  KOKKOS_INLINE_FUNCTION
  node_state_counts() : value() {}

  uint32_t value[NUM_STATES];
};

struct hash_list_sanity_type
{
  KOKKOS_INLINE_FUNCTION
  hash_list_sanity_type()
    : duplicate_keys_errors(0)
    , unordered_list_errors(0)
    , incorrect_hash_index_errors(0)
  {}

  uint32_t duplicate_keys_errors;
  uint32_t unordered_list_errors;
  uint32_t incorrect_hash_index_errors;
};

struct sanity_type
{
  sanity_type()
    : state_count()
    , hash_list()
    , free_node_count(0)
    , total_errors(0)
  {}

  node_state_counts state_count;
  hash_list_sanity_type hash_list;
  uint32_t free_node_count;
  uint32_t total_errors;
};




template <class MapData>
struct count_node_states_functor
{
  typedef typename MapData::device_type device_type;
  typedef typename device_type::size_type size_type;
  typedef typename MapData::node_type node_type;
  typedef node_state_counts value_type;

  MapData  map;

  count_node_states_functor(MapData arg_map, value_type & value)
    : map(arg_map)
  {
    parallel_reduce( map.capacity(), *this, value);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & dst)
  {
    dst = value_type();
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & dst, const volatile value_type & src)
  {
    for (int i=0; i<NUM_STATES; ++i) {
      dst.value[i] += src.value[i];
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & dst) const
  {
    node_state index = node_atomic::state(map.get_node(i).atomic);
    ++dst.value[ index < NUM_STATES ? index : INVALID];
  }
};


template <class MapData>
struct check_hash_list_functor
{
  typedef typename MapData::device_type device_type;
  typedef typename device_type::size_type size_type;
  typedef typename MapData::node_type node_type;
  typedef hash_list_sanity_type value_type;

  MapData map;

  check_hash_list_functor(MapData arg_map, value_type & value)
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
    const uint64_t * prev_atomic = &map.hashes[i].value;

    uint32_t incorrect_hash_index_errors = 0;
    uint32_t duplicate_keys_errors = 0;
    uint32_t unordered_list_errors = 0;

    //traverse the list
    while ( node_atomic::next(*prev_atomic) != node_atomic::invalid_next) {
      const uint64_t * curr_atomic = &map.get_node(node_atomic::next(*prev_atomic)).atomic.value;

      const uint32_t curr_index = node_atomic::next(*prev_atomic);
      const uint32_t next_index = node_atomic::next(*curr_atomic);

      //check that the key hashes to this index
      const uint32_t hash_value = map.key_hash(map.get_node(curr_index).value.first);
      const uint32_t hash_index = hash_value%map.hashes.size();

      if ( static_cast<uint32_t>(i) != hash_index) {
        ++incorrect_hash_index_errors;
      }

      if (next_index != node_atomic::invalid_next) {
        //check that the list is ordered and has no duplicates
        const bool key_less = map.key_compare( map.get_node(curr_index).value.first, map.get_node(next_index).value.first );
        const bool key_greater = map.key_compare( map.get_node(next_index).value.first, map.get_node(curr_index).value.first );
        const bool key_equal = !key_less && !key_greater;

        if (key_equal) {
          ++duplicate_keys_errors;
        }
        else if (key_greater) {
          ++unordered_list_errors;
        }
      }

      prev_atomic = curr_atomic;
    }

    errors.incorrect_hash_index_errors += incorrect_hash_index_errors;
    errors.duplicate_keys_errors += duplicate_keys_errors;
    errors.unordered_list_errors += unordered_list_errors;
  }
};

template <class MapData>
struct remove_pending_delete_keys_functor
{
  typedef typename MapData::device_type device_type;
  typedef typename device_type::size_type size_type;
  typedef typename MapData::node_type node_type;
  typedef typename MapData::node_block_type node_block_type;

  node_block_type * node_blocks;
  node_atomic     * hashes;

  remove_pending_delete_keys_functor( MapData arg_map )
    : node_blocks( const_cast<node_block_type *>(arg_map.node_blocks.ptr_on_device()) )
    , hashes( const_cast<node_atomic *>(arg_map.hashes.ptr_on_device()) )
  {
    parallel_for( arg_map.hashes.size(), *this);
    device_type::fence();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  node_type & get_node(uint32_t i) const
  {
    return node_blocks[i>>node_block_type::shift].nodes[i&node_block_type::mask];
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i) const
  {
    uint64_t * prev_atomic = &hashes[i].value;

    while (node_atomic::next(*prev_atomic) != node_atomic::invalid_next) {
      uint64_t * curr_atomic = &get_node( node_atomic::next(*prev_atomic)).atomic.value;
      uint64_t prev = *prev_atomic;
      uint64_t curr = *curr_atomic;
      if (node_atomic::state(curr) == PENDING_DELETE) {
        //remove the node
        *prev_atomic = node_atomic::make_atomic( node_atomic::next(curr), node_atomic::state(prev) );
        *curr_atomic = node_atomic::make_atomic( node_atomic::invalid_next, UNUSED );
      }
      else {
        prev_atomic = curr_atomic;
      }
    }
  }
};

template <typename Key, typename T, typename Device, typename Compare, typename Hash>
struct map_data
{
  typedef map_data<Key,T,Device,Compare,Hash> self_type;

  typedef typename remove_const<Key>::type key_type;
  typedef typename add_const<Key>::type const_key_type;

  typedef typename remove_const<T>::type mapped_type;
  typedef typename add_const<T>::type const_mapped_type;

  typedef Device device_type;
  typedef Compare compare_type;
  typedef Hash hash_type;

  typedef map_data< key_type, mapped_type, Device, Compare, Hash>              insertable_map_type;
  typedef map_data< const_key_type, mapped_type, Device, Compare, Hash>        modifiable_map_type;
  typedef map_data< const_key_type, const_mapped_type, Device, Compare, Hash>  const_map_type;

  static const bool has_const_key_type = is_const<Key>::value;
  static const bool has_void_mapped_type = is_same<T,void>::value;
  static const bool has_const_mapped_type = has_void_mapped_type || is_const<T>::value;
  static const bool is_const_map = has_const_key_type && has_const_mapped_type;


  typedef pair<const_key_type, mapped_type> value_type;

  typedef typename if_c< is_const_map, value_type const *, value_type *>::type pointer;
  typedef value_type const * const_pointer;

  typedef node<value_type> node_type;
  typedef node_block<node_type> node_block_type;


  typedef uint32_t size_type;

  typedef typename if_c<   has_const_key_type
                         , View< const node_atomic *, device_type, MemoryTraits<RandomRead> >
                         , View< node_atomic *, device_type >
                       >::type hash_view;

  typedef typename if_c<   is_const_map
                         , View< const node_block_type *, device_type, MemoryTraits<RandomRead> >
                         , View< node_block_type *, device_type >
                       >::type node_block_view;


  typedef View< bool, device_type > bool_view;

  map_data(  uint32_t num_nodes
                     , compare_type compare
                     , hash_type hash
                    )
    : node_blocks("unordered_map_nodes", find_hash_size(static_cast<uint32_t>((num_nodes+node_block_type::size-1u)/node_block_type::size)))
    , hashes("unordered_map_hashes", find_hash_size(capacity()) )
    , no_failed_inserts("unordered_map_no_failed_inserts")
    , key_compare(compare)
    , key_hash(hash)
  {
    typedef Kokkos::DeepCopy< typename device_type::memory_space, Kokkos::HostSpace > deep_copy;
    const bool true_ = true;
    deep_copy(no_failed_inserts.ptr_on_device(), &true_, sizeof(bool));
  }

  template <typename MMapType>
  KOKKOS_INLINE_FUNCTION
  map_data( const MMapType & m)
    : node_blocks(m.node_blocks)
    , hashes(m.hashes)
    , no_failed_inserts(m.no_failed_inserts)
    , key_compare(m.key_compare)
    , key_hash(m.key_hash)
  {}

  template <typename MMapType>
  KOKKOS_INLINE_FUNCTION
  map_data & operator=( const MMapType & m)
  {
    node_blocks = m.node_blocks;
    hashes = m.hashes;
    no_failed_inserts = m.no_failed_inserts;
    key_compare = m.key_compare;
    key_hash = m.key_hash;

    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  uint32_t capacity() const
  {
    return node_blocks.size() * node_block_type::size;
  }

  node_state_counts get_node_states() const
  {
    node_state_counts result;
    count_node_states_functor<const_map_type>(*this, result);
    device_type::fence();
    return result;
  }

  hash_list_sanity_type check_hash_sanity() const
  {
    hash_list_sanity_type result;
    check_hash_list_functor<const_map_type>(*this, result);
    device_type::fence();
    return result;
  }

  void check_sanity() const
  {
    sanity_type result;

    count_node_states_functor<const_map_type>(*this, result.state_count);
    check_hash_list_functor<const_map_type>(*this, result.hash_list);

    device_type::fence();

    std::ostringstream out;

    const bool failed = failed_inserts();

    if (failed) {
      //out << "Error: " << failed << " failed insertions\n";
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
    if (result.state_count.value[ INVALID ] > 0u) {
      out << "Error: found " << result.state_count.value[ INVALID] << " invalid nodes \n";
      ++result.total_errors;
    }

    // state_count[ PENDING_INSERT ] == 0
    if (result.state_count.value[ PENDING_INSERT ] > 0u) {
      out << "Error: found " << result.state_count.value[ PENDING_INSERT] << " pending insert nodes (should always be 0)\n";
      ++result.total_errors;
    }

    if (result.total_errors > 0u) {
      out << "  Details:\n";
      out << "    Node States: ";
      out << " {UNUSED: " << result.state_count.value[ UNUSED ] << "}";
      out << " {USED: " << result.state_count.value[ USED ] << "}";
      out << " {PENDING_DELETE: " << result.state_count.value[ PENDING_DELETE ] << "}";
      out << " {PENDING_INSERT: " << result.state_count.value[ PENDING_INSERT ] << "}";
      out << " {INVALID: " << result.state_count.value[ INVALID ] << "}\n";
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
    remove_pending_delete_keys_functor<self_type> remove_keys(*this);
  }

  bool failed_inserts() const
  {
    typedef Kokkos::DeepCopy< Kokkos::HostSpace, typename device_type::memory_space > deep_copy;
    bool result = false;
    deep_copy(&result, no_failed_inserts.ptr_on_device(), sizeof(bool));
    return !result;
  }

  KOKKOS_INLINE_FUNCTION
  uint32_t find_node_index( const key_type & k) const
  {
    const uint32_t hash_value = key_hash(k);
    const uint32_t hash_index = hash_value % hashes.size();

    uint64_t prev = hashes[hash_index];

    uint32_t index = node_atomic::invalid_next;
    do {
      const uint32_t curr_index = node_atomic::next(prev);

      if ( curr_index != node_atomic::invalid_next ) {
        const node_type & curr_node = get_node(curr_index);
        const uint64_t curr = get_node(curr_index).atomic;

        const bool curr_greater = key_compare( k, curr_node.value.first);
        const bool curr_less =  key_compare( curr_node.value.first, k);
        const bool curr_equal = !curr_less && !curr_greater;

        if (curr_greater) {
          index = node_atomic::invalid_next;
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

  typedef typename if_c< is_const_map, const node_type &, node_type &>::type get_node_type;
  KOKKOS_FORCEINLINE_FUNCTION
  get_node_type get_node(uint32_t i) const
  {
    return node_blocks[i>>node_block_type::shift].nodes[i&node_block_type::mask];
  }

  // Data members
  node_block_view node_blocks;
  hash_view       hashes;
  bool_view       no_failed_inserts;
  compare_type    key_compare;
  hash_type       key_hash;
};












}}} // namespace Kokkos::Impl::UnorderedMap

#endif //KOKKOS_UNORDERED_MAP_IMPL_HPP

