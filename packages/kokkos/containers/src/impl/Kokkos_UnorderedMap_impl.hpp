#ifndef KOKKOS_UNORDERED_MAP_IMPL_HPP
#define KOKKOS_UNORDERED_MAP_IMPL_HPP

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


// type used for atomic compare and swap in the unordered map
union unordered_map_atomic
{

  enum node_state
  {
      UNUSED = 0          // not used in a list
    , USED = 1            // used in a list
    , PENDING_INSERT = 2  // not used in a list, but reserved by a thread for inserting
    , MARKED_DELETED = 3  // node in the list is marked deleted
    , INVALID = 4         // the 0th node in the node view is set to invalid
  };

  struct value_type {
    volatile uint32_t next;
    volatile int32_t state;
  };

  KOKKOS_INLINE_FUNCTION
  unordered_map_atomic(uint64_t a = 0u)
    : m_atomic(a)
  {}

  KOKKOS_INLINE_FUNCTION
  static bool compare_and_swap( volatile uint64_t * address, uint64_t o, uint64_t n)
  {
    return atomic_compare_exchange_strong(address,o,n);
  }

  KOKKOS_INLINE_FUNCTION
  uint32_t next() const
  { return m_value.next; }

  KOKKOS_INLINE_FUNCTION
  node_state state() const
  { return static_cast<node_state>(m_value.state); }

  KOKKOS_INLINE_FUNCTION
  bool unused() const
  { return m_value.state == UNUSED; }

  KOKKOS_INLINE_FUNCTION
  bool used() const
  { return m_value.state == USED; }

  KOKKOS_INLINE_FUNCTION
  bool pending_insert() const
  { return m_value.state == PENDING_INSERT; }

  KOKKOS_INLINE_FUNCTION
  bool marked_deleted() const
  { return m_value.state == MARKED_DELETED; }

  KOKKOS_INLINE_FUNCTION
  bool invalid() const
  { return m_value.state == INVALID; }

  // set the state of the node
  KOKKOS_INLINE_FUNCTION
  bool atomic_set_state( unordered_map_atomic old_value, node_state s)
  {
    unordered_map_atomic new_value;
    new_value.m_value.next = old_value.m_value.next;
    new_value.m_value.state = s;

    return compare_and_swap( &m_atomic, old_value.m_atomic, new_value.m_atomic);
  }

  // set the next node in the list
  KOKKOS_INLINE_FUNCTION
  bool atomic_set_next( unordered_map_atomic old_value, uint32_t n)
  {
    unordered_map_atomic new_value;
    new_value.m_value.next = n;
    new_value.m_value.state = old_value.m_value.state;

    return compare_and_swap( &m_atomic, old_value.m_atomic, new_value.m_atomic);
  }

  // set both state and value at the same time
  KOKKOS_INLINE_FUNCTION
  bool atomic_set( unordered_map_atomic old_value, unordered_map_atomic new_value)
  {
    return compare_and_swap( &m_atomic, old_value.m_atomic, new_value.m_atomic);
  }

  volatile uint64_t m_atomic;  // value used for compare and swap
  value_type m_value; // next and state
};


template <typename ValueType>
struct unordered_map_node
{
  typedef ValueType value_type;

  // contruct a new value at the current node
  KOKKOS_INLINE_FUNCTION
  void construct_value( const value_type & v )
  { new (&value) value_type(v); }

  // destruct the value at the current node
  KOKKOS_INLINE_FUNCTION
  void destruct_value()
  { value.~value_type(); }

  unordered_map_atomic atomic;
  value_type  value;
};




template <class MapType, bool is_const_map = MapType::is_const_map>
struct init_unordered_map
{
  typedef typename MapType::device_type    device_type;
  typedef typename MapType::node_view_type node_view_type;
  typedef typename MapType::hash_view_type hash_view_type;
  typedef typename MapType::free_count_view_type free_count_view_type;
  typedef typename device_type::size_type size_type ;

  init_unordered_map(node_view_type /*nodes*/, free_count_view_type /*free_counts*/, uint32_t /*free_block_size*/ )
  {
    printf("Error: calling init_unordered_map from a const view\n");
  }
};

template <class MapType>
struct init_unordered_map<MapType, false /*not a const map*/>
{
  typedef typename MapType::device_type    device_type;
  typedef typename MapType::node_view_type node_view_type;
  typedef typename MapType::hash_view_type hash_view_type;
  typedef typename MapType::free_count_view_type free_count_view_type;
  typedef typename device_type::size_type size_type ;

  uint32_t             m_free_block_size;
  node_view_type       m_nodes;
  free_count_view_type m_free_counts;

  init_unordered_map(node_view_type nodes, free_count_view_type free_counts, uint32_t free_block_size )
    : m_free_block_size(free_block_size)
    , m_nodes(nodes)
    , m_free_counts(free_counts)
  {
    parallel_for( m_free_counts.size(), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const
  {
    //don't count 0th node
    const uint64_t offset = i * m_free_block_size;

    const uint64_t starting_offset =   ((i > 0u) && (offset < m_nodes.size()))
                                     ? offset
                                     : ((i == 0u) ? 1u : m_nodes.size());
    const uint64_t ending_offset =   (((i+1) * m_free_block_size ) <= m_nodes.size())
                                   ? ((i+1) * m_free_block_size )
                                   : m_nodes.size();

    m_free_counts[i] = static_cast<uint32_t>(ending_offset - starting_offset);

    if (i==0u) {
      m_nodes[0].atomic.m_value.state = Impl::unordered_map_atomic::INVALID;
    }
  }
};



template <class MapType>
struct unordered_map_size_functor
{
  typedef typename MapType::device_type    device_type;
  typedef typename MapType::node_view_type node_view_type;
  typedef typename MapType::hash_view_type hash_view_type;
  typedef typename MapType::free_count_view_type free_count_view_type;
  typedef typename device_type::size_type size_type ;
  typedef uint32_t value_type;

  node_view_type m_nodes;

  unordered_map_size_functor(node_view_type nodes, value_type & value)
    : m_nodes(nodes)
  {
    parallel_reduce( m_nodes.size(), *this, value);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & dst)
  {
    dst = 0;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & dst, const volatile value_type & src)
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & dst) const
  {
    if ( m_nodes[i].atomic.used() ) {
      ++dst;
    }
  }
};



template <class MapType>
struct unordered_map_count_unused_functor
{
  typedef typename MapType::device_type    device_type;
  typedef typename MapType::node_view_type node_view_type;
  typedef typename MapType::hash_view_type hash_view_type;
  typedef typename MapType::free_count_view_type free_count_view_type;
  typedef typename device_type::size_type size_type ;
  typedef uint32_t value_type;

  free_count_view_type m_free_counts;

  unordered_map_count_unused_functor(free_count_view_type free_counts, value_type & value)
    : m_free_counts(free_counts)
  {
    parallel_reduce( m_free_counts.size(), *this, value);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & dst)
  {
    dst = 0;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & dst, const volatile value_type & src)
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & dst) const
  {
    dst += m_free_counts[i];
  }
};




template <class MapType, bool is_const_map = MapType::is_const_map>
struct unordered_map_delete_marked_keys_functor
{
  typedef typename MapType::device_type    device_type;
  typedef typename MapType::node_view_type node_view_type;
  typedef typename MapType::hash_view_type hash_view_type;
  typedef typename MapType::free_count_view_type free_count_view_type;
  typedef typename device_type::size_type size_type ;
  typedef uint32_t value_type;

  unordered_map_delete_marked_keys_functor(  node_view_type //nodes
      , free_count_view_type //free_counts
      , hash_view_type //hashes
      , uint32_t //free_block_size
      , value_type & //num_deleted
      )
  {
    printf("Error: calling unordered_map_delete_marked_keys from a const view\n");
  }
};


template <class MapType>
struct unordered_map_delete_marked_keys_functor<MapType, false /*not a const map*/>
{
  typedef typename MapType::device_type    device_type;
  typedef typename MapType::node_view_type node_view_type;
  typedef typename MapType::hash_view_type hash_view_type;
  typedef typename MapType::free_count_view_type free_count_view_type;
  typedef typename device_type::size_type size_type ;
  typedef uint32_t value_type;
  uint32_t             m_free_block_size;
  node_view_type       m_nodes;
  free_count_view_type m_free_counts;
  hash_view_type       m_hashes;

  unordered_map_delete_marked_keys_functor(  node_view_type nodes
      , free_count_view_type free_counts
      , hash_view_type hashes
      , uint32_t free_block_size
      , value_type & num_deleted
      )
    : m_free_block_size(free_block_size)
      , m_nodes(nodes)
      , m_free_counts(free_counts)
      , m_hashes(hashes)
  {
    parallel_reduce( m_hashes.size(), *this, num_deleted);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & dst)
  {
    dst = 0;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & dst, const volatile value_type & src)
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & dst) const
  {
    Impl::unordered_map_atomic * prev_atomic = &m_hashes[i];
    Impl::unordered_map_atomic * curr_atomic = &m_nodes[prev_atomic->next()].atomic;

    while (!curr_atomic->invalid()) {
      if(curr_atomic->marked_deleted()) {

        prev_atomic->m_value.next = curr_atomic->m_value.next;
        curr_atomic->m_value.state = Impl::unordered_map_atomic::UNUSED;

        uint32_t free_block = prev_atomic->next() / m_free_block_size;

        volatile uint32_t * addr = &m_free_counts[free_block];
        atomic_fetch_add( addr, 1u);
        ++dst;
      }
      prev_atomic = curr_atomic;
      curr_atomic = &m_nodes[prev_atomic->next()].atomic;
    }
  }
};

template <class MapType>
struct unordered_map_count_duplicate_keys_functor
{
  typedef typename MapType::device_type    device_type;
  typedef typename MapType::node_view_type node_view_type;
  typedef typename MapType::hash_view_type hash_view_type;
  typedef typename MapType::free_count_view_type free_count_view_type;
  typedef typename device_type::size_type size_type ;
  typedef uint32_t value_type;

  node_view_type       m_nodes;
  hash_view_type       m_hashes;

  unordered_map_count_duplicate_keys_functor(  node_view_type nodes
      , hash_view_type hashes
      , value_type & num_deleted
      )
    : m_nodes(nodes)
      , m_hashes(hashes)
  {
    parallel_reduce( m_hashes.size(), *this, num_deleted);
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & dst)
  {
    dst = 0;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & dst, const volatile value_type & src)
  { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & errors) const
  {
    Impl::unordered_map_atomic * prev_atomic = &m_hashes[i];
    Impl::unordered_map_atomic * curr_atomic = &m_nodes[prev_atomic->next()].atomic;
    Impl::unordered_map_atomic * next_atomic = &m_nodes[curr_atomic->next()].atomic;

    while (!curr_atomic->invalid() && !next_atomic->invalid()) {
      if (m_nodes[prev_atomic->next()].value.first == m_nodes[curr_atomic->next()].value.first)
      {
        ++errors;
        printf("Duplicate key %d\n", m_nodes[prev_atomic->next()].value.first);
      }

      prev_atomic = curr_atomic;
      curr_atomic = next_atomic;
      next_atomic = &m_nodes[curr_atomic->next()].atomic;
    }
  }
};



}} // namespace Kokkos::Impl

#endif //KOKKOS_UNORDERED_MAP_IMPL_HPP

