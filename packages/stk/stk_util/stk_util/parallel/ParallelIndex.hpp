#ifndef stk_util_parallel_ParallelIndex_hpp
#define stk_util_parallel_ParallelIndex_hpp

#include <utility>
#include <vector>
#include <algorithm>

#include <stdint.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/util/PairIter.hpp>

namespace stk {
namespace util {

template <class K, class P>
struct ParallelIndexDecomp
{
  typedef K Key ;
  typedef P Proc ;

  Proc operator()( const Proc proc_size, const Key key ) const {
    return ( key >> 8 ) % proc_size ;
  };
};


/** \brief Parallel cross-reference index for a collection of
 *  'Key' keys.
 *
 *  Each processor constructs a ParallelIndex with its
 *  local collection of keys.
 *  The resulting ParallelIndex may be queried for
 *  - which other processors input the same local keys or
 *  - which processors submitted an arbitrary set of keys.
 */
template <class K = uint64_t, class P = unsigned, class D = ParallelIndexDecomp<K, P> >
class ParallelIndex {
public:
  typedef K Key ;
  typedef P Proc ;
  typedef D Decomp ;
  typedef std::pair<Key, Proc> KeyProc ;

#ifndef DOXYGEN_COMPILE
  struct LessKeyProc {
    bool operator()( const KeyProc & lhs, const KeyProc & rhs ) const {
      return lhs < rhs ;
    }

    bool operator()( const KeyProc & lhs, const Key rhs ) const {
      return lhs.first < rhs ;
    }
  };

  struct EqualKeyProc {
    bool operator()( const KeyProc & lhs, const KeyProc & rhs ) const {
      return lhs == rhs ;
    }
  };
#endif /* DOXYGEN_COMPILE */


  /** \brief  Construct with locally-submitted keys */
  ParallelIndex( ParallelMachine comm, 
                 const std::vector<Key> & local )
    : m_comm( comm ), 
      m_key_proc(), 
      m_decomp()
  {
    const unsigned p_size = parallel_machine_size( comm );

    CommAll all( comm );

    pack_map( all, local );

    all.allocate_buffers( p_size / 4, false );
  
    pack_map( all, local );

    all.communicate();

    unpack_map( all, m_key_proc );

    sort_unique( m_key_proc );
  }

  ~ParallelIndex()
  {}
  

  /** \brief  Query which other processors submitted the
   *          same keys that the local processor submitted.
   */
  void query(std::vector<KeyProc> & global ) const
  {
    const unsigned p_size = parallel_machine_size( m_comm );

    CommAll all( m_comm );

    pack_query( all, m_key_proc );

    all.allocate_buffers( p_size / 4, false );

    pack_query( all, m_key_proc );

    all.communicate();

    unpack_query( all, global );

    sort_unique( global );
  }


  /** \brief  Query which processors submitted the given keys.
   *          The local processor is in the output if it
   *          submitted a queried key.
   */
  void query(const std::vector<Key> & local, 
             std::vector<KeyProc> & global ) const
  {
    const unsigned p_size = parallel_machine_size( m_comm );

    std::vector<KeyProc> tmp ;

    {
      CommAll all( m_comm );

      pack_map( all, local );

      all.allocate_buffers( p_size / 4, false );

      pack_map( all, local );

      all.communicate();

      unpack_map( all, tmp ); // { ( key, querying_processor ) }

      sort_unique( tmp );
    }

    {
      CommAll all( m_comm );

      pack_query( all, m_key_proc, tmp );

      all.allocate_buffers( p_size / 4, false );

      pack_query( all, m_key_proc, tmp );

      all.communicate();

      unpack_query( all, global );

      sort_unique( global );
    }
  }

private:
  void sort_unique( std::vector<KeyProc> & key_proc ) const
  {
    typename std::vector<KeyProc>::iterator i = key_proc.begin();
    typename  std::vector<KeyProc>::iterator j = key_proc.end();

    std::sort( i, j, LessKeyProc() );
    i = std::unique( i, j, EqualKeyProc() );
    key_proc.erase( i, j );
  }

  void pack_map( CommAll & all, const std::vector<Key> & local ) const
  {
    const unsigned p_size = all.parallel_size();

    typename std::vector<Key>::const_iterator i ;

    for ( i = local.begin() ; i != local.end() ; ++i ) {
      const Key value = *i ;
      const unsigned proc = m_decomp( p_size, value );
      CommBuffer & buf = all.send_buffer(proc);
      buf.pack<Key>( value );
    }
  }

  void unpack_map( CommAll & all, std::vector< KeyProc > & key_proc ) const
  {
    const unsigned p_size = all.parallel_size();

    unsigned count = 0 ;
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      count += all.recv_buffer( p ).capacity() / sizeof(Key);
    }

    key_proc.clear();
    key_proc.reserve( count );

    KeyProc value ;

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = all.recv_buffer( p );
      value.second = p ;
      while ( buf.remaining() ) {
        buf.unpack<Key>( value.first );
        key_proc.push_back( value );
      }
    }
  }

  void pack_query( CommAll & all, const std::vector< KeyProc > & key_proc ) const
  {
    KeyProc value ;

    typename std::vector< KeyProc >::const_iterator i ;

    for ( i = key_proc.begin() ; i != key_proc.end() ; ) {
      value.first = i->first ;

      const typename std::vector< KeyProc >::const_iterator i_beg = i ;

      for ( ; i != key_proc.end() && value.first == i->first ; ++i )
        ;

      const typename std::vector< KeyProc >::const_iterator i_end = i ;
    
      for ( i = i_beg ; i != i_end ; ++i ) {
        CommBuffer & buf = all.send_buffer( i->second );

        typename std::vector< KeyProc >::const_iterator j ;

        for ( j = i_beg ; j != i_end ; ++j ) {
          if ( j != i ) {
            value.second = j->second ;
            buf.pack<KeyProc>( value );
          }
        }
      }
    }
  }

  void pack_query( CommAll & all, 
                   const std::vector< KeyProc > & key_proc_map, 
                   const std::vector< KeyProc > & query ) const
  {
    KeyProc value ;

    for (typename std::vector< KeyProc >::const_iterator i = query.begin() ; i != query.end() ; ) {
      value.first = i->first ;

      typename std::vector< KeyProc >::const_iterator key_begin = std::lower_bound( key_proc_map.begin(), key_proc_map.end(), value.first, LessKeyProc() );

      typename std::vector< KeyProc >::const_iterator key_end = key_begin;
      while ( key_end != key_proc_map.end() && key_end->first == value.first)
        ++key_end;

      for ( ; i != query.end() && value.first == i->first ; ++i ) {
        CommBuffer & buf = all.send_buffer( i->second );

        for ( typename std::vector< KeyProc >::const_iterator j = key_begin ; j != key_end ; ++j ) {
          value.second = j->second ;
          buf.pack<KeyProc>( value );
        }
      }
    }
  }

  void unpack_query( CommAll & all, std::vector< KeyProc > & key_proc ) const
  {
    const unsigned p_size = all.parallel_size();

    KeyProc entry ;

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = all.recv_buffer( p );
      while ( buf.remaining() ) {
        buf.unpack<KeyProc>( entry );
        key_proc.push_back( entry );
      }
    }
  }

private:
  ParallelIndex();
  ParallelIndex( const ParallelIndex & );
  ParallelIndex & operator = ( const ParallelIndex & );

private:
  ParallelMachine               m_comm ;
  std::vector<KeyProc>          m_key_proc ;
  Decomp                        m_decomp;
};

} // namespace util
} // namespace stk

//----------------------------------------------------------------------

#endif

