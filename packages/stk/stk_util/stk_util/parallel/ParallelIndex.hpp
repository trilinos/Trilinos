// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_util_parallel_ParallelIndex_hpp
#define stk_util_parallel_ParallelIndex_hpp

#include <utility>
#include <vector>
#include <algorithm>

#include <stdint.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>
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
    CommSparse sparse( comm );

    pack_map( sparse, local );

    sparse.allocate_buffers();
  
    pack_map( sparse, local );

    sparse.communicate();

    unpack_map( sparse, m_key_proc );

    sort_unique( m_key_proc );
  }

  ~ParallelIndex()
  {}
  

  /** \brief  Query which other processors submitted the
   *          same keys that the local processor submitted.
   */
  void query(std::vector<KeyProc> & global ) const
  {
    CommSparse sparse( m_comm );

    pack_query( sparse, m_key_proc );

    sparse.allocate_buffers();

    pack_query( sparse, m_key_proc );

    sparse.communicate();

    unpack_query( sparse, global );

    sort_unique( global );
  }


  /** \brief  Query which processors submitted the given keys.
   *          The local processor is in the output if it
   *          submitted a queried key.
   */
  void query(const std::vector<Key> & local, 
             std::vector<KeyProc> & global ) const
  {
    std::vector<KeyProc> tmp ;

    {
      CommSparse sparse( m_comm );

      pack_map( sparse, local );

      sparse.allocate_buffers();

      pack_map( sparse, local );

      sparse.communicate();

      unpack_map( sparse, tmp ); // { ( key, querying_processor ) }

      sort_unique( tmp );
    }

    {
      CommSparse sparse( m_comm );

      pack_query( sparse, m_key_proc, tmp );

      sparse.allocate_buffers();

      pack_query( sparse, m_key_proc, tmp );

      sparse.communicate();

      unpack_query( sparse, global );

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

  void pack_map( CommSparse & sparse, const std::vector<Key> & local ) const
  {
    const unsigned p_size = sparse.parallel_size();

    typename std::vector<Key>::const_iterator i ;

    for ( i = local.begin() ; i != local.end() ; ++i ) {
      const Key value = *i ;
      const unsigned proc = m_decomp( p_size, value );
      CommBuffer & buf = sparse.send_buffer(proc);
      buf.pack<Key>( value );
    }
  }

  void unpack_map( CommSparse & sparse, std::vector< KeyProc > & key_proc ) const
  {
    const unsigned p_size = sparse.parallel_size();

    unsigned count = 0 ;
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      count += sparse.recv_buffer( p ).capacity() / sizeof(Key);
    }

    key_proc.clear();
    key_proc.reserve( count );

    KeyProc value ;

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = sparse.recv_buffer( p );
      value.second = p ;
      while ( buf.remaining() ) {
        buf.unpack<Key>( value.first );
        key_proc.push_back( value );
      }
    }
  }

  void pack_query( CommSparse & sparse, const std::vector< KeyProc > & key_proc ) const
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
        CommBuffer & buf = sparse.send_buffer( i->second );

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

  void pack_query( CommSparse & sparse, 
                   const std::vector< KeyProc > & key_proc_map, 
                   const std::vector< KeyProc > & query_in ) const
  {
    KeyProc value ;

    for (typename std::vector< KeyProc >::const_iterator i = query_in.begin() ; i != query_in.end() ; ) {
      value.first = i->first ;

      typename std::vector< KeyProc >::const_iterator key_begin = std::lower_bound( key_proc_map.begin(), key_proc_map.end(), value.first, LessKeyProc() );

      typename std::vector< KeyProc >::const_iterator key_end = key_begin;
      while ( key_end != key_proc_map.end() && key_end->first == value.first)
        ++key_end;

      for ( ; i != query_in.end() && value.first == i->first ; ++i ) {
        CommBuffer & buf = sparse.send_buffer( i->second );

        for ( typename std::vector< KeyProc >::const_iterator j = key_begin ; j != key_end ; ++j ) {
          value.second = j->second ;
          buf.pack<KeyProc>( value );
        }
      }
    }
  }

  void unpack_query( CommSparse & sparse, std::vector< KeyProc > & key_proc ) const
  {
    const unsigned p_size = sparse.parallel_size();

    KeyProc entry ;

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = sparse.recv_buffer( p );
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

