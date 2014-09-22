/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_util/parallel/DistributedIndex.hpp>
#include <algorithm>                    // for sort, lower_bound, min, etc
#include <iterator>                     // for distance
#include <limits>                       // for numeric_limits
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll, CommBuffer
#include <stk_util/util/RadixSort.hpp>  // for radix_sort_unsigned
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_rank, etc
#include "stk_util/stk_config.h"        // for STK_HAS_MPI


namespace stk {
namespace parallel {

//----------------------------------------------------------------------

namespace {

struct KeyProcLess {

  bool operator()( const DistributedIndex::KeyProc & lhs ,
                   const DistributedIndex::KeyType & rhs ) const
  { return lhs.first < rhs ; }

};

template <typename Vector>
void sort_unique( Vector & vec )
{
  typename Vector::iterator  i = vec.begin() , j = vec.end() ;

  std::sort( i , j );
  i = std::unique( i , j );
  vec.erase( i , j );
}


// reserve vector size (current size + rev_buffer remaining)
template < typename Vector >
inline void reserve_for_recv_buffer( const CommAll& all, const DistributedIndex::ProcType& comm_size, Vector & v)
{
  unsigned num_remote = 0;
  for (DistributedIndex::ProcType p = 0 ; p < comm_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    num_remote += buf.remaining() / sizeof(typename Vector::value_type);
  }
  v.reserve(v.size() + num_remote);
}

// unpack buffer into vector
template < typename Vector >
inline void unpack_recv_buffer( const CommAll& all, const DistributedIndex::ProcType& comm_size, Vector & v)
{
  typedef typename Vector::value_type value_type;
  reserve_for_recv_buffer(all, comm_size, v);
  for (DistributedIndex::ProcType p = 0 ; p < comm_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    while ( buf.remaining() ) {
      value_type kp;
      buf.unpack( kp );
      v.push_back( kp );
    }
  }
}

// unpack buffer into vector, where pair.second is the processor
template < typename VectorProcPair >
inline void unpack_with_proc_recv_buffer( const CommAll& all, const DistributedIndex::ProcType& comm_size, VectorProcPair & v)
{
  typedef typename VectorProcPair::value_type pair_type;
  reserve_for_recv_buffer(all, comm_size, v);
  for ( DistributedIndex::ProcType p = 0 ; p < comm_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    pair_type kp;
    kp.second = p;
    while ( buf.remaining() ) {
      buf.unpack( kp.first );
      v.push_back( kp );
    }
  }
}

} // namespace <unnamed>

//----------------------------------------------------------------------

enum { DISTRIBUTED_INDEX_CHUNK_BITS = 8 }; ///< Each chunk is 256 keys

enum { DISTRIBUTED_INDEX_CHUNK_SIZE =
       size_t(1) << DISTRIBUTED_INDEX_CHUNK_BITS };

DistributedIndex::ProcType
DistributedIndex::to_which_proc( const DistributedIndex::KeyType & key ) const
{
  for ( size_t i = 0 ; i < m_span_count ; ++i ) {
    if ( m_key_span[i].first <= key && key <= m_key_span[i].second ) {
       const KeyType offset = key - m_key_span[i].first ;
       return ( offset >> DISTRIBUTED_INDEX_CHUNK_BITS ) % m_comm_size ;
    }
  }
  return m_comm_size ;
}

//----------------------------------------------------------------------

DistributedIndex::~DistributedIndex() {}

DistributedIndex::DistributedIndex (
  ParallelMachine comm ,
  const KeySpanVector & partition_bounds )
  : m_comm( comm ),
    m_comm_rank( parallel_machine_rank( comm ) ),
    m_comm_size( parallel_machine_size( comm ) ),
    m_span_count(0),
    m_key_span(),
    m_removed_keys(),
    m_key_usage()
{
  unsigned info[2] ;
  info[0] = partition_bounds.size();
  info[1] = 0 ;

  // Check each span for validity

  for ( KeySpanVector::const_iterator
        i = partition_bounds.begin() ; i != partition_bounds.end() ; ++i ) {
    if ( i->second < i->first ||
         ( i != partition_bounds.begin() && i->first <= (i-1)->second ) ) {
      info[1] = 1 ;
    }
  }

#if defined( STK_HAS_MPI )
  if (m_comm_size > 1) {
    // BABBLE_STK_PARALLEL_COMM(comm, "    calling MPI_Bcast from DistributedIndex");
    MPI_Bcast( info , 2 , MPI_UNSIGNED , 0 , comm );
  }

  if ( 0 < info[0] ) {
    m_key_span.resize( info[0] );
    if ( 0 == parallel_machine_rank( comm ) ) {
      m_key_span.assign(partition_bounds.begin(),partition_bounds.end());
    }
    if (m_comm_size > 1) {
      // BABBLE_STK_PARALLEL_COMM(comm, "    calling MPI_Bcast from DistributedIndex");
      MPI_Bcast( (m_key_span.empty() ? NULL : & m_key_span[0]), info[0] * sizeof(KeySpan), MPI_BYTE, 0, comm );
    }
  }
#else
 m_key_span.assign(partition_bounds.begin(),partition_bounds.end());
#endif

  if ( info[1] ) {
    std::ostringstream msg ;
    msg << "sierra::parallel::DistributedIndex ctor( comm , " ;

    for ( KeySpanVector::const_iterator
          i = partition_bounds.begin() ; i != partition_bounds.end() ; ++i ) {
      msg << " ( min = " << i->first << " , max = " << i->second << " )" ;
    }
    msg << " ) contains invalid span of keys" ;
    throw std::runtime_error( msg.str() );
  }

  m_span_count = info[0] ;

  if ( 0 == m_span_count ) {
    m_key_span.push_back(
      KeySpan( std::numeric_limits<KeyType>::min(),
               std::numeric_limits<KeyType>::max() ) );
    m_span_count = 1 ;
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template <typename Vector>
bool is_sorted_and_unique( const Vector & v )
{
  for (size_t i=1u, size = v.size(); i<size; ++i) {
    if ( v[i-1] >= v[i] ) return false;
  }
  return true;
}

template <typename KeyUsageVector, typename KeyRequestVector>
void query_pack_to_usage(
  const KeyUsageVector & key_usage ,
  const KeyRequestVector & request ,
  CommAll & all )
{
  typedef typename KeyUsageVector::const_iterator   usage_iterator_type;
  typedef typename KeyRequestVector::const_iterator request_iterator_type;

  usage_iterator_type i = key_usage.begin();
  request_iterator_type k = request.begin();

  for ( ; k != request.end() && i != key_usage.end() ; ++k ) {

    for ( ; i != key_usage.end() && i->first < *k ; ++i ) {}

    usage_iterator_type j = i ;
    for ( ; j != key_usage.end() && j->first == *k ; ++j ) {}

    for ( usage_iterator_type jsend = i ; jsend != j ; ++jsend ) {

      for ( usage_iterator_type jinfo = i ; jinfo != j ; ++jinfo ) {
        DistributedIndex::ProcType proc = jsend->second;
        DistributedIndex::KeyProc key_proc(*jinfo);
        all.send_buffer( proc ).pack<DistributedIndex::KeyProc>( key_proc );
      }
    }
  }
}

template <typename KeyUsageVector, typename KeyRequestVector>
void query_pack( const KeyUsageVector & key_usage ,
                 const KeyRequestVector & request ,
                 CommAll & all )
{
  typedef typename KeyUsageVector::const_iterator   usage_iterator_type;
  typedef typename KeyRequestVector::const_iterator request_iterator_type;

  usage_iterator_type i = key_usage.begin();

  for ( request_iterator_type k = request.begin(); k != request.end() && i != key_usage.end() ; ++k ) {
    for ( ; i != key_usage.end() && i->first < k->first ; ++i ) {}
    for ( usage_iterator_type j = i; j != key_usage.end() && j->first == k->first ; ++j ) {
      DistributedIndex::ProcType proc = k->second;
      DistributedIndex::KeyProc key_proc(*j);
      all.send_buffer( proc ).pack<DistributedIndex::KeyProc>( key_proc );
    }
  }
}

} // unnamed namespace

void DistributedIndex::query( const KeyProcVector & request , KeyProcVector & sharing_of_keys ) const
{
  sharing_of_keys.clear();

  CommAll all( m_comm );

  query_pack( m_key_usage , request , all ); // Sizing

  BABBLE_STK_PARALLEL_COMM(m_comm, "          query calling allocate_buffers")
  all.allocate_buffers( m_comm_size / 4 , false );

  query_pack( m_key_usage , request , all ); // Packing

  BABBLE_STK_PARALLEL_COMM(m_comm, "          query calling communicate")
  all.communicate();

  unpack_recv_buffer(all, m_comm_size, sharing_of_keys);

  std::sort( sharing_of_keys.begin() , sharing_of_keys.end() );
}

void DistributedIndex::query( KeyProcVector & sharing_of_local_keys ) const
{
  query( m_key_usage , sharing_of_local_keys );
}

void DistributedIndex::query( const KeyTypeVector & keys , KeyProcVector & sharing_keys ) const
{
  KeyProcVector request ;

  {
    bool bad_key = false ;
    CommAll all( m_comm );

    for ( KeyTypeVector::const_iterator
          k = keys.begin() ; k != keys.end() ; ++k ) {
      const ProcType p = to_which_proc( *k );

      if ( p < m_comm_size ) {
        all.send_buffer( p ).pack<KeyType>( *k );
      }
      else {
        bad_key = true ;
      }
    }

    // Error condition becomes global:
    BABBLE_STK_PARALLEL_COMM(m_comm, "          query calling allocate_buffers")
    bad_key = all.allocate_buffers( m_comm_size / 4 , false , bad_key );

    if ( bad_key ) {
      throw std::runtime_error("stk::parallel::DistributedIndex::query given a key which is out of range");
    }

    for ( KeyTypeVector::const_iterator
          k = keys.begin() ; k != keys.end() ; ++k ) {
      all.send_buffer( to_which_proc( *k ) ).pack<KeyType>( *k );
    }

    BABBLE_STK_PARALLEL_COMM(m_comm, "          query calling communicate")
    all.communicate();

    unpack_with_proc_recv_buffer(all, m_comm_size, request);
  }

  sort_unique( request );

  query( request , sharing_keys );
}

void DistributedIndex::query_to_usage( const KeyTypeVector & keys ,  KeyProcVector & sharing_keys ) const
{
  BABBLE_STK_PARALLEL_COMM(m_comm, "          enter query_to_usage");
  KeyTypeVector request ;

  {
    bool bad_key = false ;
    CommAll all( m_comm );

    for ( KeyTypeVector::const_iterator
          k = keys.begin() ; k != keys.end() ; ++k ) {
      const ProcType p = to_which_proc( *k );

      if ( p < m_comm_size ) {
        all.send_buffer( p ).pack<KeyType>( *k );
      }
      else {
        bad_key = true ;
      }
    }

    // Error condition becomes global:

    BABBLE_STK_PARALLEL_COMM(m_comm, "          query_to_usage calling allocate_buffers")
    bad_key = all.allocate_buffers( m_comm_size / 4 , false , bad_key );

    if ( bad_key ) {
      throw std::runtime_error("stk::parallel::DistributedIndex::query given a key which is out of range");
    }

    for ( KeyTypeVector::const_iterator
          k = keys.begin() ; k != keys.end() ; ++k ) {
      all.send_buffer( to_which_proc( *k ) ).pack<KeyType>( *k );
    }

    BABBLE_STK_PARALLEL_COMM(m_comm, "          query_to_usage calling communicate")
    all.communicate();

    unpack_recv_buffer(all, m_comm_size, request);
  }

  sort_unique( request );

  {
    CommAll all( m_comm );

    query_pack_to_usage( m_key_usage , request , all ); // Sizing

    BABBLE_STK_PARALLEL_COMM(m_comm, "          query_to_usage calling allocate_buffers")
    all.allocate_buffers( m_comm_size / 4 , false );

    query_pack_to_usage( m_key_usage , request , all ); // Packing

    BABBLE_STK_PARALLEL_COMM(m_comm, "          query_to_usage calling communicate")
    all.communicate();

    unpack_recv_buffer(all, m_comm_size, sharing_keys);

    std::sort( sharing_keys.begin() , sharing_keys.end() );
  }
  // BABBLE_STK_PARALLEL_COMM(m_comm, "        exit query_to_usage");
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

struct RemoveKeyProc {

  bool operator()( const DistributedIndex::KeyProc & kp ) const
  { return kp.second < 0 ; }

  static void mark( DistributedIndex::KeyProcVector & key_usage ,
                    const DistributedIndex::KeyProc & kp )
  {
    DistributedIndex::KeyProcVector::iterator
      i = std::lower_bound( key_usage.begin(),
                            key_usage.end(), kp.first, KeyProcLess() );

    // Iterate over the span of KeyProcs with matching key until an exact match
    // is found. We have to do it this way because marking a KeyProc unsorts it
    // in the key_usage vector, so we cannot look up KeyProcs directly once marking
    // has begun.
    while ( i != key_usage.end() && kp.first == i->first && kp.second != i->second) { ++i ; }

    if ( i != key_usage.end() && kp == *i ) {
      i->second = -1 ;
    }
  }

  static void clean( DistributedIndex::KeyProcVector & key_usage )
  {
    DistributedIndex::KeyProcVector::iterator end =
      std::remove_if( key_usage.begin() , key_usage.end() , RemoveKeyProc() );
    key_usage.erase( end , key_usage.end() );
  }
};

}

void DistributedIndex::update_keys( KeyTypeVector::const_iterator add_new_keys_begin,
                                    KeyTypeVector::const_iterator add_new_keys_end )
{
  KeyTypeVector removed_keys;
  removed_keys.swap(m_removed_keys);

  sort_unique(removed_keys);

//  profile_memory_usage<DistributedIndex>("before update_keys", m_comm, m_comm_rank);

  update_keys(add_new_keys_begin, add_new_keys_end, removed_keys.begin(), removed_keys.end());

//  profile_memory_usage<DistributedIndex>("after update_keys", m_comm, m_comm_rank);

}

void DistributedIndex::update_keys(
  KeyTypeVector::const_iterator add_new_keys_begin ,
  KeyTypeVector::const_iterator add_new_keys_end ,
  KeyTypeVector::const_iterator remove_existing_keys_begin,
  KeyTypeVector::const_iterator remove_existing_keys_end )
{
  if (!m_removed_keys.empty()) {
    throw std::runtime_error("");
  }

  UnsignedVector count_remove( m_comm_size , 0 );
  UnsignedVector count_add(    m_comm_size , 0 );

  size_t local_bad_input = 0 ;

  // Iterate over keys being removed and keep a count of keys being removed
  // from other processes
  for ( KeyTypeVector::const_iterator
        i = remove_existing_keys_begin;
        i != remove_existing_keys_end; ++i ) {
    const ProcType p = to_which_proc( *i );
    if ( m_comm_size <= p ) {
      // Key is not within one of the span:
      ++local_bad_input ;
    }
    else if ( p != m_comm_rank ) {
      ++( count_remove[ p ] );
    }
  }

  // Iterate over keys being added and keep a count of keys being added
  // to other processes
  for ( KeyTypeVector::const_iterator
        i = add_new_keys_begin;
        i != add_new_keys_end; ++i ) {
    const ProcType p = to_which_proc( *i );
    if ( p == m_comm_size ) {
      // Key is not within one of the span:
      ++local_bad_input ;
    }
    else if ( p != m_comm_rank ) {
      ++( count_add[ p ] );
    }
  }

  CommAll all( m_comm );

  // Sizing and add_new_keys bounds checking:

  // For each process, we are going to send the number of removed keys,
  // the removed keys, and the added keys. It will be assumed that any keys
  // beyond the number of removed keys will be added keys.
//  unsigned max_add = 0;
  for ( int p = 0 ; p < m_comm_size ; ++p ) {
    if ( count_remove[p] || count_add[p] ) {
      CommBuffer & buf = all.send_buffer( p );
      buf.skip<unsigned long>( 1 );
      buf.skip<KeyType>( count_remove[p] );
      buf.skip<KeyType>( count_add[p] );
//      if (count_add[p] > max_add) max_add = count_add[p];
    }
  }

  // Allocate buffers and perform a global OR of error_flag
  const bool symmetry_flag = false ;
  const bool error_flag = 0 < local_bad_input ;

  BABBLE_STK_PARALLEL_COMM(m_comm, "          update keys calling allocate_buffers");
  bool global_bad_input =
    all.allocate_buffers( m_comm_size / 4, symmetry_flag , error_flag );

  if ( global_bad_input ) {
    std::ostringstream msg ;

    if ( 0 < local_bad_input ) {
      msg << "stk::parallel::DistributedIndex::update_keys ERROR Given "
          << local_bad_input << " of " << (add_new_keys_end - add_new_keys_begin)
          << " add_new_keys outside of any span" ;
    }

    throw std::runtime_error( msg.str() );
  }

  // Packing:

  // Pack the remove counts for each process
  for ( int p = 0 ; p < m_comm_size ; ++p ) {
    if ( count_remove[p] || count_add[p] ) {
      all.send_buffer( p ).pack<unsigned long>( count_remove[p] );
    }
  }

  // Pack the removed keys for each process
  for ( KeyTypeVector::const_iterator
        i = remove_existing_keys_begin;
        i != remove_existing_keys_end; ++i ) {
    const ProcType p = to_which_proc( *i );
    if ( p != m_comm_rank ) {
      all.send_buffer( p ).pack<KeyType>( *i );
    }
  }

  // Pack the added keys for each process
  for ( KeyTypeVector::const_iterator
        i = add_new_keys_begin;
        i != add_new_keys_end; ++i ) {
    const ProcType p = to_which_proc( *i );
    if ( p != m_comm_rank ) {
      all.send_buffer( p ).pack<KeyType>( *i );
    }
  }

  BABBLE_STK_PARALLEL_COMM(m_comm, "          update keys calling communicate");
  // Communicate keys
  all.communicate();

  //------------------------------
  // Remove for local keys

  for ( KeyTypeVector::const_iterator
        i = remove_existing_keys_begin;
        i != remove_existing_keys_end; ++i ) {
    const ProcType p = to_which_proc( *i );
    if ( p == m_comm_rank ) {
      RemoveKeyProc::mark( m_key_usage , KeyProc( *i , p ) );
    }
  }

  // Unpack the remove key and find it.
  // Set the process to a negative value for subsequent removal.

  for ( int p = 0 ; p < m_comm_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    if ( buf.remaining() ) {
      unsigned long remove_count = 0 ;

      KeyProc kp ;

      kp.second = p ;

      buf.unpack<unsigned long>( remove_count );

      for ( ; 0 < remove_count ; --remove_count ) {
        buf.unpack<KeyType>( kp.first );

        RemoveKeyProc::mark( m_key_usage , kp );
      }
    }
  }

  RemoveKeyProc::clean( m_key_usage );

  //------------------------------
  // Append for local keys

  // Add new_keys going to this proc to local_key_usage
  KeyProcVector new_key_usage ;
  new_key_usage.reserve(add_new_keys_end - add_new_keys_begin);
  for ( KeyTypeVector::const_iterator
        i = add_new_keys_begin;
        i != add_new_keys_end; ++i ) {

    const ProcType p = to_which_proc( *i );
    if ( p == m_comm_rank ) {
      new_key_usage.push_back( KeyProc( *i , p ) );
    }
  }

  // Unpack and append for remote keys:

  unpack_with_proc_recv_buffer(all, m_comm_size, new_key_usage);

  m_key_usage.insert(m_key_usage.end(), new_key_usage.begin(), new_key_usage.end());
  std::sort(m_key_usage.begin(), m_key_usage.end());

  // Unique m_key_usage
  m_key_usage.erase(std::unique( m_key_usage.begin(), m_key_usage.end()),
                    m_key_usage.end() );

  // Check invariant that m_key_usage is sorted
#ifndef NDEBUG
  if (!is_sorted_and_unique(m_key_usage)) {
    throw std::runtime_error("Sorted&unique invariant violated!");
  }
#endif
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void DistributedIndex::generate_new_global_key_upper_bound(
  const std::vector<size_t>  & requests ,
        KeyTypeVector & global_key_upper_bound ) const
{
  bool bad_request = m_span_count != requests.size();

  std::ostringstream error_msg ;

  error_msg
  << "sierra::parallel::DistributedIndex::generate_new_keys_global_counts( " ;


  UnsignedVector
    local_counts( m_span_count + 1 , 0 ),
    global_counts( m_span_count + 1 , 0 );

  // Count unique keys in each span and add requested keys for
  // final total count of keys needed.

  // Append the error check to this communication to avoid
  // and extra reduction operation.
  local_counts[ m_span_count ] = m_span_count != requests.size();

  if ( m_span_count == requests.size() ) {

    for ( size_t i = 0 ; i < m_span_count ; ++i ) {
      local_counts[i] = requests[i] ;
    }

    KeyProcVector::const_iterator j = m_key_usage.begin();

    for ( size_t i = 0 ; i < m_span_count && j != m_key_usage.end() ; ++i ) {
      const KeyType key_span_last = m_key_span[i].second ;
      size_t count = 0 ;
      while ( j != m_key_usage.end() && j->first <= key_span_last ) {
        const KeyType key = j->first ;
        while ( j != m_key_usage.end() && key == j->first ) { ++j ; }
        ++count ;
      }
      local_counts[i] += count ;
    }
  }

#if defined( STK_HAS_MPI )
  if (m_comm_size > 1) {
    BABBLE_STK_PARALLEL_COMM(m_comm, "    calling MPI_Allreduce from generate_new_global_key_upper_bound");
    MPI_Allreduce( (local_counts.empty() ? NULL : & local_counts[0]) , (global_counts.empty() ? NULL : & global_counts[0]) ,
                   m_span_count + 1 , MPI_UNSIGNED_LONG ,
                   MPI_SUM , m_comm );
  }
  else {
    global_counts = local_counts ;
  }
#else
  global_counts = local_counts ;
#endif

  bad_request = global_counts[m_span_count] != 0 ;

  if ( bad_request ) {
    if ( m_span_count != requests.size() ) {
      error_msg << " requests.size() = " << requests.size()
                << " != " << m_span_count << " )" ;
    }
  }

  if ( ! bad_request ) {
    for ( unsigned i = 0 ; i < m_span_count ; ++i ) {
      const size_t span_available =
        ( 1 + m_key_span[i].second - m_key_span[i].first );

      const size_t span_requested = global_counts[i];

      if ( span_available < span_requested ) {
        bad_request = true ;
        error_msg << " global_sum( (existing+request)[" << i << "] ) = "
                  << span_requested
                  << " > global_sum( span_available ) = "
                  << span_available ;
      }
    }
  }

  if ( bad_request ) {
    throw std::runtime_error( error_msg.str() );
  }

  // Determine the maximum generated key

  global_key_upper_bound.resize( m_span_count );

  for ( size_t i = 0 ; i < m_span_count ; ++i ) {
     global_key_upper_bound[i] = m_key_span[i].first + global_counts[i] - 1 ;
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void DistributedIndex::generate_new_keys_local_planning(
  const KeyTypeVector  & key_global_upper_bound ,
  const std::vector<size_t>   & requests_local ,
        LongVector    & new_request ,
        KeyTypeVector  & requested_keys ,
        KeyTypeVector  & contrib_keys ) const
{
  new_request.assign( m_span_count , static_cast<long>(0) );

  contrib_keys.clear();

  KeyProcVector::const_iterator j = m_key_usage.begin();

  for ( size_t i = 0 ; i < m_span_count ; ++i ) {
    // The maximum generated key from any process will
    // not exceed this value.
    const KeyType key_upper_bound = key_global_upper_bound[i] ;

    const size_t init_size = contrib_keys.size();

    const size_t chunk_inc = m_comm_size * DISTRIBUTED_INDEX_CHUNK_SIZE ;

    const size_t chunk_rsize = m_comm_rank * DISTRIBUTED_INDEX_CHUNK_SIZE ;

    for ( KeyType key_begin = m_key_span[i].first +
                  chunk_rsize ;
          key_begin <= key_upper_bound ; key_begin += chunk_inc ) {

      // What is the first key of the chunk
      KeyType key_iter = key_begin ;

      // What is the last key belonging to this process' chunk
      const KeyType key_last =
        std::min( key_begin + DISTRIBUTED_INDEX_CHUNK_SIZE - 1 , key_upper_bound );

      // Jump into the sorted used key vector to
      // the key which may be contributed

      j = std::lower_bound( j, m_key_usage.end(), key_iter, KeyProcLess() );
      // now know:  j == m_key_usage.end() OR
      //            key_iter <= j->first

      for ( ; key_iter <= key_last ; ++key_iter ) {
        if ( j == m_key_usage.end() || key_iter < j->first ) {
          // The current attempt 'key_iter' is not used, contribute it.
          contrib_keys.push_back( key_iter );
        }
        else { // j != m_key_usage.end() && key_iter == j->first
          // The current attempt 'key_iter' is already used,
          // increment the used-iterator to its next key value.
          while ( j != m_key_usage.end() && key_iter == j->first ) {
            ++j ;
          }
        }
      }
    }

    // Determine which local keys will be contributed,
    // keeping what this process could use from the contribution.
    // This can reduce the subsequent communication load when
    // donating keys to another process.

    const size_t this_contrib = contrib_keys.size() - init_size ;

    // How many keys will this process keep:
    const size_t keep = std::min( requests_local[i] , this_contrib );

    // Take the kept keys from the contributed key vector.
    requested_keys.insert( requested_keys.end() ,
                           contrib_keys.end() - keep ,
                           contrib_keys.end() );

    contrib_keys.erase( contrib_keys.end() - keep ,
                        contrib_keys.end() );

    // New request is positive for needed keys or negative for donated keys
    new_request[i] = requests_local[i] - this_contrib ;
  }
}

//----------------------------------------------------------------------

void DistributedIndex::generate_new_keys_global_planning(
  const LongVector    & new_request ,
        LongVector    & my_donations ) const
{
  my_donations.assign( m_comm_size * m_span_count , static_cast<long>(0) );

  // Gather the global request plan for receiving and donating keys
  // Positive values for receiving, negative values for donating.

  LongVector new_request_global( m_comm_size * m_span_count );

#if defined( STK_HAS_MPI )

  if (m_comm_size > 1) { // Gather requests into per-process spans

    // There is a possible bug in MPI_Allgather, for Intel 12;  use MPI_Gather instead
#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1200)
      {
        // MPI doesn't do 'const' in its interface, but the send buffer is const
        void * send_buf = const_cast<long int*>( (new_request.empty() ? NULL : & new_request[0]) );
        void * recv_buf = (new_request_global.empty() ? NULL : & new_request_global[0]) ;
        for (int root = 0; root < m_comm_size; ++root)
          {
            BABBLE_STK_PARALLEL_COMM(m_comm, "    calling MPI_Gather from generate_new_keys_global_planning");
            MPI_Gather( send_buf , m_span_count , MPI_LONG ,
                        recv_buf , m_span_count , MPI_LONG , root, m_comm );
          }
      }
#else
      {
        // MPI doesn't do 'const' in its interface, but the send buffer is const
        void * send_buf = const_cast<long int*>( (new_request.empty() ? NULL : & new_request[0]) );
        void * recv_buf = (new_request_global.empty() ? NULL : & new_request_global[0]) ;
        BABBLE_STK_PARALLEL_COMM(m_comm, "    calling MPI_Allgather from generate_new_keys_global_planning");
        MPI_Allgather( send_buf , m_span_count , MPI_LONG,
                       recv_buf , m_span_count , MPI_LONG , m_comm );
      }
#endif

  }
  else {
    new_request_global.assign(new_request.begin(),new_request.end()) ;
  }
#else
  new_request_global.assign(new_request.begin(),new_request.end()) ;
#endif

  // Now have the global receive & donate plan.
  //--------------------------------------------------------------------
  // Generate my donate plan from the global receive & donate plan.

  for ( unsigned i = 0 ; i < m_span_count ; ++i ) {

    if ( new_request[i] < 0 ) { // This process is donating on this span
      long my_total_donate = - new_request[i] ;

      long previous_donate = 0 ;

      // Count what previous processes have donated:
      for ( int p = 0 ; p < m_comm_rank ; ++p ) {
        const long new_request_p = new_request_global[ p * m_span_count + i ] ;
        if ( new_request_p < 0 ) {
          previous_donate -= new_request_p ;
        }
      }

      // What the donation count will be with my donation:
      long end_donate = previous_donate + my_total_donate ;

      long previous_receive = 0 ;

      // Determine my donation to other processes (one to many).

      for ( int p = 0 ; p < m_comm_size && 0 < my_total_donate ; ++p ) {

        const long new_request_p = new_request_global[ p * m_span_count + i ];

        if ( 0 < new_request_p ) { // Process 'p' receives keys

          // Accumulation of requests:

          previous_receive += new_request_p ;

          if ( previous_donate < previous_receive ) {
            // I am donating to process 'p'
            const long n = std::min( previous_receive , end_donate )
                           - previous_donate ;

            my_donations[ p * m_span_count + i ] = n ;
            previous_donate += n ;
            my_total_donate -= n ;
          }
        }
      }
    }
  }
}

//--------------------------------------------------------------------

void DistributedIndex::generate_new_keys(
  const std::vector<size_t>                 & requests ,
        std::vector< KeyTypeVector > & requested_keys )
{
  //--------------------------------------------------------------------
  // Develop the plan:

  KeyTypeVector global_key_upper_bound ;
  LongVector    new_request ;
  LongVector    my_donations ;
  KeyTypeVector contrib_keys ;
  KeyTypeVector new_keys ;

  // Verify input and generate global sum of
  // current key usage and requested new keys.
  // Throw a parallel consistent exception if the input is bad.

  generate_new_global_key_upper_bound( requests , global_key_upper_bound );

  // No exception thrown means all inputs are good and parallel consistent

  // Determine which local keys will be contributed,
  // keeping what this process could use from the contribution.
  // This can reduce the subsequent communication load when
  // donating keys to another process.

  generate_new_keys_local_planning( global_key_upper_bound ,
                                    requests ,
                                    new_request ,
                                    new_keys ,
                                    contrib_keys );

  // Determine where this process will be donating 'contrib_keys'
  generate_new_keys_global_planning( new_request, my_donations );

  // Due to using an upper bound as opposed to an exact maximum
  // the contrib_keys is likely to contain more keys that are needed.
  // Remove unneeded keys.

  // Backwards to erase from the end
  for ( size_t i = m_span_count ; 0 < i ; ) {
    --i ;
    size_t count = 0 ;
    for ( int p = 0 ; p < m_comm_size ; ++p ) {
      count += my_donations[ p * m_span_count + i ];
    }
    KeyTypeVector::iterator j_beg = contrib_keys.begin();
    KeyTypeVector::iterator j_end = contrib_keys.end();
    j_beg = std::lower_bound( j_beg , j_end , m_key_span[i].first );
    j_end = std::upper_bound( j_beg , j_end , m_key_span[i].second );
    const size_t n = std::distance( j_beg , j_end );
    if ( count < n ) {
      contrib_keys.erase( j_beg + count , j_end );
    }
  }

  // Plan is done, communicate the new keys.
  //--------------------------------------------------------------------
  // Put key this process is keeping into the index.
  m_key_usage.reserve(m_key_usage.size() + new_keys.size());
  for ( KeyTypeVector::iterator i = new_keys.begin();
        i != new_keys.end() ; ++i ) {
    m_key_usage.push_back( KeyProc( *i , m_comm_rank ) );
  }

  //--------------------------------------------------------------------

  CommAll all( m_comm );

  // Sizing

  for ( size_t i = 0 ; i < m_span_count ; ++i ) {
    for ( int p = 0 ; p < m_comm_size ; ++p ) {
      const size_t n_to_p = my_donations[ p * m_span_count + i ];
      if ( 0 < n_to_p ) {
        all.send_buffer(p).skip<KeyType>( n_to_p );
      }
    }
  }

  BABBLE_STK_PARALLEL_COMM(m_comm, "          generate_new_keys calling allocate_buffers")
  all.allocate_buffers( m_comm_size / 4 , false );

  // Packing

  {
    size_t n = 0 ;
    for ( size_t i = 0 ; i < m_span_count ; ++i ) {
      for ( int p = 0 ; p < m_comm_size ; ++p ) {
        const size_t n_to_p = my_donations[ p * m_span_count + i ];
        if ( 0 < n_to_p ) {
          all.send_buffer(p).pack<KeyType>( & contrib_keys[n] , n_to_p );
          for ( size_t k = 0 ; k < n_to_p ; ++k , ++n ) {
            m_key_usage.push_back( KeyProc( contrib_keys[n] , p ) );
          }
        }
      }
    }
  }

  std::sort( m_key_usage.begin() , m_key_usage.end() );

  BABBLE_STK_PARALLEL_COMM(m_comm, "          generate_new_keys calling communicate")
  all.communicate();

  // Unpacking
  unpack_recv_buffer( all, m_comm_size, new_keys);

  stk::util::radix_sort_unsigned((new_keys.empty() ? NULL : &new_keys[0]), new_keys.size());

  requested_keys.resize( m_span_count );

  {
    KeyTypeVector::iterator i_beg = new_keys.begin();
    for ( size_t i = 0 ; i < m_span_count ; ++i ) {
      KeyTypeVector::iterator i_end = i_beg + requests[i] ;
      requested_keys[i].assign( i_beg , i_end );
      i_beg = i_end ;
    }
  }
}

//----------------------------------------------------------------------

} // namespace util
} // namespace stk


