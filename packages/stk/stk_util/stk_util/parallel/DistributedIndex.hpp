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

#ifndef stk_util_parallel_DistributedIndex_hpp
#define stk_util_parallel_DistributedIndex_hpp

#include <stk_util/stk_config.h>
#include <stdint.h>                     // for uint64_t
#include <cstddef>                      // for size_t
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <utility>                      // for pair
#include <vector>                       // for vector

class UnitTestSTKParallelDistributedIndex ;

namespace stk {
namespace parallel {

/** \brief Parallel cross-reference index for a collection of keys.
 *
 * \section  Parallel operations
 *
 *  All methods use parallel collective communication operations.
 *
 *  Each processor constructs a DistributedIndex with its
 *  local collection of keys.
 *  The resulting DistributedIndex may be queried for
 *  - which other processors input the same local keys or
 *  - which processors submitted an arbitrary set of keys.
 *
 * \section  Partitioned key space
 *
 *  When construction the key space is partitioned into N spans
 *  defined by a minimum and maximum key value for each span.
 */
class DistributedIndex {
public:
  typedef uint64_t                    KeyType ;
  typedef int                         ProcType ;
  typedef std::pair<KeyType,KeyType>  KeySpan ;
  typedef std::pair<KeyType,ProcType> KeyProc ;

  typedef std::vector<unsigned long> UnsignedVector;
  typedef std::vector<long>          LongVector;
  typedef std::vector<KeySpan>       KeySpanVector;
  typedef std::vector<KeyType>       KeyTypeVector;
  typedef std::vector<KeyProc>       KeyProcVector;

  /*----------------------------------------*/

  ~DistributedIndex();

  /** \brief  Construct a parallel index with a parititioning of the key space.
   *
   *  To guarantee parallel consistency process zero broadcasts
   *  the partition bounds and all other processes accept those bounds.
   *
   *  The vector of spans must be well-ordered and not overlap; e.g.,
   *    partition_spans[i].first  <= partition_spans[i].second
   *    partition_spans[i].second <  partition_spans[i+1].first
   */
  DistributedIndex( ParallelMachine comm ,
                    const KeySpanVector & partition_spans );

  /*----------------------------------------*/
  /** \brief  Query with which process the local added keys are shared. */
  void query( KeyProcVector & sharing_of_local_keys ) const ;

  /** \brief  Query which processors added the given keys.
   *          The local processor is in the output if it
   *          submitted a queried key.
   */
  void query( const KeyTypeVector & keys, KeyProcVector & sharing_of_keys ) const ;

  /** \brief  Query which processors added the given keys.
   *          The results of the query are pushed to the processes
   *          on which the keys are used.
   */
  void query_to_usage( const KeyTypeVector & keys , KeyProcVector & sharing_of_keys ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Update a parallel index with new and changed keys.
   *          FIRST: Remove this process' participation in the existing keys.
   *          SECOND: Add this process' participation in the new keys.
   */
  void update_keys( KeyTypeVector::const_iterator  add_new_keys_begin, KeyTypeVector::const_iterator  add_new_keys_end ,
                    KeyTypeVector::const_iterator remove_existing_keys_begin, KeyTypeVector::const_iterator remove_existing_keys_end);

  void update_keys( KeyTypeVector::const_iterator add_new_keys_begin, KeyTypeVector::const_iterator add_new_keys_end);

  void register_removed_key( KeyType removed_key )
  { m_removed_keys.push_back(removed_key); }

  /** \brief  Request a collection of unused keys.
   *
   *  Each process inputs its independent request for keys which
   *  do not currenly appear in the distributed index.
   *  The output keys are grouped by partition_spans.
   *  The output keys are guaranteed to be unique among all processes.
   *  The output keys are added into the distributed index.
   *
   *  Multiple request should be bundled to reduce
   *  parallel communication costs.
   *  The output 'requested_keys' are sorted according to the policy.
   *
   *  The the 'first' member of the requests are the lower bound
   *  value for the keys.
   *
   *  \throw  Throw an exception on all process if any request
   *          cannot be satisfied.
   */
  void generate_new_keys(
    const std::vector<size_t>                 & requests ,
          std::vector< KeyTypeVector > & requested_keys );

private:

  /*------------------------------------------------------------------*/
  /** \brief  An internally used method to count communication needs.
   */
  void generate_new_global_key_upper_bound(
    const std::vector<size_t>  & requests , KeyTypeVector & global_key_upper_bound ) const;


  /** \brief  An internally used method to determine which
   *          keys will be kept and which will be donated.
   */
  void generate_new_keys_local_planning(
    const KeyTypeVector & global_key_upper_bound ,
    const std::vector<size_t>  & requests_local ,
          LongVector    & new_requests ,
          KeyTypeVector & requested_keys ,
          KeyTypeVector & contrib_keys ) const ;

  /** \brief  An internally used method to determine to which
   *          other processes keys will be donated.
   */
  void generate_new_keys_global_planning(
    const LongVector    & new_request ,
          LongVector    & my_donations ) const ;

  /** \brief  An internally used method to query usage of keys. */
  void query( const KeyProcVector & request ,
                    KeyProcVector & sharing_of_keys ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  An internally used method to determine to which
   *          process is responsible for knowing about the key's usage.
   */
  ProcType to_which_proc( const KeyType & key ) const ;

  /*------------------------------------------------------------------*/
  /*  Disable default construction and copies. */

  DistributedIndex();
  DistributedIndex( const DistributedIndex & );
  DistributedIndex & operator = ( const DistributedIndex & );

  /*------------------------------------------------------------------*/

  ParallelMachine      m_comm ;      ///< Distributed communicator
  ProcType             m_comm_rank ; ///< This process' rank
  ProcType             m_comm_size ; ///< The number of processes
  size_t               m_span_count ;///< Number of spans of keys
  KeySpanVector m_key_span ;  ///< (min,max) for N span
  KeyTypeVector m_removed_keys;

  /*  Unit testing of internal methods requires the unit test to have
   *  access to those internal methods.
   */
  friend class ::UnitTestSTKParallelDistributedIndex ;

protected:
  KeyProcVector m_key_usage ; ///< Index for all key usage

};

//----------------------------------------------------------------------

} // namespace parallel
} // namespace stk

//----------------------------------------------------------------------

#endif

