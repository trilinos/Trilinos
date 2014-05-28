/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_parallel_DistributedIndex_hpp
#define stk_util_parallel_DistributedIndex_hpp

#include <stdint.h>
#include <utility>
#include <vector>
#include <cstddef>
#include <stk_util/parallel/Parallel.hpp>

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
                    const std::vector<KeySpan> & partition_spans );

  /*----------------------------------------*/
  /** \brief  Query with which process the local added keys are shared. */
  void query( std::vector<KeyProc> & sharing_of_local_keys ) const ;

  /** \brief  Query which processors added the given keys.
   *          The local processor is in the output if it
   *          submitted a queried key.
   */
  void query( const std::vector<KeyType> & keys , 
              std::vector<KeyProc> & sharing_of_keys ) const ;

  /** \brief  Query which processors added the given keys.
   *          The results of the query are pushed to the processes
   *          on which the keys are used.
   */
  void query_to_usage( const std::vector<KeyType> & keys , 
                       std::vector<KeyProc> & sharing_of_keys ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Update a parallel index with new and changed keys.
   *          FIRST: Remove this process' participation in the existing keys.
   *          SECOND: Add this process' participation in the new keys.
   */
  void update_keys( const std::vector<KeyType> & add_new_keys ,
                    const std::vector<KeyType> & remove_existing_keys );

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
          std::vector< std::vector<KeyType> > & requested_keys );

private:

  /*------------------------------------------------------------------*/
  /** \brief  An internally used method to count communication needs.
   */
  void generate_new_global_key_upper_bound(
    const std::vector<size_t>  & requests ,
          std::vector<DistributedIndex::KeyType> & global_key_upper_bound ) const;


  /** \brief  An internally used method to determine which
   *          keys will be kept and which will be donated.
   */
  void generate_new_keys_local_planning(
    const std::vector<DistributedIndex::KeyType> & global_key_upper_bound ,
    const std::vector<size_t>  & requests_local ,
          std::vector<long>    & new_requests ,
          std::vector<KeyType> & requested_keys ,
          std::vector<KeyType> & contrib_keys ) const ;

  /** \brief  An internally used method to determine to which
   *          other processes keys will be donated.
   */
  void generate_new_keys_global_planning(
    const std::vector<long>    & new_request ,
          std::vector<long>    & my_donations ) const ;

  /** \brief  An internally used method to query usage of keys. */
  void query( const std::vector<KeyProc> & request ,
                    std::vector<KeyProc> & sharing_of_keys ) const ;

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
  std::vector<KeySpan> m_key_span ;  ///< (min,max) for N span
  std::vector<KeyProc> m_key_usage ; ///< Index for all key usage

  /*  Unit testing of internal methods requires the unit test to have
   *  access to those internal methods.
   */
  friend class ::UnitTestSTKParallelDistributedIndex ;
};

//----------------------------------------------------------------------

} // namespace parallel
} // namespace stk

//----------------------------------------------------------------------

#endif

