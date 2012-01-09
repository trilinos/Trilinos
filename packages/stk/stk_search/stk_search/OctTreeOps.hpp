/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_util_util_OctTreeOps_hpp
#define stk_util_util_OctTreeOps_hpp

#include <utility>
#include <vector>
#include <map>
#include <set>
#include <list>

#include <iostream>
#include <stdexcept>

#include <TPI.h>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/OctTree.hpp>

namespace stk {
namespace search {

//----------------------------------------------------------------------
/**  A recursive kernel used within the oct_tree_partitioning algorithms.
 *   Exposed to support unit testing.
 */

void oct_tree_partition_private(
  const unsigned p_first ,
  const unsigned p_end ,
  const unsigned depth ,
  const double   tolerance ,
  float * const weights ,
  const unsigned cuts_length ,
  OctTreeKey * const cuts );

//----------------------------------------------------------------------
/** Generate an oct-tree covering of a small box within a global box.
 *  The cartesian space is mapped to an oct-tree via Hilbert space
 *  filling curve.  The covering consists of 1..8 oct-tree cells,
 *  the 'covering' array must be dimensioned to at least eight.
 *  Returns true for a "good" small box: it a non-negative volume
 *  and is fully contained within the global box.
 */
bool hsfc_box_covering(
    const float * const global_box ,
    const float * const small_box ,
    OctTreeKey  * const covering ,
    unsigned    &       number,
    double              scale
    );

template <class DomainBoundingBox, class RangeBoundingBox>
void search_tree_statistics( stk::ParallelMachine  arg_comm ,
                             const std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > >  & s ,
                             unsigned * const data )
{
  typedef std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;

  const unsigned huge = std::numeric_limits<unsigned>::max();
  unsigned avg[2] = { 0 , 0 };
  unsigned max[2] = { 0 , 0 };
  unsigned min[2] ;
  min[0] = min[1] = huge ;

  typename SearchTree::const_iterator i ;

  unsigned rcnt = 0;
  unsigned dcnt = 0;
  for ( i = s.begin() ; i != s.end() ; ++i ) {
    const typename SearchTree::value_type  & inode = *i ;
    const unsigned d_size = inode.second.first.size();
    const unsigned r_size = inode.second.second.size();

    if (d_size > 0) {
      avg[0] += d_size ;
      if ( d_size < min[0] ) { min[0] = d_size ; }
      if ( max[0] < d_size ) { max[0] = d_size ; }
      ++dcnt;
    }

    if (r_size > 0) {
      avg[1] += r_size ;
      if ( r_size < min[1] ) { min[1] = r_size ; }
      if ( max[1] < r_size ) { max[1] = r_size ; }
      ++rcnt;
    }
  }

  // Average for this processor

  if (dcnt)
    avg[0] = ( avg[0] + dcnt - 1 ) / dcnt ;

  if ( rcnt ) {
    avg[1] = ( avg[1] + rcnt - 1 ) / rcnt ;
  }

  if ( min[0] == huge ) { min[0] = 0 ; }
  if ( min[1] == huge ) { min[1] = 0 ; }

  all_reduce( arg_comm, ReduceMin<2>( min ) );
  all_reduce( arg_comm, ReduceMax<2>( max ) );
  all_reduce( arg_comm, ReduceSum<2>( avg ) );

  const unsigned p_size = parallel_machine_size( arg_comm );

  // Average among all processors:

  avg[0] = ( avg[0] + p_size - 1 ) / p_size ;
  avg[1] = ( avg[1] + p_size - 1 ) / p_size ;

  data[0] = min[0] ;
  data[1] = max[0] ;
  data[2] = avg[0] ;

  data[3] = min[1] ;
  data[4] = max[1] ;
  data[5] = avg[1] ;
}

//----------------------------------------------------------------------
/** Global bounds for a set of boxes.
 *  The lower bound is the minimum of all boxes, decreased by epsilon.
 *  The upper bound is the maximum of all boxes, increased by epsilon.
 *  Thus all input boxes are fully contained within the global box.
 */
template <class DomainBoundingBox, class RangeBoundingBox>
void box_global_bounds(
  ParallelMachine            arg_comm ,
  const size_t               arg_domain_boxes_number ,
  const DomainBoundingBox * const arg_domain_boxes ,
  const size_t               arg_range_boxes_number ,
  const RangeBoundingBox  * const arg_range_boxes ,
  float        * const arg_global_box )
{
  enum { Dim = 3 };

#if defined(__INTEL_COMPILER) && (__INTEL_COMPILER == 1210)
#pragma warning disable 191
#endif
  const bool symmetric = static_cast<const void * const>(arg_range_boxes) == static_cast<const void * const>(arg_domain_boxes ) ||
    arg_range_boxes == NULL;


  Copy<Dim>( arg_global_box ,         std::numeric_limits<float>::max() );
  Copy<Dim>( arg_global_box + Dim , - std::numeric_limits<float>::max() );

  //------------------------------------
  // Trivial loop threading possible:

  for ( size_t i = 0 ; i < arg_domain_boxes_number ; ++i ) {
    const DomainBoundingBox & box = arg_domain_boxes[i];
    for (int j=0; j<Dim; ++j) {
      if (arg_global_box[j] > static_cast<float>(box.lower(j)) ) {
        arg_global_box[j] = static_cast<float>(box.lower(j));
      }
      if (arg_global_box[j+Dim] < static_cast<float>(box.upper(j)) ) {
        arg_global_box[j+Dim] = static_cast<float>(box.upper(j));
      }
    }
  }

  if ( ! symmetric ) {
    for ( size_t i = 0 ; i < arg_range_boxes_number ; ++i ) {
      const RangeBoundingBox & box = arg_range_boxes[i];
      for (int j=0; j<Dim; ++j) {
        if (arg_global_box[j] > static_cast<float>(box.lower(j)) ) {
          arg_global_box[j] = static_cast<float>(box.lower(j));
        }
        if (arg_global_box[j+Dim] < static_cast<float>(box.upper(j)) ) {
          arg_global_box[j+Dim] = static_cast<float>(box.upper(j));
        }
      }
    }
  }

  //------------------------------------

  all_reduce( arg_comm , ReduceMin<Dim>( arg_global_box ) );
  all_reduce( arg_comm , ReduceMax<Dim>( arg_global_box + Dim ) );

  if (arg_domain_boxes_number == 0 &&
      arg_range_boxes_number  == 0)
    return;

  // Scale up and down by epsilon

  const double eps = std::numeric_limits<float>::epsilon();

  for ( unsigned i = 0 ; i < Dim ; ++i ) {
    float upper = arg_global_box[i+Dim] ;
    float lower = arg_global_box[i] ;

    double delta = eps * ( upper - lower );
    delta = delta > 0 ? delta : eps;

    while ( upper <= arg_global_box[i+Dim] ||
            arg_global_box[i] <= lower ) {
      upper = static_cast<float>( upper + delta );
      lower = static_cast<float>( lower - delta );
      delta *= 2 ;
    }

    arg_global_box[i+Dim] = static_cast<float>(upper);
    arg_global_box[i]     = static_cast<float>(lower);
  }
}

//----------------------------------------------------------------------

template< class S >
struct SetInsertBuffer {
  enum { N = 128 };
  S      & m_set ;
  unsigned m_lock ;
  unsigned m_iter ;
  typename S::value_type m_buffer[ N ] ;

  void overflow();
  void operator()( const typename S::value_type & v );

  SetInsertBuffer( S & s , unsigned l )
    : m_set(s), m_lock(l), m_iter(0) {}

  ~SetInsertBuffer() { overflow(); }
};

template<class S>
void SetInsertBuffer<S>::overflow()
{
  try {
    TPI_Lock( m_lock );
    while ( m_iter ) {
      m_set.insert( m_buffer[ --m_iter ] );
    }
    TPI_Unlock( m_lock );
  }
  catch(...) {
    TPI_Unlock( m_lock );
    throw ;
  }
}

template<class S>
void SetInsertBuffer<S>::operator()( const typename S::value_type & v )
{
  m_buffer[ m_iter ] = v ;
  if ( N == ++m_iter ) { overflow(); }
}


/*
template <class DomainBoundingBox, class RangeBoundingBox>
void proximity_search_symmetric(
  const typename std::map< stk::OctTreeKey,
                           std::pair< std::list< DomainBoundingBox >,
                                      std::list< RangeBoundingBox > > >::const_iterator i_beg ,
  const typename std::map< stk::OctTreeKey,
                           std::pair< std::list< DomainBoundingBox >,
                                      std::list< RangeBoundingBox > > >::const_iterator i_end ,
  SetInsertBuffer< std::set< std::pair< typename DomainBoundingBox::Key,
                                        typename RangeBoundingBox::Key > > > & arg_out )
{
  typedef typename DomainBoundingBox::Key DomainKey;
  typedef typename RangeBoundingBox::Key RangeKey;
  typedef std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;

  typename SearchTree::const_iterator j ;

  typename std::list<DomainBoundingBox>::const_iterator id;
  typename std::list<RangeBoundingBox>::const_iterator ir;

  const typename SearchTree::value_type  & inode = *i_beg ;
  const std::list<DomainBoundingBox> & domain_outer = inode.second.first ;

  const typename std::list<DomainBoundingBox>::const_iterator
    beg_dom_out = domain_outer.begin(),
    end_dom_out = domain_outer.end();

  // Outer cell vs. itself

  for ( id = beg_dom_out ; id != end_dom_out ; ++id ) {
    const DomainBoundingBox & d = *id ;
    for ( ir = id ; ++ir != end_dom_out ; ) {
      const RangeBoundingBox & r = *ir ;
      if ( d.intersect(r) ) {
        const DomainKey & dip = d.key ;
        const RangeKey & rip = r.key ;
        if ( dip < rip ) {
          std::pair<DomainKey,RangeKey> tmp( dip , rip );
          arg_out( tmp );
        }
        else {
          std::pair<Key,Key> tmp( rip , dip );
          arg_out( tmp );
        }
      }
    }
  }

  // Outer cell searching inner cells.
  // Outer cell always precedes inner cells
  // Iterate forward until the cell is not contained.

  const stk::OctTreeKey & outer_key = inode.first ;

  for ( j = i_beg ; ++j != i_end && outer_key.intersect( (*j).first ) ; ) {

    const typename SearchTree::value_type & jnode = *j ;

    const std::list<BoundingBox> & domain_inner = jnode.second.first ;

    const typename std::list<BoundingBox>::const_iterator
      beg_dom_inn = domain_inner.begin(),
      end_dom_inn = domain_inner.end();

    // Check domain_outer vs. domain_inner,
    // skip if the same box.

    for ( id = beg_dom_out ; id != end_dom_out ; ++id ) {
      const BoundingBox & d = *id ;
      for ( ir = beg_dom_inn ; ir != end_dom_inn ; ++ir ) {
        const BoundingBox & r = *ir ;
        if ( d.key != r.key ) {
          if ( d.intersect(r) ) {
            const Key & dip = d.key ;
            const Key & rip = r.key ;
            if ( dip < rip ) {
              std::pair<Key,Key> tmp( dip , rip );
              arg_out( tmp );
            }
            else {
              std::pair<Key,Key> tmp( rip , dip );
              arg_out( tmp );
            }
          }
        }
      }
    }
  }
}
*/

template <class DomainBoundingBox, class RangeBoundingBox>
void proximity_search_asymmetric(
  const typename std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > >::const_iterator i_beg ,
  const typename std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > >::const_iterator i_end ,
  SetInsertBuffer< std::set< std::pair< typename DomainBoundingBox::Key,  typename RangeBoundingBox::Key > > > & arg_out )
{
  typedef typename DomainBoundingBox::Key DomainKey;
  typedef typename RangeBoundingBox::Key RangeKey;
  typedef std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;

  typename SearchTree::const_iterator j ;

  typename std::list<DomainBoundingBox>::const_iterator id;
  typename std::list<RangeBoundingBox>::const_iterator ir;

  const typename SearchTree::value_type & inode = *i_beg ;

  const std::list<DomainBoundingBox> & domain_outer = inode.second.first ;
  const std::list<RangeBoundingBox> & range_outer  = inode.second.second ;

  const typename std::list<DomainBoundingBox>::const_iterator
    beg_dom_out = domain_outer.begin(),
    end_dom_out = domain_outer.end();
  const typename std::list<RangeBoundingBox>::const_iterator
    beg_ran_out = range_outer.begin(),
    end_ran_out = range_outer.end();

  // domain_outer vs. range_outer

  for ( id = beg_dom_out ; id != end_dom_out ; ++id ) {
    const DomainBoundingBox & d = *id ;
    for ( ir = beg_ran_out ; ir != end_ran_out ; ++ir ) {
      const RangeBoundingBox & r = *ir ;
      if ( d.intersect(r) ) {
        std::pair<DomainKey,RangeKey> tmp( d.key , r.key );
        arg_out( tmp );
      }
    }
  }

  // Outer cell searching inner cells.
  // Outer cell always precedes inner cells
  // Iterate forward until the cell is not contained.

  const stk::OctTreeKey & outer_key = inode.first ;

  for ( j = i_beg ; ++j != i_end && outer_key.intersect( (*j).first ) ; ) {

    const typename SearchTree::value_type & jnode = *j ;

    const std::list<RangeBoundingBox> & range_inner  = jnode.second.second ;

    const typename std::list<RangeBoundingBox>::const_iterator
      beg_ran_inn = range_inner.begin(),
      end_ran_inn = range_inner.end();

    // Check domain_outer vs. range_inner

    for ( ir = beg_ran_inn ; ir != end_ran_inn ; ++ir ) {
      const RangeBoundingBox & r = *ir ;
      for ( id = beg_dom_out ; id != end_dom_out ; ++id ) {
        const DomainBoundingBox & d = *id ;
        if ( d.intersect(r) ) {
          std::pair<DomainKey,RangeKey> tmp( d.key , r.key );
          arg_out( tmp );
        }
      }
    }

    // Check domain_inner vs. range_outer if non-symmetric
    const std::list<DomainBoundingBox> & domain_inner = jnode.second.first ;

    const typename std::list<DomainBoundingBox>::const_iterator
      beg_dom_inn = domain_inner.begin(),
      end_dom_inn = domain_inner.end();

    for ( id = beg_dom_inn ; id != end_dom_inn ; ++id ) {
      const DomainBoundingBox & d = *id ;
      for ( ir = beg_ran_out ; ir != end_ran_out ; ++ir ) {
        const RangeBoundingBox & r = *ir ;
        if ( d.intersect(r) ) {
          std::pair<DomainKey,RangeKey> tmp( d.key , r.key );
          arg_out( tmp );
        }
      }
    }
  }
}

unsigned processor( const stk::OctTreeKey * const cuts_b ,
                    const stk::OctTreeKey * const cuts_e ,
                    const stk::OctTreeKey & key );

template <class DomainBoundingBox, class RangeBoundingBox>
void pack(
  CommAll & comm_all ,
  const stk::OctTreeKey * const cuts_b ,
  const std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > & send_tree ,
        std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > * recv_tree )
{
  typedef std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;

  const unsigned p_rank = comm_all.parallel_rank();
  const unsigned p_size = comm_all.parallel_size();
  const stk::OctTreeKey * const cuts_e = cuts_b + p_size ;

  typename SearchTree::const_iterator i ;

  for ( i = send_tree.begin() ; i != send_tree.end() ; ++i ) {
    const stk::OctTreeKey & key = (*i).first ;

    unsigned p = processor( cuts_b , cuts_e , key );

    do {
      if ( p != p_rank ) {
        CommBuffer & buf = comm_all.send_buffer(p);

        const std::list< DomainBoundingBox > & domain = (*i).second.first ;
        const std::list< RangeBoundingBox > & range  = (*i).second.second ;

        typename std::list< DomainBoundingBox >::const_iterator jd ;
        typename std::list< RangeBoundingBox >::const_iterator jr ;

        const unsigned dsize = domain.size();
        const unsigned rsize = range.size();

        buf.pack<unsigned>( key.value() , stk::OctTreeKey::NWord );
        buf.pack<unsigned>( dsize );
        buf.pack<unsigned>( rsize );

        for ( jd = domain.begin() ; jd != domain.end() ; ++jd ) {
          const DomainBoundingBox & box = *jd ;
          buf.pack<DomainBoundingBox>( box );
        }

        for ( jr = range.begin() ; jr != range.end() ; ++jr ) {
          const RangeBoundingBox & box = *jr ;
          buf.pack<RangeBoundingBox>( box );
        }
      }
      else if ( recv_tree ) {
        // Copy contents of the send node
        (*recv_tree)[ key ] = (*i).second ;
      }

      // If the cut keys are at a finer granularity than
      // this key then this key may overlap more than one
      // processor's span.  Check for overlap with the
      // beginning key of the next processor.

      ++p ;

    } while( p < p_size && key.intersect( cuts_b[p] ) );
  }
}

template <class DomainBoundingBox, class RangeBoundingBox>
void unpack(
  CommAll & comm_all ,
  std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > & tree )
{
  typedef std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;

  unsigned domain_size(0) ;
  unsigned range_size(0) ;
  unsigned value[ stk::OctTreeKey::NWord ];
  stk::OctTreeKey key ;
  DomainBoundingBox domain_box ;
  RangeBoundingBox range_box ;

  const unsigned p_size = comm_all.parallel_size();

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = comm_all.recv_buffer(p);

    while ( buf.remaining() ) {
      buf.unpack<unsigned>( value , stk::OctTreeKey::NWord );
      buf.unpack<unsigned>( domain_size );
      buf.unpack<unsigned>( range_size );

      // Insert key, get domain and range

      key.set_value( value );

      typename SearchTree::mapped_type & node = tree[ key ];

      std::list< DomainBoundingBox > & domain = node.first ;
      std::list< RangeBoundingBox > & range  = node.second ;

      for ( unsigned j = 0 ; j < domain_size ; ++j ) {
        buf.unpack<DomainBoundingBox>( domain_box );
        domain.push_back( domain_box );
      }

      for ( unsigned j = 0 ; j < range_size ; ++j ) {
        buf.unpack<RangeBoundingBox>( range_box );
        range.push_back( range_box );
      }
    }
  }
}

template <class DomainBoundingBox, class RangeBoundingBox>
bool communicate(
  stk::ParallelMachine arg_comm ,
  const stk::OctTreeKey * const arg_cuts ,
  const std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > & send_tree ,
        std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > & recv_tree ,
  const bool local_flag )
{
  const unsigned p_size = parallel_machine_size( arg_comm );

  // Communicate search_tree members

  CommAll comm_all( arg_comm );

  // Sizing pass for pack
  pack<DomainBoundingBox, RangeBoundingBox>( comm_all , arg_cuts , send_tree , NULL );

  // If more than 25% then is dense
  const bool global_flag =
    comm_all.allocate_buffers( p_size / 4 , false , local_flag );

  // Actual packing pass, copy local entries too
  pack<DomainBoundingBox, RangeBoundingBox>( comm_all , arg_cuts , send_tree , & recv_tree );

  comm_all.communicate();

  unpack<DomainBoundingBox, RangeBoundingBox>( comm_all , recv_tree );

  return global_flag ;
}


template <class DomainBoundingBox, class RangeBoundingBox>
void communicate(
  stk::ParallelMachine arg_comm ,
  const std::set< std::pair< typename DomainBoundingBox::Key,  typename RangeBoundingBox::Key > > & send_relation ,
        std::set< std::pair< typename DomainBoundingBox::Key,  typename RangeBoundingBox::Key > > & recv_relation )
{
  typedef typename DomainBoundingBox::Key DomainKey;
  typedef typename RangeBoundingBox::Key RangeKey;
  typedef std::pair<DomainKey, RangeKey> ValueType ;

  CommAll comm_all( arg_comm );

  const unsigned p_rank = comm_all.parallel_rank();
  const unsigned p_size = comm_all.parallel_size();

  typename std::set< ValueType >::const_iterator i ;

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) {
    const ValueType & val = *i ;
    if ( val.first.proc == p_rank || val.second.proc == p_rank ) {
      recv_relation.insert( val );
    }
    if ( val.first.proc != p_rank ) {
      CommBuffer & buf = comm_all.send_buffer( val.first.proc );
      buf.skip<ValueType>( 1 );
    }
    if ( val.second.proc != p_rank && val.second.proc != val.first.proc ) {
      CommBuffer & buf = comm_all.send_buffer( val.second.proc );
      buf.skip<ValueType>( 1 );
    }
  }

  // If more than 25% messages then is dense

  comm_all.allocate_buffers( p_size / 4 , false );

  for ( i = send_relation.begin() ; i != send_relation.end() ; ++i ) {
    const ValueType & val = *i ;
    if ( val.first.proc != p_rank ) {
      CommBuffer & buf = comm_all.send_buffer( val.first.proc );
      buf.pack<ValueType>( val );
    }
    if ( val.second.proc != p_rank && val.second.proc != val.first.proc ) {
      CommBuffer & buf = comm_all.send_buffer( val.second.proc );
      buf.pack<ValueType>( val );
    }
  }

  comm_all.communicate();

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = comm_all.recv_buffer( p );
    while ( buf.remaining() ) {
      ValueType val ;
      buf.unpack<ValueType>( val );
      recv_relation.insert( val );
    }
  }
}

//----------------------------------------------------------------------
// Partition a search tree among processors

template <class DomainBoundingBox, class RangeBoundingBox>
void oct_tree_partition(
  stk::ParallelMachine        arg_comm ,
  const std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > & arg_tree ,
  const double                arg_tolerance ,
  std::vector< stk::OctTreeKey > & arg_cuts )
{
  typedef std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;

  enum { tree_depth  = 4 };
  enum { tree_size   = OctTreeSize< tree_depth >::value };
  enum { tree_size_2 = tree_size * 2 };

  const unsigned p_size = parallel_machine_size( arg_comm );
  const stk::OctTreeKey k_null ;

  arg_cuts.assign( p_size , k_null );

  float local_count[  tree_size_2 ];
  float global_count[ tree_size_2 ];

  for ( unsigned i = 0 ; i < tree_size_2 ; ++i ) {
    local_count[i] = 0.0 ;
  }

  for ( typename SearchTree::const_iterator i =  arg_tree.begin() ;
                                   i != arg_tree.end() ; ++i ) {

    const stk::OctTreeKey & key = (*i).first ;

    const std::list< DomainBoundingBox > & domain = (*i).second.first ;
    const std::list< RangeBoundingBox > & range  = (*i).second.second ;

    const unsigned depth   = key.depth();
    const unsigned ordinal = oct_tree_offset( tree_depth , key );
    const unsigned num_d   = domain.size();
    const unsigned num_r   = range.size();
    const unsigned number  = num_d + num_r ;

    if ( depth <= 4 ) { // Values for this node:
      local_count[ 2 * ordinal ] += number ;
    }
    else { // Values for a deeper node
      local_count[ 2 * ordinal + 1 ] += number ;
    }
  }

  all_reduce_sum( arg_comm , local_count , global_count , tree_size_2 );

  oct_tree_partition_private( 0, p_size, tree_depth,
                              arg_tolerance, global_count,
                              p_size, & arg_cuts[0]);
}

template <class DomainBoundingBox, class RangeBoundingBox>
class ProximitySearch {
public:
  typedef typename DomainBoundingBox::Key DomainKey;
  typedef typename RangeBoundingBox::Key RangeKey;
  typedef std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;
  typedef void (*proximity_search_work_routine)(
    const typename SearchTree::const_iterator i_beg ,
    const typename SearchTree::const_iterator i_end ,
    SetInsertBuffer< std::set< std::pair<DomainKey, RangeKey> > > & );

  enum { NLOCKS = 2 };
  enum { GET_LOCK = 0 };
  enum { PUT_LOCK = 1 };

  bool m_symmetric;

  std::set< std::pair<DomainKey, RangeKey> > & m_relation ;
  typename SearchTree::const_iterator m_tree_iter ;
  typename SearchTree::const_iterator m_tree_end ;

  ~ProximitySearch() {}

  ProximitySearch(
    bool symmetric ,
    const SearchTree & search_tree ,
    std::set< std::pair<DomainKey, RangeKey> > & relation );
  void iterate_tree();

private:
  ProximitySearch();
  ProximitySearch( const ProximitySearch & );
  ProximitySearch & operator = ( const ProximitySearch & );
};

template <class DomainBoundingBox, class RangeBoundingBox>
void proximity_search_work( TPI_Work * work )
{
  ProximitySearch<DomainBoundingBox, RangeBoundingBox> * p = (ProximitySearch<DomainBoundingBox, RangeBoundingBox> *) work->info ;
  p->iterate_tree();
}


template <class DomainBoundingBox, class RangeBoundingBox>
void ProximitySearch<DomainBoundingBox, RangeBoundingBox>::iterate_tree()
{
  enum { N_WORK = 32 };

  try {
    SetInsertBuffer< std::set< std::pair<DomainKey, RangeKey> > >
      tmp( m_relation , PUT_LOCK );

    const typename SearchTree::const_iterator i_tree_end = m_tree_end ;

    unsigned n_work = N_WORK ;

    while ( n_work ) {

      n_work = 0 ;

      typename SearchTree::const_iterator i_beg , i_end ;

      // Get work:

      try {
        TPI_Lock( GET_LOCK );

        i_end = i_beg = m_tree_iter ;

        while ( n_work < N_WORK && i_tree_end != i_end ) {
          ++i_end ; ++n_work ;
        }

        m_tree_iter = i_end ;

        TPI_Unlock( GET_LOCK );
      }
      catch( ... ) {
        TPI_Unlock( GET_LOCK );
	throw ;
      }

      // Perform work:

    //  if (m_symmetric) {
    //    for ( ; i_beg != i_end ; ++i_beg ) {
    //      proximity_search_symmetric<DomainBoundingBox,RangeBoundingBox>( i_beg , i_tree_end , tmp );
    //    }
    //  } else {
        for ( ; i_beg != i_end ; ++i_beg ) {
          proximity_search_asymmetric<DomainBoundingBox, RangeBoundingBox>( i_beg , i_tree_end , tmp );
        }
    //  }
    }
  }
  catch ( const std::exception & x ) {
    std::cerr << x.what() << std::endl ;
    std::cerr.flush();
  }
  catch ( ... ) {
    std::cerr << "ProximitySearch::iterate_tree FAILED" << std::endl ;
    std::cerr.flush();
  }
}


template <class DomainBoundingBox, class RangeBoundingBox>
ProximitySearch<DomainBoundingBox, RangeBoundingBox>::ProximitySearch(
  bool                  symmetric ,
  const SearchTree &    search_tree ,
  std::set< std::pair<DomainKey, RangeKey> > & relation )
  : m_symmetric(symmetric),
    m_relation( relation ),
    m_tree_iter( search_tree.begin() ),
    m_tree_end(  search_tree.end() )
{
  if ( m_tree_iter != m_tree_end ) {

    TPI_work_subprogram worker = (TPI_work_subprogram) proximity_search_work<DomainBoundingBox, RangeBoundingBox>;
    TPI_Run_threads(worker, this, NLOCKS );

    if ( m_tree_iter != m_tree_end ) {
      std::string msg("stk::proximity_search FAILED to complete" );
      throw std::runtime_error(msg);
    }
  }
}


//----------------------------------------------------------------------
/** Search for intersection of domain boxes with range boxes
 *  within a given global bounding box.
 *  Output vector of matches with a domain or range box on
 *  the local processor.
 *
 *  If 'arg_cuts' is given it will be used for the parallel
 *  search.  If 'arg_cuts == NULL' then a balanced internal
 *  partitioning will be generated.
 *
 *  The search_tree_stats are for the local search:
 *    [0] = minimum search tree domain cell size
 *    [1] = maximum search tree domain cell size
 *    [2] = average search tree domain cell size
 *    [3] = minimum search tree range cell size
 *    [4] = maximum search tree range cell size
 *    [5] = average search tree range cell size
 *  These statistics require an extra communication to gather.
 *
 *  Returns 'true' if all small boxes on all processors
 *  had non-negative volumes and were fully contained within
 *  the global box.
 */

template <class DomainBoundingBox, class RangeBoundingBox>
bool oct_tree_proximity_search(
  ParallelMachine            arg_comm ,
  const float        * const arg_global_box ,
  const size_t               arg_domain_boxes_number ,
  const DomainBoundingBox * const arg_domain_boxes ,
  const size_t               arg_range_boxes_number ,
  const RangeBoundingBox * const arg_range_boxes ,
  const OctTreeKey   * const arg_cuts ,
  std::vector< std::pair< typename DomainBoundingBox::Key,  typename RangeBoundingBox::Key > > & arg_relation ,
  unsigned * const arg_search_tree_stats = NULL )
{
  typedef typename DomainBoundingBox::Key DomainKey;
  typedef typename RangeBoundingBox::Key RangeKey;
  typedef std::map< stk::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;

  enum { Dim = 3 };

  //const bool symmetric = static_cast<const void * const>(arg_range_boxes) == static_cast<const void *const>(arg_domain_boxes ) ||
  //  arg_range_boxes == NULL;
  const bool symmetric = false;

  const unsigned p_size = parallel_machine_size( arg_comm );
  const unsigned p_rank = parallel_machine_rank( arg_comm );

  //----------------------------------------------------------------------
  // Search tree defined by oct-tree covering for boxes

  bool local_violations = false ;
  bool global_violations = false ;

  SearchTree search_tree ;

  {
    stk::OctTreeKey covering[8] ;
    unsigned   number = 0 ;

    double scale = arg_global_box[0+Dim] - arg_global_box[0];
    for ( unsigned i = 1 ; i < Dim ; ++i ) {
      double tst_scale = arg_global_box[i+Dim] - arg_global_box[i];
      if (tst_scale > scale) scale = tst_scale;
    }
    if (scale > 0.0) // Not an error. Could arise with a point bounding box in the range/domain...
      scale = 1.0 / scale;
    else
      scale = 1.0;

    for ( size_t i = 0 ; i < arg_domain_boxes_number ; ++i ) {

      DomainBoundingBox tmp( arg_domain_boxes[i] );

      tmp.key.proc = p_rank ;

      float box[6];

      for (int x=0; x<Dim; ++x) {
        box[x] = tmp.lower(x);
        box[x+Dim] = tmp.upper(x);
      }

      const bool valid =
        hsfc_box_covering( arg_global_box, box, covering, number, scale );

      if ( ! valid ) { local_violations = true ; }

      for ( unsigned k = 0 ; k < number ; ++k ) {
        const stk::OctTreeKey key = covering[k] ;
        search_tree[key].first.push_back(tmp);
      }
    }

    if ( ! symmetric ) {
      for ( size_t i = 0 ; i < arg_range_boxes_number ; ++i ) {

        RangeBoundingBox tmp( arg_range_boxes[i] );

        tmp.key.proc = p_rank ;
        float box[6];

        for (int x=0; x<Dim; ++x) {
          box[x] = tmp.lower(x);
          box[x+Dim] = tmp.upper(x);
        }

        const bool valid =
          hsfc_box_covering( arg_global_box, box, covering, number, scale );

        if ( ! valid ) { local_violations = true ; }

        for ( unsigned k = 0 ; k < number ; ++k ) {
          const stk::OctTreeKey key = covering[k] ;
          search_tree[key].second.push_back(tmp);
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // Use a set to provide a unique and sorted result.

  std::set< std::pair<DomainKey, RangeKey> > tmp_relation ;

  if ( p_size == 1 ) {

    global_violations = local_violations ;

    if ( arg_search_tree_stats ) {
      search_tree_statistics( arg_comm , search_tree ,
                              arg_search_tree_stats );
    }
    ProximitySearch<DomainBoundingBox, RangeBoundingBox>( symmetric, search_tree, tmp_relation);
  }
  else {
    // Communicate search_tree members

    SearchTree local_tree ;

    std::set< std::pair<DomainKey, RangeKey> > local_relation ;

    if ( arg_cuts ) {
      global_violations =
        communicate<DomainBoundingBox, RangeBoundingBox>( arg_comm , arg_cuts , search_tree , local_tree ,
                           local_violations );
    }
    else {
      const double tolerance = 0.001 ;

      std::vector< stk::OctTreeKey > cuts ;

      oct_tree_partition( arg_comm , search_tree , tolerance , cuts );

      global_violations =
        communicate<DomainBoundingBox, RangeBoundingBox>(arg_comm , & cuts[0] , search_tree , local_tree ,
                          local_violations );
    }

    // Local proximity search with received members

    if ( arg_search_tree_stats ) {
      search_tree_statistics( arg_comm , local_tree ,
                              arg_search_tree_stats );
    }

    ProximitySearch<DomainBoundingBox, RangeBoundingBox>( symmetric, local_tree, local_relation);

    // Communicate relations back to domain and range processors

    communicate<DomainBoundingBox, RangeBoundingBox>( arg_comm , local_relation , tmp_relation );
  }

  arg_relation.clear();
  arg_relation.reserve( tmp_relation.size() );

  typename std::set< std::pair<DomainKey, RangeKey> >::iterator ir ;
  for ( ir = tmp_relation.begin() ; ir != tmp_relation.end() ; ++ir ) {
    arg_relation.push_back( *ir );
  }
  return global_violations ;
}

//----------------------------------------------------------------------

} // namespace search
} // namespace stk

#endif

