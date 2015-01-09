/*------------------------------------------------------------------------*/
/*      stk : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @file
 * @author H. Carter Edwards
 * @date   January 2007
 */

#include <iostream>
#include <map>
#include <set>
#include <list>
#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <TPI.h>

#include <stk_search/SearchTypes.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_search/OctTreeOps.hpp>

namespace stk_classic {
namespace search {
namespace {

inline unsigned int log2(unsigned int x)
{
    unsigned int l=0;
    if(x >= 1<<16) { x>>=16; l|=16; }
    if(x >= 1<< 8) { x>>= 8; l|= 8; }
    if(x >= 1<< 4) { x>>= 4; l|= 4; }
    if(x >= 1<< 2) { x>>= 2; l|= 2; }
    if(x >= 1<< 1) {         l|= 1; }
    return l;
}

}


//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool hsfc_box_covering( const float     * const global_box ,
                        const float     * const small_box ,
                        stk_classic::OctTreeKey * const covering ,
                        unsigned        &       number,
			double                  scale)
{
  enum { Dimension    = 3 };
  enum { Combinations = 8 };

  const double min = std::numeric_limits<float>::epsilon();
  const double max = 1.0 - min ;

  // Determine the unit-box bounds and bisection depth for the box

  double ubox_low[ Dimension ] ;
  double ubox_up[  Dimension ] ;

  bool valid = true ;

  // Determine unit box and is maximum length
  // The box is bounded by [eps,1-eps].

  double unit_size = 0.0 ;

  for ( unsigned i = 0 ; i < Dimension ; ++i ) {

    const float global_low = global_box[i] ;
    const float global_up  = global_box[i+Dimension] ;
    const float small_low  = small_box[i] ;
    const float small_up   = small_box[i+Dimension] ;

    if ( small_up < global_low ) {
      // Entirely less than 'min'
      ubox_low[i] = ubox_up[i] = min ;
      valid = false ;
    }
    else if ( global_up < small_low ) {
      // Entirely greater than 'max'
      ubox_low[i] = ubox_up[i] = max ;
      valid = false ;
    }
    else {
      double unit_low = ( small_low - global_low ) * scale ;
      double unit_up  = ( small_up  - global_low ) * scale ;

      if ( unit_low < min ) {
        unit_low = min ;
        valid = false ;
      }

      if ( max < unit_up ) {
        unit_up = max ;
        valid = false ;
      }

      if ( unit_up < unit_low ) {
        // A negative volume, consider it a point at the lower
        unit_up = unit_low ;
        valid = false ;
      }
      else {
        const double tmp_size = unit_up - unit_low ;
        if ( unit_size < tmp_size ) { unit_size = tmp_size ; }
      }

      ubox_low[i] = unit_low ;
      ubox_up[i]  = unit_up ;
    }
  }

  // Depth is determined by smallest cell depth
  // that could contain the small_box

  unsigned depth = stk_classic::OctTreeKey::MaxDepth ;

  if ( 0 < unit_size ) {
    const unsigned size_inv = static_cast<unsigned>(1.0 / unit_size);
    depth = log2(size_inv);
    if (depth > stk_classic::OctTreeKey::MaxDepth) depth = stk_classic::OctTreeKey::MaxDepth;
  }

  // Determine the oct-tree nodes for each key

  const unsigned shift    = stk_classic::OctTreeKey::BitsPerWord - depth ;
  const unsigned num_cell = 1 << depth ;

  // At most two cells in each axis at this depth

  unsigned coord_low[ Dimension ];
  unsigned coord_up[  Dimension ];

  for ( unsigned i = 0 ; i < Dimension ; ++i ) {
    const unsigned low = static_cast<unsigned>( ubox_low[i] * num_cell );
    const unsigned up  = static_cast<unsigned>( ubox_up[i]  * num_cell );

    if ( low + 1 < up ) {
      std::string msg("stk_classic::hsfc_box_covering FAILED : depth determination logic error");
      throw std::logic_error( msg );
    }

    coord_low[i] = low << shift ;
    coord_up[i]  = up  << shift ;
  }

  unsigned n = 0 ;

  // Combination 0. No duplicate possible, so pull out of loop.
  covering[n] = hsfc3d( depth , coord_low );
  ++n ;
  
  for ( unsigned i = 1 ; i < Combinations ; ++i ) {

    const bool duplicate = 
      ( ( i & 01 ) && coord_up[0] == coord_low[0] ) ||
      ( ( i & 02 ) && coord_up[1] == coord_low[1] ) ||
      ( ( i & 04 ) && coord_up[2] == coord_low[2] ) ;

    if ( ! duplicate ) {
      unsigned coord[3] ;

      coord[0] = ( i & 01 ) ? coord_up[0] : coord_low[0] ;
      coord[1] = ( i & 02 ) ? coord_up[1] : coord_low[1] ;
      coord[2] = ( i & 04 ) ? coord_up[2] : coord_low[2] ;

      covering[n] = hsfc3d( depth , coord );

      ++n ;
    }
  }

  number = n ;

  return valid ;
}


namespace {

//------------------------------------------//----------------------------------------------------------------------
// Reset the accumulated node weights to only include
// those nodes in the range [ k_first , k_last ]

void accumulate_weights(
  const stk_classic::OctTreeKey &k_node_p ,
  const stk_classic::OctTreeKey &k_first_p ,
  const unsigned ord_end ,
  const unsigned depth ,
        float * const weights )
{
  stk_classic::OctTreeKey k_node (k_node_p);
  stk_classic::OctTreeKey k_first(k_first_p);
  const unsigned ord_node_2 = 2 * oct_tree_offset( depth , k_node );

  if ( k_node.depth() < depth ) {

    double w = 0 ;

    const unsigned d1 = k_node.depth() + 1 ;

    unsigned i = k_first.index( d1 );

    if ( i ) {
      k_node.set_index( d1 , i );

      const unsigned ord = oct_tree_offset( depth , k_node );
      const unsigned ord_2 = ord * 2 ;

      accumulate_weights( k_node , k_first , ord_end , depth , weights );

      // Counts of this node and all of its descending nodes
      w += weights[ord_2] + weights[ ord_2 + 1 ] ;

      k_first = stk_classic::OctTreeKey(); // Done with the lower bound
    }

    for ( ++i ; i <= 8 ; ++i ) {

      k_node.set_index( d1 , i );

      const unsigned ord = oct_tree_offset( depth , k_node );
      const unsigned ord_2 = ord * 2 ;

      if ( ord < ord_end ) {
        accumulate_weights( k_node, k_first , ord_end , depth , weights );

        // Counts of this node and all of its descending nodes
        w += weights[ord_2] + weights[ ord_2 + 1 ] ;
      }
    }

    // Descending node weight

    weights[ ord_node_2 + 1 ] = static_cast<float>(w); 
  }
}

//----------------------------------------------------------------------

void oct_key_split(
  const stk_classic::OctTreeKey & key ,
  const unsigned     upper_ord ,
        stk_classic::OctTreeKey & key_upper )
{
  // Split key at key.depth() + 1

  unsigned d = key.depth();

  key_upper = key ;

  if ( upper_ord == 1 ) { // key_upper gets it all
    while ( d && 1 == key_upper.index(d) ) {
      key_upper.clear_index(d);
      --d ;
    }
  }
  else if ( 8 < upper_ord ) { // key_upper get none of it, Increment key_upper

    unsigned i = 0 ;
    while ( d && 8 == ( i = key_upper.index(d) ) ) {
      key_upper.clear_index(d);
      --d ;
    }
    if ( d ) { key_upper.set_index( d , i + 1 ); }
  }
  else {
    key_upper.set_index( d + 1 , upper_ord );
  }
}

//----------------------------------------------------------------------

void partition( 
  const stk_classic::OctTreeKey & k_first ,
  const unsigned     i_end ,
  const stk_classic::OctTreeKey & key ,
  const unsigned     depth ,
  const float      * weights ,
  const double tolerance ,
  const double target_ratio ,
  double w_lower ,
  double w_upper ,
  stk_classic::OctTreeKey & k_upper )
{
  const unsigned ord_node = oct_tree_offset( depth , key );
  const float * const w_node = weights + ord_node * 2 ;

  const unsigned d1 = key.depth() + 1 ;

  // Add weights from nested nodes and their descendents
  // Try to achieve the ratio.

  const unsigned i_first = k_first.index( d1 );

  unsigned i = ( i_first ) ? i_first : 1 ;
  unsigned j = 8 ;
  {
    stk_classic::OctTreeKey k_upp = key ;
    k_upp.set_index( d1 , j );
    while ( i_end <= oct_tree_offset( depth , k_upp ) ) {
      k_upp.set_index( d1 , --j );
    }
  }

  w_lower += w_node[0] ;
  w_upper += w_node[0] ;

  // At the maximum depth?

  if ( key.depth() == depth ) {
    // Assume weight from unrepresented nested nodes is
    // evenly distributed among the nodes in the span [i,j]

    const unsigned n = 1 + j - i ;

    const double val = static_cast<double>(w_node[1]) / static_cast<double>(n);

    // val = val_lower + val_upper
    // ( w_lower + val_lower ) / ( w_upper + val_upper ) == target_ratio

    const double val_lower =
      ( target_ratio * ( w_upper + val ) - w_lower ) /
      ( target_ratio + 1 ) ;

    if ( 0 < val_lower ) {
      // How much of the range does the lower portion get?
      // Roundoff instead of merely truncating:
      i += static_cast<unsigned>( 0.5 + ( n * val_lower ) / val );

      // Can only get up to the maximum
      if ( j < i ) { i = j ; }
    }
    oct_key_split( key , i , k_upper );
  }
  else {

//    while ( i != j ) {
    while ( i < j ) {
      stk_classic::OctTreeKey ki = key ; ki.set_index( d1 , i );
      stk_classic::OctTreeKey kj = key ; kj.set_index( d1 , j );

      const float * const vi = weights + 2 * oct_tree_offset( depth , ki );
      const float * const vj = weights + 2 * oct_tree_offset( depth , kj );

      const double vali = vi[0] + vi[1] ;
      const double valj = vj[0] + vj[1] ;

      if ( 0 < vali && 0 < valj ) {

        // Choose between ( w_lower += vali ) vs. ( w_upper += valj )
        // Knowing that the skipped value will be revisited.

        if ( ( w_lower + vali ) < target_ratio * ( w_upper + valj ) ) {
          // Add to 'w_lower' and will still need more later
          w_lower += vali ;
          ++i ;
        }
        else {
           // Add to 'w_upper' and will still need more later
          w_upper += valj ;
          --j ;
        }
      }
      else {
        if ( vali <= 0.0 ) { ++i ; }
        if ( valj <= 0.0 ) { --j ; }
      }
    }

    // If 'i' has not incremented then 'k_first' is still in force
    stk_classic::OctTreeKey nested_k_first ;
    if ( i_first == i ) { nested_k_first = k_first ; }

    // Split node nested[i] ?
    stk_classic::OctTreeKey ki = key ; ki.set_index( d1 , i );

    const float * const vi = weights + 2 * oct_tree_offset( depth , ki );
    const double vali = vi[0] + vi[1] ;

    double diff = 0.0 ;

    if ( vali <= 0.0 ) {
      diff = 0.0 ; // Nothing can be done.  Give 'i' to the upper range
    }
    else if ( w_lower < w_upper * target_ratio ) {
      // Try adding to w_lower
      diff = static_cast<double> (w_lower + vali) / static_cast<double>(w_upper) - target_ratio ;
      ++i ;
    }
    else {
      // Try adding to w_upper
      diff = static_cast<double>(w_lower) / static_cast<double>(w_upper + vali) - target_ratio ;
    }

    if ( - tolerance < diff && diff < tolerance ) {
      oct_key_split( key , i , k_upper );
    }
    else {
      partition( nested_k_first , i_end , ki ,
                 depth , weights ,
                 tolerance , target_ratio ,
                 w_lower , w_upper , k_upper );
    }
  }
}

} // namespace <empty>

unsigned processor( const stk_classic::OctTreeKey * const cuts_b ,
                    const stk_classic::OctTreeKey * const cuts_e ,
                    const stk_classic::OctTreeKey & key )
{
  const stk_classic::OctTreeKey * const cuts_p = std::upper_bound( cuts_b , cuts_e , key );

  if ( cuts_p == cuts_b ) {
    std::string msg("stk_classic::processor FAILED: Bad cut-key array");
    throw std::runtime_error(msg);
  }

  return ( cuts_p - cuts_b ) - 1 ;
}

//----------------------------------------------------------------------

void oct_tree_partition_private(
  const unsigned p_first ,
  const unsigned p_end ,
  const unsigned depth ,
  const double   tolerance ,
  float * const weights ,
  const unsigned cuts_length ,
  stk_classic::OctTreeKey * const cuts )
{
  // split tree between [ p_first , p_end )
  const unsigned p_size  = p_end - p_first ;
  const unsigned p_upper = ( p_end + p_first ) / 2 ;

  const double target_fraction =
    static_cast<double> ( p_upper - p_first ) / static_cast<double>(p_size);

  const double target_ratio = target_fraction / ( 1.0 - target_fraction );

  // Determine k_lower and k_upper such that
  //
  // Weight[ k_first , k_lower ] / Weight [ k_upper , k_last ] == target_ratio
  //
  // Within a tollerance

  const stk_classic::OctTreeKey k_first = cuts[ p_first ];

  const unsigned i_end   =
    p_end < cuts_length ? oct_tree_offset( depth , cuts[ p_end ] )
                        : oct_tree_size( depth );

  // Walk the tree [ k_first , k_last ] and accumulate weight

  accumulate_weights( stk_classic::OctTreeKey() , k_first , i_end , depth , weights );

  stk_classic::OctTreeKey k_root ;
  stk_classic::OctTreeKey & k_upper = cuts[ p_upper ] ;

  unsigned w_lower = 0 ;
  unsigned w_upper = 0 ;

  partition( k_first, i_end, k_root ,
             depth, weights,
             tolerance, target_ratio,
             w_lower, w_upper, k_upper );

  const bool nested_lower_split = p_first + 1 < p_upper ;
  const bool nested_upper_split = p_upper + 1 < p_end ;

  // If splitting both lower and upper, and a thread is available
  // then one of the next two calls could be a parallel thread
  // with a local copy of the shared 'weights' array.

  if ( nested_lower_split ) {
    oct_tree_partition_private( p_first, p_upper, depth,
                                tolerance, weights, cuts_length, cuts );
  }

  if ( nested_upper_split ) {
    oct_tree_partition_private( p_upper, p_end, depth,
                                tolerance, weights, cuts_length, cuts );
  }
}

} // namespace search
} // namespace stk_classic
