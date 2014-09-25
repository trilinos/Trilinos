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

namespace stk {
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
                        stk::OctTreeKey * const covering ,
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

  unsigned depth = stk::OctTreeKey::MaxDepth ;

  if ( 0 < unit_size ) {
    const unsigned size_inv = static_cast<unsigned>(1.0 / unit_size);
    depth = log2(size_inv);
    if (depth > stk::OctTreeKey::MaxDepth) depth = stk::OctTreeKey::MaxDepth;
  }

  // Determine the oct-tree nodes for each key

  const unsigned shift    = stk::OctTreeKey::BitsPerWord - depth ;
  const unsigned num_cell = 1 << depth ;

  // At most two cells in each axis at this depth

  unsigned coord_low[ Dimension ];
  unsigned coord_up[  Dimension ];

  for ( unsigned i = 0 ; i < Dimension ; ++i ) {
    const unsigned low = static_cast<unsigned>( ubox_low[i] * num_cell );
    const unsigned up  = static_cast<unsigned>( ubox_up[i]  * num_cell );

    if ( low + 1 < up ) {
      std::string msg("stk::hsfc_box_covering FAILED : depth determination logic error");
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

//----------------------------------------------------------------------

void getNonZeroWeightsAndOffsets( float const * const weights, unsigned tree_size, float &totalWeight, std::vector<float> &nonZeroWeights, std::vector<unsigned> &ordinalsOfNonZeroWeights)
{
    for (size_t i=0;i<tree_size*2;i+=2)
    {
        unsigned index=i/2;
        float totalWeightAtOrdinal = weights[i]+weights[i+1];
        if ( totalWeightAtOrdinal > 0 )
        {
            nonZeroWeights.push_back(totalWeightAtOrdinal);
            totalWeight += totalWeightAtOrdinal;
            ordinalsOfNonZeroWeights.push_back(index);
        }
    }
}

void partition_oct_tree(unsigned numProcsLocal, unsigned depth, const float * const weights, unsigned cuts_length, stk::OctTreeKey *cuts)
{
    std::vector<float> nonZeroWeights;
    std::vector<unsigned> offsetsOfNonZeroWeights;

    unsigned tree_size = stk::oct_tree_size(depth);
    float totalWeight=0;
    getNonZeroWeightsAndOffsets(weights, tree_size, totalWeight, nonZeroWeights, offsetsOfNonZeroWeights);

    // Following code is taking the total weight on all nodes of the tree and calculating an average weight per processor.
    // Using that it is determining the ordinal of the node on where to cut, and then switching the ordinal to a key
    // that is being stored in the cuts array. Keys are what are used outside this function for the "physical" tree.
    // Note: Partition assumes depth 4 while physical tree can be multiple depths.
    // the weights array is 2*tree_size where tree_size is associated with depth 4 (unless the enum for depth has changed)
    // weights[2*ordinal] = the number of boxes this node has. weights[2*ordinal+1] = the number of boxes

    float targetWeightPerProc = totalWeight/numProcsLocal;
    unsigned indexOfLastNodeInTree = tree_size - 1;
    float totalAccumulatedWeight = 0;

    unsigned procCounter = 0;
    cuts[procCounter] = stk::OctTreeKey();
    procCounter++;

    unsigned weightIndex = 0;
    while ( totalAccumulatedWeight < totalWeight && procCounter < numProcsLocal)
    {
        float accumulatedWeightForProc = 0;
        unsigned offset = 0;
        while ( accumulatedWeightForProc < targetWeightPerProc && weightIndex < nonZeroWeights.size() )
        {
            accumulatedWeightForProc += nonZeroWeights[weightIndex];
            offset = offsetsOfNonZeroWeights[weightIndex] + 1;
            if ( offset >= tree_size )
            {
                offset = tree_size-1;
            }
            weightIndex++;
        }
        stk::search::calculate_key_using_offset(depth, offset, cuts[procCounter]);
        totalAccumulatedWeight += accumulatedWeightForProc;
        procCounter++;
    }

    for (unsigned i=procCounter; i<numProcsLocal; i++)
    {
        stk::search::calculate_key_using_offset(depth, indexOfLastNodeInTree, cuts[i]);
    }
}

unsigned processor( const stk::OctTreeKey * const cuts_b ,
                    const stk::OctTreeKey * const cuts_e ,
                    const stk::OctTreeKey & key )
{
  const stk::OctTreeKey * const cuts_p = std::upper_bound( cuts_b , cuts_e , key );

  if ( cuts_p == cuts_b ) {
    std::string msg("stk::processor FAILED: Bad cut-key array");
    throw std::runtime_error(msg);
  }

  return ( cuts_p - cuts_b ) - 1 ;
}

//----------------------------------------------------------------------


void calculate_key_using_offset(unsigned depth, unsigned offset, stk::OctTreeKey &kUpper)
{
    stk::OctTreeKey localKey;
    for (int depthLevel = depth-1; depthLevel>=0; depthLevel--)
    {
        unsigned max = oct_tree_size(depthLevel);
        unsigned index = (offset-1)/max;
        if ( index > 0 && offset > 0)
        {
            offset -= (max*index + 1);
            unsigned bitToSetIndexOn = depth-depthLevel;
            index++;
            localKey.set_index(bitToSetIndexOn, index);
        }
        else if ( offset > 0 )
        {
            offset--;
            localKey.set_index(depth-depthLevel, 1);
        }
    }
    kUpper = localKey;
}

} // namespace search
} // namespace stk
