// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_GreedyMWM.hpp
    \brief Greedy Maximal Weight Matching
 */

#ifndef ZOLTAN2_GREEDYMWM_HPP__
#define ZOLTAN2_GREEDYMWM_HPP__
namespace Zoltan2 {

////////////////////////////////////////////////////////////////////
//
// Greedy algorithm for maximum-weight matching.
// This is an 1/2-approximation, but requires a sort.
// Linear-time approximation algorithms should be considered
// in the future, e.g., the Path Growing Algorithm by
// Drake & Hogardy, and the Suitor algorithm by Manne & Halappanavar.
// The algorithm runs in serial; the graph must be gathered to 
// the process before entry.
////////////////////////////////////////////////////////////////////

#include <vector>
#include <algorithm>

// This struct should be local to GreedyMWM(), but compiler complains.
template <typename vtx_t, typename wgt_t> 
struct GMWM_triplet{
  vtx_t i;
  vtx_t j;
  wgt_t val;
};

template <typename vtx_t, typename wgt_t>
static bool compare_triplets(GMWM_triplet<vtx_t,wgt_t> a,
                             GMWM_triplet<vtx_t,wgt_t> b)
{
  return (a.val > b.val); // descending order
}


template <typename vtx_t, typename wgt_t>
vtx_t GreedyMWM(
  int *idx,         // Index into compressed sparse edge list; 
                    // idx[i] is index into adj of first edge for vertex i
  vtx_t *adj,       // Edge list in compressed sparse format
  wgt_t *wgt,       // weights for each edge
  vtx_t tnVtx,      // number of vertices
  vtx_t *match      // output result:  vtx i matches with vtx match[i]
)
{
  typedef GMWM_triplet<vtx_t, wgt_t> triplet_t;
  vtx_t nmatch=0;
  std::vector<triplet_t> edges(idx[tnVtx]);

  // Make vector of triplets (edges)
  size_t k=0;
  for (int i=0; i<tnVtx; i++){
    for (int jj=idx[i]; jj<idx[i+1]; jj++){
      int j = adj[jj];
      if (i<=j){ // We need each edge only once.
        edges[k].i = i;
        edges[k].j = j;
        edges[k].val = wgt[k];
      }
      k++;
    }
  }

  // Sort edges by weight
  std::sort (edges.begin(), edges.end(), compare_triplets<vtx_t,wgt_t>);

  // Greedy loop over sorted edges
  // std::cout << "After sort:" << std::endl;
  for (typename std::vector<triplet_t>::iterator it=edges.begin();
       it!=edges.end(); ++it){

    // std::cout << "edge (" << it->i << ", " << it->j << ", " 
    //           << it->val << ")" << std::endl;

    if ((match[it->i] == it->i) && (match[it->j] == it->j )){
      match[it->i] = it->j;
      match[it->j] = it->i;
      nmatch++;
    }
  }
  return nmatch;
}
}


#endif
