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

/*! \file Zoltan2_PartitioningSolution.cpp
    \brief Helper functions for partitioning solution, including simple remapping of part numbers to minimize migration cost
 */
#include "Zoltan2_PartitioningSolution.hpp"

namespace Zoltan2 {


////////////////////////////////////////////////////////////////////
long measure_stays(
  partId_t *remap,
  int *idx,
  partId_t *adj,
  long *wgt,
  partId_t nGlobalParts,
  partId_t maxPartsWithEdges
)
{
// Return the weight of objects staying with a given remap.
// If remap is NULL, compute weight of objects staying with given partition
  long staying = 0;
  for (partId_t i = 0; i < maxPartsWithEdges; i++) {
    partId_t k = (remap ? remap[i] : i);
    for (partId_t j = idx[i]; j < idx[i+1]; j++) {
      if (k == (adj[j]-maxPartsWithEdges)) {
        staying += wgt[j];
        break;
      }
    }
  }
  return staying;
}

////////////////////////////////////////////////////////////////////
//
// Greedy algorithm for maximum-weight matching.
// This is an 1/2-approximation, but requires a sort.
// We could also use the Path Growing Algorithm by
// Drake & Hogardy, which runs in linear time.

#include <vector>
#include <algorithm>

// This struct should be local to matching(), but compiler complains.
typedef struct {
  partId_t i;
  partId_t j;
  long val;
} triplet; // edge (i,j,val)

static bool compare_triplets(triplet a, triplet b)
{
  return (a.val > b.val); // descending order
}


partId_t matching(
  int *idx,
  partId_t *adj,
  long *wgt,
  partId_t tnVtx,
  partId_t *match
)
{
  partId_t nmatch=0;
  std::vector<triplet> edges(idx[tnVtx]);

  // Make vector of triplets (edges)
  size_t k=0;
  for (int i=0; i<tnVtx; i++){
    for (int jj=idx[i]; jj<idx[i+1]; jj++){
      int j = adj[jj];
      if (i<=j){ // We only need each edge once.
        edges[k].i = i;
        edges[k].j = j;
        edges[k].val = wgt[k];
      }
      k++;
    }
  }

  // Sort edges by weight
  std::sort (edges.begin(), edges.end(), compare_triplets);

  // Greedy loop over sorted edges
  // std::cout << "After sort:" << std::endl;
  for (std::vector<triplet>::iterator it=edges.begin(); it!=edges.end(); ++it){
    // std::cout << "edge (" << it->i << ", " << it->j << ", " << it->val << ")" << std::endl
;

    if ((match[it->i] == it->i) && (match[it->j] == it->j )){
      match[it->i] = it->j;
      match[it->j] = it->i;
      nmatch++;
    }
  }
  return nmatch;
}
}


