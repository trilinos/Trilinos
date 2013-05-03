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

/*! \file Zoltan2_Remap.hpp
    \brief Simple remapping of part numbers to minimize migration cost
 */
#include "Teuchos_CommHelpers.hpp"

#ifndef _ZOLTAN2_REMAP_HPP_
#define _ZOLTAN2_REMAP_HPP_


namespace Zoltan2 {


////////////////////////////////////////////////////////////////////
static size_t measure_stays(
  partId_t *remap,
  int *idx,
  partId_t *adj,
  long *wgt,
  partId_t nGlobalParts
)
{
// Return the weight of objects staying with a given remap.
// If remap is NULL, compute weight of objects staying with given partition
  size_t staying = 0;
  for (partId_t i = 0; i < nGlobalParts; i++) {
    partId_t k = (remap ? remap[i] : i);
    for (int j = idx[i]; j < idx[i+1]; j++) {
      if (adj[j] == k) {
        staying += wgt[j];
        break;
      }
    }
  }
  return staying;
}

////////////////////////////////////////////////////////////////////
static partId_t matching(
  int *idx,
  partId_t *adj,
  long *wgt,
  partId_t tnVtx,
  partId_t *match
)
{
  // TODO:  Add Erik's code here
  return 1;
}

////////////////////////////////////////////////////////////////////
// Remap a new part assignment vector for maximum overlap with an input
// part assignment.
//
// Assumptions for this version:
//   input part assignment == processor rank for every local object.
//   assuming nGlobalParts <= num ranks
// TODO:  Write a version that takes the input part number as input;
//        this change requires input parts in adapters to be provided in
//        the Adapter.
// TODO:  For repartitioning, compare to old remapping results; see Zoltan1.

template <typename Adapter>
static void RemapParts(
  ArrayRCP<partId_t> &parts, // Array of computed part assignments; one per obj
  partId_t nGlobalParts,       // Requested number of parts
  RCP<Comm<int> > &comm        // the Problem's communicator
)
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;

  size_t len = parts.size();

  partId_t me = comm->getRank();
  int np = comm->getSize();

  // Build edges of a bipartite graph with nGlobalParts vertices in each set.
  // Weight edge[parts[i]] by the number of objects that are going from
  // this rank to parts[i].
  // We use a map, assuming the number of unique parts in the parts array
  // is small to keep the binary search efficient.
  // TODO We use the count of objects to move; should change to SIZE of objects
  // to move; need SIZE function in Adapter.

  std::map<partId_t, size_t> edges;
  size_t lstaying = 0;  // Total num of local objects staying if we keep the
                        // current mapping. TODO:  change to SIZE of local objs
  size_t gstaying = 0;  // Total num of objects staying in the current partition

  if (me < nGlobalParts) {
    for (size_t i = 0; i < len; i++) {
      edges[parts[i]]++;                // TODO Use obj size instead of count
      if (parts[i] == me) lstaying++;    // TODO Use obj size instead of count
    }
  }
  else {
    // If me >= nGlobalPart, all of me's data has to move to new part #.
    // No need to include in the matching.
  }

  Teuchos::reduceAll<int, size_t>(*comm, Teuchos::REDUCE_SUM, 1,
                                  &lstaying, &gstaying);
//TODO  if (gstaying == Adapter::getGlobalNumObjs()) return;  // Nothing to do

  partId_t *remap = new partId_t[nGlobalParts];
  for (partId_t i = 0; i < nGlobalParts; i++) remap[i] = -1;

  int nedges = edges.size();

  // Accumulate on rank 0:  2* (number of edges on each rank)
  int *idx = NULL;
  int *sizes = NULL;
  if (me == 0) {
    idx = new int[np+nGlobalParts+1];
    sizes = new int[np];
  }
  if (np > 1)
    Teuchos::gather<int, int>(&nedges, 1, sizes, 1, 0, *comm);
  else
    sizes[0] = nedges;

  // prefix sum to build the idx array
  partId_t maxPartsWithEdges = 0;
  if (me == 0) {
    idx[0] = 0;
    for (int i = 0; i < np; i++) {
      idx[i+1] += sizes[i];
      if (sizes[i]) maxPartsWithEdges = i;
    }
    maxPartsWithEdges++;
  }

  // prepare to send edges
  partId_t *bufv = new partId_t[nedges];
  size_t *bufw = new size_t[nedges];
  // Create buffer with edges (me, part[i]) and weight edges[parts[i]].
  int cnt = 0;
  for (std::map<partId_t, size_t>::iterator it = edges.begin();
       it != edges.end(); it++) {
    bufv[cnt] = it->first;  // target part
    bufw[cnt] = it->second; // weight
    cnt++;
  }

  // Prepare to receive edges on rank 0
  partId_t *adj;
  long *wgt;
  if (me == 0) {
    adj = new partId_t[2*idx[np]];  // need 2x space to symmetrize later
    wgt = new long[2*idx[np]];  // need 2x space to symmetrize later
  }

  Teuchos::gatherv<int, partId_t>(bufv, cnt, adj, sizes, idx, 0, *comm);
  Teuchos::gatherv<int, long>(bufw, cnt, wgt, sizes, idx, 0, *comm);
  delete [] bufv;
  delete [] bufw;

  // Now have constructed graph on rank 0.
  // Call the matching algorithm

  int doRemap;
  if (me == 0) {
    // We have the "LHS" vertices of the bipartite graph; need to create
    // "RHS" vertices and add symmetric edges.
    for (int i = 0; i < np; i++)
      sizes[i] = 0;  // Reuse of sizes array assumes nGlobalParts <= np

    for (int i = 0; i < idx[np]; i++) {
      sizes[adj[i]]++;  // Number of edges with adj[i] as target
      adj[i] += maxPartsWithEdges;  // New RHS vertex number;
                                        // offset by num LHS vertices
    }

    // Build idx for RHS vertices
    partId_t tnVtx = maxPartsWithEdges + nGlobalParts;  // total # vertices
    for (partId_t i = maxPartsWithEdges; i < tnVtx; i++) {
      idx[i+1] = idx[i] + sizes[i-maxPartsWithEdges];
    }

    // Add edges from RHS vertices to LHS vertices
    for (int i = 0; i < np; i++)
      sizes[i] = 0;  // Reuse of sizes array assumes nGlobalParts <= np

    for (partId_t i = 0; i < maxPartsWithEdges; i++) {
      for (int j = idx[i]; j < idx[i+1]; j++) {
        partId_t tgt = adj[j];
        partId_t stgt = tgt - maxPartsWithEdges;
        adj[idx[tgt]+sizes[stgt]] = i;
        wgt[idx[tgt]+sizes[stgt]] = wgt[j+1];
        sizes[stgt]++;
      }
    }

    // Perform matching on the graph
    partId_t *match = new partId_t[tnVtx];
    for (partId_t i = 0; i < tnVtx; i++) match[i] = i;
    partId_t nmatches = matching(idx, adj, wgt, tnVtx, match);

    // Process the matches
    bool *used = new bool[nGlobalParts];
    for (partId_t i = 0; i < nGlobalParts; i++) used[i] = false;
    if (nmatches) {

      // First, process all matched parts
      for (partId_t i = 0; i < nGlobalParts; i++) {
        partId_t tmp = i + maxPartsWithEdges;
        if (match[tmp] != tmp) {
          remap[i] = match[tmp];
          used[match[tmp]] = true;
        }
      }

      // Second, process unmatched parts; keep same part number if possible
      for (partId_t i = 0; i < nGlobalParts; i++) {
        if (remap[i] > -1) continue;
        if (!used[i]) {
          remap[i] = i;
          used[i] = true;
        }
      }

      // Third, process unmatched parts; give them the next unused part
      for (partId_t i = 0, uidx = 0; i < nGlobalParts; i++) {
        if (remap[i] > -1) continue;
        while (used[uidx]) uidx++;
        remap[i] = uidx;
        used[uidx] = true;
      }
    }

    delete [] match;
    delete [] used;

    size_t newgstaying = measure_stays(remap, idx, adj, wgt, nGlobalParts);
    doRemap = (newgstaying > gstaying);
    cout << "gstaying " << gstaying << " measure(input) "
         << measure_stays(NULL, idx, adj, wgt, nGlobalParts)
         << " newgstaying " << newgstaying
         << " doRemap " << doRemap << endl;
  }
  delete [] idx;
  delete [] sizes;
  delete [] adj;
  delete [] wgt;

  Teuchos::broadcast<int, int>(*comm, 0, 1, &doRemap);

  if (doRemap) {
    Teuchos::broadcast<int, partId_t>(*comm, 0, nGlobalParts, remap);
    for (size_t i = 0; i < len; i++) {
      parts[i] = remap[parts[i]];
    }
  }

  delete [] remap;  // TODO May want to keep for repartitioning as in Zoltan
}


}

#endif
