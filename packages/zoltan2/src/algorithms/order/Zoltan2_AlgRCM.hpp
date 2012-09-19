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
#ifndef _ZOLTAN2_ALGRCM_HPP_
#define _ZOLTAN2_ALGRCM_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <queue>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgRCM.hpp
//! \brief RCM ordering of a graph (serial, local graph only)


namespace Zoltan2{

template <typename Adapter>
int AlgRCM(
  const RCP<GraphModel<Adapter> > &model, 
  const RCP<OrderingSolution<typename Adapter::gid_t,
                             typename Adapter::lno_t> > &solution,
  const RCP<Teuchos::ParameterList> &pl,
  const RCP<Teuchos::Comm<int> > &comm
) 
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr= 0;

  HELLO;

  // Check size of communicator: serial only.
  // TODO: Remove this test when RCM works on local graph.
  if (comm->getSize() > 1){
    throw std::runtime_error("RCM currently only works in serial.");
  }

  const size_t nVtx = model->getLocalNumVertices();
  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
  for (size_t i=0; i<nVtx; i++){
    perm[i] = -1;
  }

  // Get local graph.
  ArrayView<const lno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<StridedData<lno_t, scalar_t> > wgts;

  //model->getLocalEdgeList(edgeIds, offsets, wgts); // BUGGY!
  // Use global graph for now. This only works in serial!
  ArrayView<const int> procIds;
  model->getEdgeList( edgeIds, procIds, offsets, wgts);

  //cout << "Debug: Local graph from getLocalEdgeList" << endl;
  //cout << "edgeIds: " << edgeIds << endl;
  //cout << "offsets: " << offsets << endl;

  // TODO: Find min-degree (or pseudo-peripheral) root vertex.
  lno_t root = 0;

  // Do BFS from root
  std::queue<lno_t> Q;
  size_t count = 0; // CM label, reversed later
  size_t next = 0;

  while (count < nVtx-1){ // Some vertex remains unlabelled

    // Label connected component starting at root
    Q.push(root);
    //cout << "Debug: perm[" << root << "] = " << count << endl;
    perm[root] = count++;

    while (Q.size()){
      // Get a vertex from the queue
      lno_t v = Q.front();
      Q.pop();
      //cout << "Debug: v= " << v << ", offsets[v] = " << offsets[v] << endl;

      // Add unmarked nbors to queue
      // TODO: If edge weights, sort nbors by decreasing weight,
      // TODO: Else, sort nbors by increasing degree
      for (lno_t ptr = offsets[v]; ptr < offsets[v+1]; ++ptr){
        lno_t nbor = edgeIds[ptr];
        if (perm[nbor] == -1){
          //cout << "Debug: perm[" << nbor << "] = " << count << endl;
          perm[nbor] = count++; // Label as we push on Q
          Q.push(nbor);
        }
      }
    }

    // Find an unmarked vertex, use as new root
    while ((next < nVtx) && (perm[next] != -1)) next++;
    root = next;
  }

  // Reverse labels for RCM
  bool reverse = true; // TODO: Make parameter
  if (reverse) {
    lno_t temp;
    for (size_t i=0; i < nVtx/2; ++i) {
      // Swap (perm[i], perm[nVtx-i])
      temp = perm[i];
      perm[i] = perm[nVtx-1-i];
      perm[nVtx-1-i] = temp;
    }
  }

  return ierr;
}

}
#endif
