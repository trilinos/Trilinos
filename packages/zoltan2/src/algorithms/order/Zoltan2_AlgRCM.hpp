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

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <Zoltan2_Sort.hpp>
#include <queue>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgRCM.hpp
//! \brief RCM ordering of a graph (serial, local graph only)


namespace Zoltan2{

template <typename Adapter>
class AlgRCM : public Algorithm<Adapter>
{
  private:

  const RCP<GraphModel<Adapter> > model;
  const RCP<Teuchos::ParameterList> &pl;
  const RCP<Teuchos::Comm<int> > &comm;

  public:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::zgid_t zgid_t;
  typedef typename Adapter::scalar_t scalar_t;

  AlgRCM(
    const RCP<GraphModel<Adapter> > &model__,
    const RCP<Teuchos::ParameterList> &pl__,
    const RCP<Teuchos::Comm<int> > &comm__
  ) : model(model__), pl(pl__), comm(comm__)
  {
  }

  int order(const RCP<OrderingSolution<zgid_t, lno_t> > &solution)
  {
    int ierr= 0;

    HELLO;
  
    // Get local graph.
    ArrayView<const lno_t> edgeIds;
    ArrayView<const lno_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > wgts;
  
    const size_t nVtx = model->getLocalNumVertices();
    model->getLocalEdgeList(edgeIds, offsets, wgts); 
    const int numWeightsPerEdge = model->getNumWeightsPerEdge();
    if (numWeightsPerEdge > 1){
      throw std::runtime_error("Multiple weights not supported.");
    }
  
#if 0
    // Debug
    cout << "Debug: Local graph from getLocalEdgeList" << endl;
    cout << "rank " << comm->getRank() << ": nVtx= " << nVtx << endl;
    cout << "rank " << comm->getRank() << ": edgeIds: " << edgeIds << endl;
    cout << "rank " << comm->getRank() << ": offsets: " << offsets << endl;
#endif
  
    // RCM constructs invPerm, not perm
    ArrayRCP<lno_t> invPerm = solution->getPermutationRCP(true);
  
    // Check if there are actually edges to reorder.
    // If there are not, then just use the natural ordering.
    if (offsets[nVtx] == 0) {
      for (size_t i = 0; i < nVtx; ++i) {
        invPerm[i] = i;
      }
      solution->setHaveInverse(true);
      return 0;
    }
  
    // Set the label of each vertex to invalid.
    Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    for (size_t i = 0; i < nVtx; ++i) {
      invPerm[i] = INVALID;
    }
  
    // Loop over all connected components.
    // Do BFS within each component.
    lno_t root = 0;
    std::queue<lno_t> Q;
    size_t count = 0; // CM label, reversed later
    size_t next = 0;  // next unmarked vertex
    Teuchos::Array<std::pair<lno_t, size_t> >  children; // children and their degrees
  
    while (count < nVtx) {
  
      // Find suitable root vertex for this component.
      // First find an unmarked vertex, use to find root in next component.
      while ((next < nVtx) && (static_cast<Tpetra::global_size_t>(invPerm[next]) != INVALID)) next++;

      // Select root method. Pseudoperipheral usually gives the best
      // ordering, but the user may choose a faster method.
      std::string root_method = pl->get("root_method", "pseudoperipheral");
      if (root_method == std::string("first"))
        root = next;
      else if (root_method == std::string("smallest_degree"))
        root = findSmallestDegree(next, nVtx, edgeIds, offsets);
      else if (root_method == std::string("pseudoperipheral"))
        root = findPseudoPeripheral(next, nVtx, edgeIds, offsets);
      else {
        // This should never happen if pl was validated.
        throw std::runtime_error("invalid root_method");
      }

      // Label connected component starting at root
      Q.push(root);
      //cout << "Debug: invPerm[" << root << "] = " << count << endl;
      invPerm[root] = count++;
  
      while (Q.size()){
        // Get a vertex from the queue
        lno_t v = Q.front();
        Q.pop();
        //cout << "Debug: v= " << v << ", offsets[v] = " << offsets[v] << endl;
  
        // Add unmarked children to list of pairs, to be added to queue.
        children.resize(0);
        for (lno_t ptr = offsets[v]; ptr < offsets[v+1]; ++ptr){
          lno_t child = edgeIds[ptr];
          if (static_cast<Tpetra::global_size_t>(invPerm[child]) == INVALID){
            // Not visited yet; add child to list of pairs.
            std::pair<lno_t,size_t> newchild;
            newchild.first = child;
            newchild.second = offsets[child+1] - offsets[child];
            children.push_back(newchild); 
          }
        }
        // Sort children by increasing degree
        // TODO: If edge weights, sort children by decreasing weight,
        SortPairs<lno_t,size_t> zort;
        zort.sort(children);

        typename Teuchos::Array<std::pair<lno_t,size_t> >::iterator it = children.begin ();
        for ( ; it != children.end(); ++it){
          // Push children on the queue in sorted order.
          lno_t child = it->first;
          invPerm[child] = count++; // Label as we push on Q
          Q.push(child);
          //cout << "Debug: invPerm[" << child << "] = " << count << endl;
        }
      }
    }
  
    // Reverse labels for RCM
    bool reverse = true; // TODO: Make parameter
    if (reverse) {
      lno_t temp;
      for (size_t i=0; i < nVtx/2; ++i) {
        // Swap (invPerm[i], invPerm[nVtx-i])
        temp = invPerm[i];
        invPerm[i] = invPerm[nVtx-1-i];
        invPerm[nVtx-1-i] = temp;
      }
    }
  
    solution->setHaveInverse(true);
    return ierr;
  }

  private:
  // Find a smallest degree vertex in component containing v
  lno_t findSmallestDegree(
    lno_t v,
    lno_t nVtx,
    ArrayView<const lno_t> edgeIds,
    ArrayView<const lno_t> offsets)
  {
    std::queue<lno_t> Q;
    Teuchos::Array<bool> mark(nVtx);

    // Do BFS and compute smallest degree as we go
    lno_t smallestDegree = nVtx;
    lno_t smallestVertex = 0;

    // Clear mark array - nothing marked yet
    for (int i=0; i<nVtx; i++)
      mark[i] = false;

    // Start from v
    Q.push(v);
    while (Q.size()){
      // Get first vertex from the queue
      v = Q.front();
      Q.pop();
      // Check degree of v
      lno_t deg = offsets[v+1] - offsets[v];
      if (deg < smallestDegree){
        smallestDegree = deg;
        smallestVertex = v;
      }
      // Add unmarked children to queue
      for (lno_t ptr = offsets[v]; ptr < offsets[v+1]; ++ptr){
        lno_t child = edgeIds[ptr];
        if (!mark[child]){
          mark[child] = true; 
          Q.push(child);
        }
      }
    }
    return smallestVertex;
  }

  // Find a pseudoperipheral vertex in component containing v
  lno_t findPseudoPeripheral(
    lno_t v,
    lno_t nVtx,
    ArrayView<const lno_t> edgeIds,
    ArrayView<const lno_t> offsets)
  {
    std::queue<lno_t> Q;
    Teuchos::Array<bool> mark(nVtx);

    // Do BFS a couple times, pick vertex last visited (furthest away)
    const int numBFS = 2;
    for (int bfs=0; bfs<numBFS; bfs++){
      // Clear mark array - nothing marked yet
      for (int i=0; i<nVtx; i++)
        mark[i] = false;
      // Start from v
      Q.push(v);
      while (Q.size()){
        // Get first vertex from the queue
        v = Q.front();
        Q.pop();
        // Add unmarked children to queue
        for (lno_t ptr = offsets[v]; ptr < offsets[v+1]; ++ptr){
          lno_t child = edgeIds[ptr];
          if (!mark[child]){
            mark[child] = true; 
            Q.push(child);
          }
        }
      }
    }
    return v;
  }
  
};
}
#endif
