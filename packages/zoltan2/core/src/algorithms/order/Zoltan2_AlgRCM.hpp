// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  const RCP<const typename Adapter::base_adapter_t> adapter;
  const RCP<Teuchos::ParameterList> pl;
  const RCP<const Teuchos::Comm<int> > comm;
  RCP<const Environment> env;
  modelFlag_t graphFlags;

  public:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;

  AlgRCM(
    const RCP<const typename Adapter::base_adapter_t> &adapter__,
    const RCP<Teuchos::ParameterList> &pl__,
    const RCP<const Teuchos::Comm<int> > &comm__,
    RCP<const Environment> &env__,
    const modelFlag_t &graphFlags__
  ) : adapter(adapter__), pl(pl__), comm(comm__), env(env__), graphFlags(graphFlags__)
  {
  }

  int globalOrder(const RCP<GlobalOrderingSolution<gno_t> > &/* solution */)
  {
    throw std::logic_error("AlgRCM does not yet support global ordering.");
  }

  int localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution)
  {
    int ierr= 0;

    HELLO;

    // Get local graph.
    ArrayView<const gno_t> edgeIds;
    ArrayView<const offset_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > wgts;

    const auto model = rcp(new GraphModel<Adapter>(adapter, env, comm, graphFlags));
    const size_t nVtx = model->getLocalNumVertices();
    model->getEdgeList(edgeIds, offsets, wgts);
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
    const ArrayRCP<lno_t> invPerm = solution->getPermutationRCP(true);
    const ArrayRCP<lno_t> tmpPerm(invPerm.size()); //temporary array used in reversing order

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
    gno_t root = 0;
    std::queue<gno_t> Q;
    size_t count = 0; // CM label, reversed later
    size_t next = 0;  // next unmarked vertex
    Teuchos::Array<std::pair<gno_t, offset_t> >  children; // children and their degrees

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
      tmpPerm[invPerm[root]] = root;

      while (Q.size()){
        // Get a vertex from the queue
        gno_t v = Q.front();
        Q.pop();
        //cout << "Debug: v= " << v << ", offsets[v] = " << offsets[v] << endl;

        // Add unmarked children to list of pairs, to be added to queue.
        children.resize(0);
        for (offset_t ptr = offsets[v]; ptr < offsets[v+1]; ++ptr){
          gno_t child = edgeIds[ptr];
          if (static_cast<Tpetra::global_size_t>(invPerm[child]) == INVALID){
            // Not visited yet; add child to list of pairs.
            std::pair<gno_t,offset_t> newchild;
            newchild.first = child;
            newchild.second = offsets[child+1] - offsets[child];
            children.push_back(newchild);
          }
        }
        // Sort children by increasing degree
        // TODO: If edge weights, sort children by decreasing weight,
        SortPairs<gno_t,offset_t> zort;
        zort.sort(children);

        typename Teuchos::Array<std::pair<gno_t,offset_t> >::iterator it = children.begin ();
        for ( ; it != children.end(); ++it){
          // Push children on the queue in sorted order.
          gno_t child = it->first;
          invPerm[child] = count++; // Label as we push on Q
          tmpPerm[invPerm[child]] = child;
          Q.push(child);
          //cout << "Debug: invPerm[" << child << "] = " << count << endl;
        }
      }
    }

    // Old row tmpPerm[i] is now the new row i.

    // Reverse labels for RCM.
    bool reverse = true; // TODO: Make parameter
    if (reverse) {
      lno_t temp;
      for (size_t i=0; i < nVtx/2; ++i) {
        // This effectively does the work of two loops:
        //    1) for (i=1; i< nVtx/2; ++i)
        //         swap of tmpPerm[i] and tmpPerm[nVtx-1-i]
        //    2) for (i=0; i < nVtx; ++i)
        //         invPerm[tmpPerm[i]] = i;
        temp = tmpPerm[i];
        invPerm[tmpPerm[nVtx-1-i]] = i;
        invPerm[temp] = nVtx-1-i;
      }

    }

    solution->setHaveInverse(true);
    return ierr;
  }

  private:
  // Find a smallest degree vertex in component containing v
  gno_t findSmallestDegree(
    gno_t v,
    lno_t nVtx,
    ArrayView<const gno_t> edgeIds,
    ArrayView<const offset_t> offsets)
  {
    std::queue<gno_t> Q;
    Teuchos::Array<bool> mark(nVtx);

    // Do BFS and compute smallest degree as we go
    offset_t smallestDegree = nVtx;
    gno_t smallestVertex = 0;

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
      offset_t deg = offsets[v+1] - offsets[v];
      if (deg < smallestDegree){
        smallestDegree = deg;
        smallestVertex = v;
      }
      // Add unmarked children to queue
      for (offset_t ptr = offsets[v]; ptr < offsets[v+1]; ++ptr){
        gno_t child = edgeIds[ptr];
        if (!mark[child]){
          mark[child] = true;
          Q.push(child);
        }
      }
    }
    return smallestVertex;
  }

  // Find a pseudoperipheral vertex in component containing v
  gno_t findPseudoPeripheral(
    gno_t v,
    lno_t nVtx,
    ArrayView<const gno_t> edgeIds,
    ArrayView<const offset_t> offsets)
  {
    std::queue<gno_t> Q;
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
        for (offset_t ptr = offsets[v]; ptr < offsets[v+1]; ++ptr){
          gno_t child = edgeIds[ptr];
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
