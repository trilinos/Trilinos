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

/*! \file Zoltan2_componentMetrics.hpp
    \brief Identify and compute the number of connected components in a processor's input
    Note that this routine works with respect to the MPI PROCESS, not with 
    respect to part numbers.  It works with the MPI Process' LOCAL graph; 
    statistics reported are for the local graph, not the global graph.
*/

#ifndef _ZOLTAN2_COMPONENTMETRICS_HPP_
#define _ZOLTAN2_COMPONENTMETRICS_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Teuchos_Comm.hpp>
#include <queue>

namespace Zoltan2
{

template <typename Adapter>
class perProcessorComponentMetrics{

public:
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::offset_t offset_t;
  typedef typename Adapter::scalar_t scalar_t;

  perProcessorComponentMetrics(const Adapter &ia,
                               const Teuchos::Comm<int> &comm);

#ifdef HAVE_ZOLTAN2_MPI
  // Wrap MPI_Comm as a Teuchos::Comm, then call Teuchos::Comm constructor
  // Uses delegating constructor feature of C++11.
  perProcessorComponentMetrics(const Adapter &ia, const MPI_Comm mpicomm) :
  perProcessorComponentMetrics(ia, 
                       Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpicomm)))
  {}
#endif

  ~perProcessorComponentMetrics() {}

  inline size_t getNumComponents() {return nComponent;}
  inline size_t getMaxComponentSize() {return maxComponentSize;}
  inline size_t getMinComponentSize() {return minComponentSize;}
  inline double getAvgComponentSize() {return avgComponentSize;}

private:
  size_t nComponent;    // number of components
  size_t maxComponentSize;  // size of largest component
  size_t minComponentSize;  // size of smalled component
  double avgComponentSize;  // average component size

  inline void markAndEnqueue(std::queue<gno_t> &q, bool *mark,
                             size_t &nUnmarkedVtx, size_t &cSize, gno_t vtx) {
    // insert vtx into the queue
    q.push(vtx);

    // mark vtx
    nUnmarkedVtx--;
    mark[vtx] = true;

    // increment component size
    cSize++;
  }
};

///////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
perProcessorComponentMetrics<Adapter>::perProcessorComponentMetrics(
  const Adapter &ia, const Teuchos::Comm<int> &comm) :
  nComponent(0), maxComponentSize(0), minComponentSize(0),
  avgComponentSize(0.)
{
  // build local graph model from input adapter
  std::bitset<NUM_MODEL_FLAGS> graphFlags;
  graphFlags.set(REMOVE_SELF_EDGES);
  graphFlags.set(BUILD_LOCAL_GRAPH);  // Local graph; 
                                      // all vertex numbering is 0 to nVtx-1

  Teuchos::RCP<const Teuchos::Comm<int> > tcomm = rcp(&comm, false);
  Teuchos::RCP<const Zoltan2::Environment> env = 
               rcp(new Zoltan2::Environment(tcomm));

  typedef typename Adapter::base_adapter_t base_adapter_t;
  Teuchos::RCP<const base_adapter_t> ria = rcp(&ia, false);
  Zoltan2::GraphModel<base_adapter_t> graph(ria, env, tcomm, graphFlags);

  // get graph from model
  const size_t nVtx = graph.getLocalNumVertices();
  ArrayView<const gno_t> adj;
  ArrayView<const offset_t> offset;
  ArrayView<StridedData<lno_t, scalar_t> > wgts;  // unused
  graph.getEdgeList(adj, offset, wgts);

  // do a simple BFS on the graph; 
  // can replace later with KokkosKernels or other TPL
  size_t nUnmarkedVtx = nVtx;
  bool *mark = new bool[nUnmarkedVtx];
  for (size_t i = 0; i < nUnmarkedVtx; i++) mark[i] = false;
  size_t startVtx = 0;

  std::queue<gno_t> q;

  // until all vertices are marked...
  while (nUnmarkedVtx > 0) {

    // Find the next component
    size_t cSize = 0;
    nComponent++;

    // Find an unmarked vertex; put it in the queue
    while (mark[startVtx]) startVtx++;
    markAndEnqueue(q, mark, nUnmarkedVtx, cSize, startVtx);

    while (!q.empty()) {
      gno_t vtx = q.front();
      q.pop();

      // Add neighbors of vtx to queue.
      for (offset_t j = offset[vtx]; j < offset[vtx+1]; j++) {
        if (!mark[adj[j]]) {
          markAndEnqueue(q, mark, nUnmarkedVtx, cSize, adj[j]);
        }
      }
    }

    // update stats
    if (nComponent == 1) {
      maxComponentSize = cSize;
      minComponentSize = cSize;
    }
    else {
      if (cSize > maxComponentSize) maxComponentSize = cSize;
      if (cSize < minComponentSize) minComponentSize = cSize;
    }
  }

  // update stats
  if (nComponent) avgComponentSize = double(nVtx) / double(nComponent);

  delete [] mark;
}


}                   // namespace Zoltan2
#endif
