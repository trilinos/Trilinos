#ifndef _ZOLTAN2_ALGRCM_HPP_
#define _ZOLTAN2_ALGRCM_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <queue>
//#define RCM


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
  const RCP<const Teuchos::Comm<int> > &comm
) 
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr= 0;

  HELLO;

  // TEST: return the identity permutation.
  const size_t nVtx = model->getLocalNumVertices();
  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
  for (lno_t i=0; i<nVtx; i++){
#ifdef RCM
    perm[i] = -1;
#else
    perm[i] = i;
#endif
  }

#ifdef RCM
  // This is the real RCM algorithm.
  // Get local graph.
  ArrayView<const lno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<const scalar_t> wgts;
  model->getLocalEdgeList(edgeIds, offsets, wgts);

  // TODO: Find pseudo-peripheral root vertex.
  lno_t root = 0;

  // Do BFS from root
  queue<lno_t> Q;
  lno_t count = nVtx-1; // start numbering from n-1 (Reverse CM)
  lno_t next = 0;

  while (count){ // Some vertex remains unlabelled

    // Label connected component starting at root
    Q.push(root);
    perm[root] = count--;

    while (Q.size()){
      // Get a vertex from the queue
      lno_t v = Q.front();
      Q.pop();

      // Add unmarked nbors to queue
      // TODO: Sort nbors by degree
      for (lno_t ptr = offsets[v]; ptr < offsets[v+1]; ++ptr){
        lno_t nbor = edgeIds[ptr];
        if (perm[nbor] != -1){
          perm[nbor] = count--; // Label as we push on Q
          Q.push(nbor);
        }
      }
    }

    // Find an unmarked vertex, use as new root
    while (perm[next] != -1) next++;
    root = next;
  }
#endif

  return ierr;
}

}
#endif
