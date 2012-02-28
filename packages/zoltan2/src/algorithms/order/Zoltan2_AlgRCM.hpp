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
  const RCP<const Teuchos::Comm<int> > &comm
) 
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr= 0;

  HELLO;

  const size_t nVtx = model->getLocalNumVertices();
  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
  for (lno_t i=0; i<nVtx; i++){
    perm[i] = -1;
  }

  // Get local graph.
  ArrayView<const lno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<const scalar_t> wgts;
  model->getLocalEdgeList(edgeIds, offsets, wgts);
  //model->getLocalEdgeList(edgeIds, offsets);

  //cout << "Debug: Local graph from getLocalEdgeList" << endl;
  //cout << "edgeIds: " << edgeIds << endl;
  //cout << "offsets: " << offsets << endl;

  // TODO: Find pseudo-peripheral root vertex.
  lno_t root = 0;

  // Do BFS from root
  std::queue<lno_t> Q;
  lno_t count = 0; // CM label, reversed later
  lno_t next = 0;

  while (count < nVtx-1){ // Some vertex remains unlabelled

    // Label connected component starting at root
    Q.push(root);
    cout << "Debug: perm[" << root << "] = " << count << endl;
    perm[root] = count++;

    while (Q.size()){
      // Get a vertex from the queue
      lno_t v = Q.front();
      Q.pop();

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
    for (lno_t i=0; i < nVtx/2; ++i) {
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
