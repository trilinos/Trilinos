#ifndef _ZOLTAN2_ALGSPECTRAL_HPP_
#define _ZOLTAN2_ALGSPECTRAL_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgSpectral.hpp
//! \brief Spectral ordering of a graph (local or global).
//! \brief Sorts the Fiedler vector of the graph Laplacian.


namespace Zoltan2{

template <typename Adapter>
int AlgSpectral(
  const RCP<GraphModel<Adapter> > &model, 
  //const RCP<Adapter> &matrixadapter, // Hack: Use matrix adapter directly
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

// TODO: Check params to do local or global ordering.
  bool localOrder = true;

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

  // Form the graph Laplacian: L = D-A
  // Create a new Tpetra matrix, but use views of existing data when possible.
  // TODO

  // TODO: Find smallest eigenvalues using Anasazi

  // TODO: Sort Fiedler vector.

  return ierr;
}

}
#endif
