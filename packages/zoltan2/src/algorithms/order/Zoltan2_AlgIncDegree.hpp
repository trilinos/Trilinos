#ifndef _ZOLTAN2_INCDEGREE_HPP_
#define _ZOLTAN2_INCDEGREE_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <algorithms>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_IncDegreeM.hpp
//! \brief Order vertices by increasing degree


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

  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
  if (perm){
    for (lno_t i=0; i<nVtx; i++){
      perm[i] = i;
    }
  }

  // Get local graph.
  ArrayView<const lno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<const scalar_t> wgts;
  model->getLocalEdgeList(edgeIds, offsets, wgts);

  // TODO: Compute degrees from offsets
  // TODO: Sort degrees.

  return ierr;
}

}
#endif
