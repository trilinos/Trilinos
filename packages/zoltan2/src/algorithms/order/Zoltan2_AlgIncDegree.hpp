#ifndef _ZOLTAN2_INCDEGREE_HPP_
#define _ZOLTAN2_INCDEGREE_HPP_

#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <algorithm>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_IncDegree.hpp
//! \brief Order vertices by increasing degree

// Comparison function for sort.
bool comp(std::pair<lno_t,lno_t> a, std::pair<lno_t,lno_t> b)
{
  return (a.first < b.first);
}

namespace Zoltan2{

template <typename Adapter>
int AlgIncDegree(
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

  // Get local graph.
  ArrayView<const lno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<const scalar_t> wgts;
  model->getLocalEdgeList(edgeIds, offsets, wgts);

  // Store degrees together with index so we can sort.
  std::vector<std::pair<lno_t, lno_t> >  degrees(nVtx);
  for (lno_t i=0; i<nVtx; i++){
    degrees[i].first  = offsets[i+1] - offsets[i];
    degrees[i].second = i;
  }

  // Sort degrees.
  std::sort(degrees.begin(), degrees.end(), comp);

  // Copy permuted indices to perm.
  for (lno_t i=0; i<nVtx; i++){
    perm[i] = degrees[i].second;
  }

  return ierr;
}

}
#endif
