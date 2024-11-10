// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_SORTEDDEGREE_HPP_
#define _ZOLTAN2_SORTEDDEGREE_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>
#include <Zoltan2_Sort.hpp>
#include <algorithm>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgSortedDegree.hpp
//! \brief Order vertices by sorted (increasing) degree.
//! \brief Sorting by decreasing degree is also possible.

namespace Zoltan2{

template <typename Adapter>
class AlgSortedDegree : public Algorithm<Adapter>
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

  AlgSortedDegree(
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
     throw std::logic_error(
       "AlgSortedDegree does not yet support global ordering.");
  }

  int localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution)
  {
    int ierr= 0;

    HELLO;

    lno_t *perm = solution->getPermutationView();
    if (perm==0){
      // Throw exception
      std::cerr << "perm is NULL" << std::endl;
      ierr = -1;
    }

    // Get local graph.
    const auto model = rcp(new GraphModel<Adapter>(adapter, env, comm, graphFlags));
    const size_t nVtx = model->getLocalNumVertices();
    ArrayView<const gno_t> edgeIds;
    ArrayView<const offset_t> offsets;
    ArrayView<StridedData<lno_t, scalar_t> > wgts;
    model->getEdgeList(edgeIds, offsets, wgts);


    // Store degrees together with index so we can sort.
    Teuchos::Array<std::pair<lno_t, offset_t> >  degrees(nVtx);
    for (lno_t i=0; i<(lno_t)nVtx; i++){
      degrees[i].first = i;
      degrees[i].second  = offsets[i+1] - offsets[i];
    }

    // Sort degrees.
    SortPairs<lno_t,offset_t> zort;
    bool inc = true; // TODO: Add user parameter
    zort.sort(degrees, inc);

    // Copy permuted indices to perm.
    for (lno_t i=0; i<(lno_t)nVtx; i++){
      perm[i] = degrees[i].first;
    }

    solution->setHavePerm(true);
    return ierr;
  }

};
}
#endif
