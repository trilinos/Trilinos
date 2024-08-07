// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGNATURAL_HPP_
#define _ZOLTAN2_ALGNATURAL_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgNatural.hpp
//! \brief Natural ordering == identity permutation.
//! \brief Mainly useful for testing "no ordering"

template <typename Adapter>
class AlgNatural : public Algorithm<Adapter>
{
  private:

  const RCP<const typename Adapter::base_adapter_t> adapter;
  const RCP<Teuchos::ParameterList> pl;
  const RCP<const Teuchos::Comm<int> > comm;
  const RCP<const Environment> env;

  public:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;

  AlgNatural(
    const RCP<const typename Adapter::base_adapter_t> &adapter__,
    const RCP<Teuchos::ParameterList> &pl__,
    const RCP<const Teuchos::Comm<int> > &comm__,
    const RCP<const Environment> &env__
  ) : adapter(adapter__), pl(pl__), comm(comm__), env(env__)
  {
  }

  int globalOrder(const RCP<GlobalOrderingSolution<gno_t> > &/* solution */) {
    throw std::logic_error("AlgNatural does not yet support global ordering.");
  }

  int localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution)
  {

    int ierr= 0;

    HELLO;

    // Local permutation only for now.

    // Set identity permutation.
    modelFlag_t modelFlags;
    const auto model = rcp(new IdentifierModel<Adapter>(adapter, env, comm, modelFlags));
    const size_t n = model->getLocalNumIdentifiers();
    lno_t *perm = solution->getPermutationView();
    if (perm){
      for (size_t i=0; i<n; i++){
        perm[i] = i;
      }
    }
    else
      // TODO: throw exception?
      ierr = -1;

    solution->setHavePerm(true);
    return ierr;

  }

};
}
#endif
