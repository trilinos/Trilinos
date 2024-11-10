// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGRANDOM_HPP_
#define _ZOLTAN2_ALGRANDOM_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


namespace Zoltan2{

////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgRandom.hpp
//! \brief Random ordering using the Fisher-Yates (Knuth) shuffle.
//! \brief TODO: Only local permutation, could add global option.

template <typename Adapter>
class AlgRandom : public Algorithm<Adapter>
{
  private:

  const RCP<const typename Adapter::base_adapter_t> adapter;
  const RCP<Teuchos::ParameterList> pl;
  const RCP<const Teuchos::Comm<int> > comm;
  const RCP<const Environment> env;

  public:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;

  AlgRandom(
    const RCP<const typename Adapter::base_adapter_t> &adapter__,
    const RCP<Teuchos::ParameterList> &pl__,
    const RCP<const Teuchos::Comm<int> > &comm__,
    const RCP<const Environment> &env__
  ) : adapter(adapter__), pl(pl__), comm(comm__), env(env__)
  {
  }

  int globalOrder(const RCP<GlobalOrderingSolution<gno_t> > &/* solution */)
  {
    throw std::logic_error("AlgRandom does not yet support global ordering.");
  }

  int localOrder(const RCP<LocalOrderingSolution<lno_t> > &solution)
  {

    int ierr= 0;

    HELLO;

    // This is the Fisher-Yates shuffle (aka Knuth shuffle).
    // References:
    //   R.A. Fisher, F. Yates, (1948) [1938], Statistical tables for biological, agricultural and medical research (3rd ed.).
    //   R. Durstenfeld, "Algorithm 235: Random permutation", CACM, vol. 7, 1964.
    //   D.E. Knuth, "The Art of Computer Programming", volume 2, 1969.

    // Start with the identity permutation.
    modelFlag_t modelFlags;
    const auto model = rcp(new IdentifierModel<Adapter>(adapter, env, comm, modelFlags));
    const size_t n = model->getLocalNumIdentifiers();
    lno_t *perm;
    perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
    if (perm){
      for (size_t i=0; i<n; i++){
        perm[i] = i;
      }
    }
    else
      // throw exception?
      ierr = -1;

    // Swap random pairs of indices in perm.
    lno_t j, temp;
    for (lno_t i=n-1; i>0; i--){
      // Choose j randomly in [0,i]
      j = rand() % (i+1);
      // Swap (perm[i], perm[j])
      temp = perm[i];
      perm[i] = perm[j];
      perm[j] = temp;
    }

    solution->setHavePerm(true);
    return ierr;

  }

};
}
#endif
