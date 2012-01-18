#ifndef _ZOLTAN2_ALGNATURAL_HPP_
#define _ZOLTAN2_ALGNATURAL_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgNatural.hpp
//! \brief Natural ordering == identity permutation.


namespace Zoltan2{

template <typename Adapter>
int AlgNatural(
  const RCP<IdentifierModel<Adapter> > &model, 
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

  // Local permutation only for now.

  // Set identity permutation.
  const size_t n = model->getLocalNumIdentifiers();
  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
  if (perm){
    for (lno_t i=0; i<n; i++){
      perm[i] = i;
    }
  }
  else
    // TODO: throw exception?
    ierr = -1;

  return ierr;
}

}
#endif
