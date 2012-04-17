#ifndef _ZOLTAN2_ALGRANDOM_HPP_
#define _ZOLTAN2_ALGRANDOM_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_OrderingSolution.hpp>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_AlgRandom.hpp
//! \brief Random ordering using the Knuth shuffle.
//! \brief TODO: Only local permutation, should add global option.


namespace Zoltan2{

template <typename Adapter>
int AlgRandom(
  const RCP<IdentifierModel<Adapter> > &model, 
  const RCP<OrderingSolution<typename Adapter::gid_t,
                             typename Adapter::lno_t> > &solution,
  const RCP<Teuchos::ParameterList> &pl,
  const RCP<Teuchos::Comm<int> > &comm
) 
{
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::scalar_t scalar_t;

  int ierr= 0;

  HELLO;

  // This is the classic Knuth shuffle, also known as Fisher-Yates shuffle.
  // References:
  //   D.E. Knuth, "The Art of Computer Programming", volume 2, 1969.
  //   R. Durstenfeld, "Algorithm 235: Random permutation", CACM, vol. 7, 1964.

  // Start with the identity permutation.
  const size_t n = model->getLocalNumIdentifiers();
  lno_t *perm;
  perm = (lno_t *) (solution->getPermutationRCP().getRawPtr());
  if (perm){
    for (lno_t i=0; i<n; i++){
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

  return ierr;
}

}
#endif
