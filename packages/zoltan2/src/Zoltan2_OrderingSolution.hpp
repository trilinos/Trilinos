// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_OrderingSolution.hpp

    \brief The solution to a ordering problem.
*/


#ifndef _ZOLTAN2_ORDERINGSOLUTION_HPP_
#define _ZOLTAN2_ORDERINGSOLUTION_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Solution.hpp>

namespace Zoltan2 {

/*! Zoltan2::OrderingSolution
*/

template <typename Adapter>
  class OrderingSolution : public Solution<Adapter>
{
public:
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lid_t lid_t;
  typedef typename Adapter::lno_t lno_t;

  //////////////////////////////////////////////
  void setPermutation(
    size_t length,   // Length of arrays
    gid_t *gids,     // GIDs
    lid_t *lids,     // LIDs
    lno_t *perm     // perm[i] = k means k is the i'th element in the perm.
  )
  {
    HELLO;

    if (gids != NULL)
      gids_ = ArrayView<gid_t>(gids, length);
    else     // gids cannot be NULL
      gids_ = ArrayView<gid_t>(Teuchos::null);
      // throw std::logic_error("invalid gids");


    if (lids != NULL)
      lids_ = ArrayView<lid_t>(lids, length);
    else     // lids may be NULL
      lids_ = ArrayView<lid_t>(Teuchos::null);

    perm_ = ArrayView<lno_t>(perm, length);
  }

  //////////////////////////////////////////////
  void getPermutation(
    size_t *length,   // returned: Length of arrays
    gid_t **gids,     // returned: GIDs
    lid_t **lids,     // returned: LIDs
    lno_t **perm     // returned: Permutation
  )
  {
    *length = perm_.size();
    *gids   = gids_.getRawPtr();

    if (lids_.getRawPtr() != (lid_t*) Teuchos::null) *lids = lids_.getRawPtr();
    else                                             *lids = (lid_t*) NULL;

    *perm  = perm_.getRawPtr();
  }

protected:
  // Ordering solution consists of GIDs, LIDs, and permutation vector(s).
  size_t nParts_;
  ArrayView<gid_t>  gids_;
  ArrayView<lid_t>  lids_;
  ArrayView<lno_t> perm_;    // zero-based local permutation
  //ArrayView<size_t> invperm_; // inverse of permutation above
};

}

#endif
