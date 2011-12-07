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

template <typename gid_t, typename lno_t>
  class OrderingSolution : public Solution<gid_t, lno_t>
{
public:

  //////////////////////////////////////////////
  void setPermutation(
    size_t length,// Length of arrays
    gid_t *gids,  // User's IDs, same order as appeared in adapter "get" method
    lno_t *perm   // perm[i] = k means k is the i'th element in the perm.
  )
  {
    HELLO;

    if (gids != NULL)
      gids_ = ArrayView<gid_t>(gids, length);
    else     // gids cannot be NULL
      gids_ = ArrayView<gid_t>(Teuchos::null);
      // throw std::logic_error("invalid gids");

    perm_ = ArrayView<lno_t>(perm, length);
  }

  //////////////////////////////////////////////
  void getPermutation(
    size_t *length,  // returned: Length of arrays
    gid_t **gids,    // returned: GIDs
    lno_t **perm     // returned: Permutation
  )
  {
    *length = perm_.size();
    *gids   = gids_.getRawPtr();

    *perm  = perm_.getRawPtr();
  }

protected:
  // Ordering solution is GIDs listed in the same order in which they were
  // provided in the input adapter's "get" method, and permutation vector(s).
  size_t nParts_;
  ArrayView<gid_t>  gids_;
  ArrayView<lno_t> perm_;    // zero-based local permutation
  //ArrayView<size_t> invperm_; // inverse of permutation above
};

}

#endif
