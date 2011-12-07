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
  // Constructor allocates memory for the solution.
  OrderingSolution(
    size_t perm_size, // TODO: Is this always equal to nlids ?
    size_t ngids
  )
  {
    HELLO;
    perm_size_ = perm_size;
    gids_   = ArrayRCP<gid_t>(ngids);
    perm_  = ArrayRCP<size_t>(perm_size_);
  }


  //////////////////////////////////////////////
  // Accessor functions, allowing algorithms to get ptrs to solution memory.
  // Algorithms can then load the memory.
  // Non-RCP versions are provided for applications to use.
  inline size_t getPermutationSize() {return perm_size_;}

  inline ArrayRCP<gid_t>  &getGidsRCP()  {return gids_;}
  inline ArrayRCP<lno_t> &getPermutationRCP() {return perm_;}

  inline ArrayRCP<gid_t>  &getGidsRCPConst()  const
  {
    return const_cast<ArrayRCP<gid_t>& > (gids_);
  }
  inline ArrayRCP<lno_t> &getPermutationRCPConst() const
  {
    return const_cast<ArrayRCP<lno_t>& > (perm_);
  }

  inline gid_t  *getGids(size_t *length)
  {
    *length = gids_.size();
    return gids_.getRawPtr();
  }
  inline lno_t *getPermutation(size_t *length)
  {
    *length = perm_.size();
    return perm_.getRawPtr();
  }

protected:
  // Ordering solution consists of GIDs, LIDs, and permutation vector(s).
  size_t perm_size_;
  ArrayView<gid_t>  gids_;
  ArrayView<lno_t> perm_;    // zero-based local permutation
  //ArrayView<size_t> invperm_; // inverse of permutation above
}

}

#endif
