// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_OrderingSolution.hpp
    \brief Defines the OrderingSolution class.
*/

#ifndef _ZOLTAN2_ORDERINGSOLUTION_HPP_
#define _ZOLTAN2_ORDERINGSOLUTION_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Solution.hpp>

namespace Zoltan2 {

/*! \brief The class containing ordering solutions and metrics.

    Template parameters:
    \li \c gid_t    data type for application global Ids
    \li \c lno_t    data type for local indices and local counts

   \todo documentation
*/

template <typename gid_t, typename lno_t>
  class OrderingSolution : public Solution
{
public:

  /*! \brief Constructor allocates memory for the solution.
   */
  OrderingSolution(
    size_t perm_size, // TODO: Is this always equal to nlids ?
    size_t ngids
  )
  {
    HELLO;
    perm_size_ = perm_size;
    gids_   = ArrayRCP<gid_t>(ngids);
    perm_  = ArrayRCP<lno_t>(perm_size_);
  }


  //////////////////////////////////////////////
  // Accessor functions, allowing algorithms to get ptrs to solution memory.
  // Algorithms can then load the memory.
  // Non-RCP versions are provided for applications to use.

  /*! \brief TODO.
   */
  inline size_t getPermutationSize() {return perm_size_;}

  /*! \brief TODO.
   */
  inline ArrayRCP<gid_t>  &getGidsRCP()  {return gids_;}

  /*! \brief TODO.
   */
  inline ArrayRCP<lno_t> &getPermutationRCP() {return perm_;}

  /*! \brief TODO.
   */
  inline ArrayRCP<gid_t>  &getGidsRCPConst()  const
  {
    return const_cast<ArrayRCP<gid_t>& > (gids_);
  }

  /*! \brief TODO.
   */
  inline ArrayRCP<lno_t> &getPermutationRCPConst() const
  {
    return const_cast<ArrayRCP<lno_t>& > (perm_);
  }

  /*! \brief TODO.
   */
  inline gid_t  *getGids(size_t *length)
  {
    *length = gids_.size();
    return gids_.getRawPtr();
  }

  /*! \brief TODO.
   */
  inline lno_t *getPermutation(size_t *length)
  {
    *length = perm_.size();
    return perm_.getRawPtr();
  }

protected:
  // Ordering solution consists of GIDs, LIDs, and permutation vector(s).
  size_t perm_size_;
  ArrayRCP<gid_t>  gids_;
  ArrayRCP<lno_t> perm_;    // zero-based local permutation
  //ArrayRCP<size_t> invperm_; // inverse of permutation above
};

}

#endif
