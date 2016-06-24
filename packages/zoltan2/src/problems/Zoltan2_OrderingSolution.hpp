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

/*! \brief The class containing ordering solutions.

    Template parameters:
    \li \c gno_t    data type for application global Ids
    \li \c lno_t    data type for local indices and local counts

The ordering solution always contains the permutation and the inverse permutation. These should be accessed through the accessor methods defined in this class, such as getPermutation(). Some ordering algorithms may compute and store other information. Currently, only serial ordering of the local data is supported.

In Zoltan2, perm[i]=j means index i in the reordered vector/matrix corresponds to index j in the old ordering. In Matlab notation, A(perm,perm) is the reordered matrix. This is consistent with SuiteSparse (AMD) and several other ordering packages. Unfortunately, this notation may conflict with a few other packages (such as Ifpack2). 

*/

template <typename lno_t, typename gno_t>
  class OrderingSolution : public Solution
{
public:

  /*! \brief Constructor allocates memory for the solution.
   */
  OrderingSolution(
    size_t perm_size // This should be equal to nlids
  ) : separatorColBlocks_(0)
  {
    HELLO;
    perm_size_  = perm_size;
    //gids_       = ArrayRCP<gno_t>(perm_size_);
    perm_       = ArrayRCP<lno_t>(perm_size_);
    invperm_    = ArrayRCP<lno_t>(perm_size_);
    separatorRange_   = ArrayRCP<lno_t>(perm_size_+1);
    separatorTree_    = ArrayRCP<lno_t>(perm_size_);

    havePerm_ = false;
    haveInverse_ = false;
    haveSeparatorRange_ = false;
    haveSeparatorTree_ = false;
  }

  /*! \brief Do we have the direct permutation?
   */
  bool havePerm() const
  {
    return havePerm_; 
  }

  /*! \brief Set havePerm (intended for ordering algorithms only)
   */
  void setHavePerm(bool status)
  {
    havePerm_ = status; 
  }


  /*! \brief Do we have the inverse permutation?
   */
  bool haveInverse() const
  {
    return haveInverse_; 
  }

  /*! \brief Set haveInverse (intended for ordering algorithms only)
   */
  void setHaveInverse(bool status)
  {
    haveInverse_ = status; 
  }
 
  /*! \brief set all separator flags.
   */
 void setHaveSeparator(bool status) {
    this->setHavePerm(status);
    this->setHaveInverse(status);
    this->setHaveSeparatorRange(status);
    this->setHaveSeparatorTree(status);
 } 
  /*! \brief Do we have the seperator range?
   */
  bool haveSeparatorRange() const
  {
    return haveSeparatorRange_; 
  }

  /*! \brief Set haveSeparatorRange (intended for ordering algorithms only)
   */
  void setHaveSeparatorRange(bool status)
  {
    haveSeparatorRange_ = status; 
  }
  
  /*! \brief Do we have the seperator tree?
   */
  bool haveSeparatorTree() const
  {
    return haveSeparatorTree_; 
  }
  
  /*! \brief Do we have the seperators?
   */
  bool haveSeparators() const
  {
    return haveSeparatorRange() && haveSeparatorTree(); 
  }

  /*! \brief Set haveSeparatorTree (intended for ordering algorithms only)
   */
  void setHaveSeparatorTree(bool status)
  {
    haveSeparatorTree_ = status; 
  }

  /*! \brief Compute direct permutation from inverse.
   */
  void computePerm()
  {
    if (haveInverse_) {
      for(size_t i=0; i<perm_size_; i++) {
        perm_[invperm_[i]] = i;
      }
      havePerm_ = true;
    }
    else {
      // TODO: throw exception
      std::cerr << "No inverse!" << std::endl;
    }
  }

  /*! \brief Compute inverse permutation.
   */
  void computeInverse()
  {
    if (havePerm_) {
      for(size_t i=0; i<perm_size_; i++) {
        invperm_[perm_[i]] = i;
      }
      havePerm_ = true;
    }
    else {
      // TODO: throw exception
      std::cerr << "No perm!" << std::endl;
    }
  }

  /*! \brief Set number of separator column blocks.
   */
  inline void setNumSeparatorBlocks(lno_t nblks) {separatorColBlocks_ = nblks;}

  //////////////////////////////////////////////
  // Accessor functions, allowing algorithms to get ptrs to solution memory.
  // Algorithms can then load the memory.
  // Non-RCP versions are provided for applications to use.

  /*! \brief Get (local) size of permutation.
   */
  inline size_t getPermutationSize() const {return perm_size_;}
  
  /*! \brief Get number of separator column blocks.
   */
  inline lno_t getNumSeparatorBlocks() const {return separatorColBlocks_;}

  /*! \brief Get (local) permuted GIDs by RCP.
   */
 // inline ArrayRCP<lno_t>  &getGidsRCP()  {return gids_;}

  /*! \brief Get (local) permutation by RCP.
   *  If inverse = true, return inverse permutation.
   *  By default, perm[i] is where new index i can be found in the old ordering.
   *  When inverse==true, perm[i] is where old index i can be found in the new ordering.
   */
  inline const ArrayRCP<lno_t> &getPermutationRCP(bool inverse=false) const
  {
    if (inverse)
      return invperm_;
    else
      return perm_;
  }
  
  /*! \brief return vertex separator variables by reference.
    */
  bool getVertexSeparator (lno_t &numBlocks,
                           lno_t *range,
                           lno_t *tree) const {

    if (this->haveSeparators()) {
      numBlocks = this->getNumSeparatorBlocks();
      range = this->getSeparatorRangeView();
      tree = this->getSeparatorTreeView();
      return true;
    }

    return false;
  }

  /*! \brief Get (local) seperator range by RCP.
   */
  inline const ArrayRCP<lno_t> &getSeparatorRangeRCP() const
  {
    return separatorRange_;
  }
  
  /*! \brief Get (local) seperator tree by RCP.
   */
  inline const ArrayRCP<lno_t> &getSeparatorTreeRCP() const
  {
    return separatorTree_;
  }

  /*! \brief Get (local) permuted GIDs by const RCP.
   */
 // inline ArrayRCP<lno_t>  &getGidsRCPConst()  const
 // {
 //   return const_cast<ArrayRCP<lno_t>& > (gids_);
 // }

  /*! \brief Get (local) permutation by const RCP.
   *  If inverse = true, return inverse permutation.
   *  By default, perm[i] is where new index i can be found in the old ordering.
   *  When inverse==true, perm[i] is where old index i can be found in the new ordering.
   */
  inline ArrayRCP<lno_t> &getPermutationRCPConst(bool inverse=false) const
  {
    if (inverse)
      return const_cast<ArrayRCP<lno_t>& > (invperm_);
    else
      return const_cast<ArrayRCP<lno_t>& > (perm_);
  }
  
  /*! \brief Get (local) seperator range by const RCP.
   */
  inline ArrayRCP<lno_t> &getSeparatorRangeRCPConst() const
  {
    return const_cast<ArrayRCP<lno_t> & > (separatorRange_);
  }
  
  /*! \brief Get (local) seperator tree by const RCP.
   */
  inline ArrayRCP<lno_t> &getSeparatorTreeRCPConst() const
  {
    return const_cast<ArrayRCP<lno_t> & > (separatorTree_);
  }

  /*! \brief Get pointer to (local) permutation.
   *  If inverse = true, return inverse permutation.
   *  By default, perm[i] is where new index i can be found in the old ordering.
   *  When inverse==true, perm[i] is where old index i can be found 
   *  in the new ordering.
   */
  inline lno_t *getPermutationView(bool inverse = false) const
  {
    if (perm_size_) {
      if (inverse)
        return invperm_.getRawPtr();
      else
        return perm_.getRawPtr();
    }
    else
      return NULL;
  }
  
  /*! \brief Get pointer to (local) separator range.
   */
  inline lno_t *getSeparatorRangeView() const
  {
    // Here, don't need to check perm_size_ before calling getRawPtr.
    // separatorRange_ always has some length, since it is allocated larger
    // than other arrays.
    return separatorRange_.getRawPtr();
  }

  /*! \brief Get pointer to (local) separator tree.
   */
  inline lno_t *getSeparatorTreeView() const
  {
    if (perm_size_)
      return separatorTree_.getRawPtr();
    else
      return NULL;
  }
  
  /*! \brief Get reference to (local) separator column block.
   */
  inline lno_t &NumSeparatorBlocks()
  {
    return separatorColBlocks_; 
  }

protected:
  // Ordering solution consists of permutation vector(s).
  // Either perm or invperm should be computed by the algorithm.
  size_t perm_size_;
  // For now, assume permutations are local. Revisit later (e.g., for Scotch)
  bool havePerm_;           // has perm_ been computed yet?
  bool haveInverse_;        // has invperm_ been computed yet?
  bool haveSeparatorRange_;          // has sepRange_ been computed yet?
  bool haveSeparatorTree_;        // has sepTree_ been computed yet?
  ArrayRCP<lno_t> perm_;    // zero-based local permutation
  ArrayRCP<lno_t> invperm_; // inverse of permutation above
  ArrayRCP<lno_t> separatorRange_;   // range iterator for separator tree
  ArrayRCP<lno_t> separatorTree_; // separator tree
  lno_t separatorColBlocks_;        // number of column blocks in separator
};

}

#endif
