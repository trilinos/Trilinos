// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    \li \c index_t    index type index_t will be lno_t (local) or gno_t (global)

The ordering solution always contains the permutation and the inverse permutation. These should be accessed through the accessor methods defined in this class, such as getPermutation(). Some ordering algorithms may compute and store other information. Currently, only serial ordering of the local data is supported.

In Zoltan2, perm[i]=j means index i in the reordered vector/matrix corresponds to index j in the old ordering. In Matlab notation, A(perm,perm) is the reordered matrix. This is consistent with SuiteSparse (AMD) and several other ordering packages. Unfortunately, this notation may conflict with a few other packages (such as Ifpack2). 

*/

template <typename index_t>
class OrderingSolution : public Solution
{
public:

  /*! \brief Constructor allocates memory for the solution.
   */
  OrderingSolution(index_t perm_size) : separatorColBlocks_(0)
  {
    HELLO;
    perm_size_        = perm_size;
    perm_             = ArrayRCP<index_t>(perm_size_);
    invperm_          = ArrayRCP<index_t>(perm_size_);
    separatorRange_   = ArrayRCP<index_t>(perm_size_+1);
    separatorTree_    = ArrayRCP<index_t>(perm_size_);

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
  /*! \brief Do we have the separator range?
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
  
  /*! \brief Do we have the separator tree?
   */
  bool haveSeparatorTree() const
  {
    return haveSeparatorTree_; 
  }
  
  /*! \brief Do we have the separators?
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
      haveInverse_ = true;
    }
    else {
      // TODO: throw exception
      std::cerr << "No perm!" << std::endl;
    }
  }

  /*! \brief Set number of separator column blocks.
   */
  inline void setNumSeparatorBlocks(index_t nblks) {separatorColBlocks_ = nblks;}

  //////////////////////////////////////////////
  // Accessor functions, allowing algorithms to get ptrs to solution memory.
  // Algorithms can then load the memory.
  // Non-RCP versions are provided for applications to use.

  /*! \brief Get (local) size of permutation.
   */
  inline size_t getPermutationSize() const {return perm_size_;}
  
  /*! \brief Get number of separator column blocks.
   */
  inline index_t getNumSeparatorBlocks() const {return separatorColBlocks_;}

  /*! \brief Get (local) permuted GIDs by RCP.
   */
 // inline ArrayRCP<index_t>  &getGidsRCP()  {return gids_;}

  /*! \brief Get (local) permutation by RCP.
   *  If inverse = true, return inverse permutation.
   *  By default, perm[i] is where new index i can be found in the old ordering.
   *  When inverse==true, perm[i] is where old index i can be found in the new ordering.
   */
  inline const ArrayRCP<index_t> &getPermutationRCP(bool inverse=false) const
  {
    if (inverse)
      return invperm_;
    else
      return perm_;
  }
  
  /*! \brief return vertex separator variables by reference.
    */
  bool getVertexSeparator (index_t &numBlocks,
                           index_t *range,
                           index_t *tree) const {

    if (this->haveSeparators()) {
      numBlocks = this->getNumSeparatorBlocks();
      range = this->getSeparatorRangeView();
      tree = this->getSeparatorTreeView();
      return true;
    }

    return false;
  }

  /*! \brief Get (local) separator range by RCP.
   */
  inline const ArrayRCP<index_t> &getSeparatorRangeRCP() const
  {
    return separatorRange_;
  }
  
  /*! \brief Get (local) separator tree by RCP.
   */
  inline const ArrayRCP<index_t> &getSeparatorTreeRCP() const
  {
    return separatorTree_;
  }

  /*! \brief Get (local) permuted GIDs by const RCP.
   */
 // inline ArrayRCP<index_t>  &getGidsRCPConst()  const
 // {
 //   return const_cast<ArrayRCP<index_t>& > (gids_);
 // }

  /*! \brief Get (local) permutation by const RCP.
   *  If inverse = true, return inverse permutation.
   *  By default, perm[i] is where new index i can be found in the old ordering.
   *  When inverse==true, perm[i] is where old index i can be found in the new ordering.
   */
  inline ArrayRCP<index_t> &getPermutationRCPConst(bool inverse=false) const
  {
    if (inverse)
      return const_cast<ArrayRCP<index_t>& > (invperm_);
    else
      return const_cast<ArrayRCP<index_t>& > (perm_);
  }
  
  /*! \brief Get separator range by const RCP.
   */
  inline ArrayRCP<index_t> &getSeparatorRangeRCPConst() const
  {
    return const_cast<ArrayRCP<index_t> & > (separatorRange_);
  }
  
  /*! \brief Get separator tree by const RCP.
   */
  inline ArrayRCP<index_t> &getSeparatorTreeRCPConst() const
  {
    return const_cast<ArrayRCP<index_t> & > (separatorTree_);
  }

  /*! \brief Get pointer to permutation.
   *  If inverse = true, return inverse permutation.
   *  By default, perm[i] is where new index i can be found in the old ordering.
   *  When inverse==true, perm[i] is where old index i can be found 
   *  in the new ordering.
   */
  inline index_t *getPermutationView(bool inverse = false) const
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
  
  /*! \brief Get pointer to separator range.
   */
  inline index_t *getSeparatorRangeView() const
  {
    // Here, don't need to check perm_size_ before calling getRawPtr.
    // separatorRange_ always has some length, since it is allocated larger
    // than other arrays.
    return separatorRange_.getRawPtr();
  }

  /*! \brief Get pointer to separator tree.
   */
  inline index_t *getSeparatorTreeView() const
  {
    if (perm_size_)
      return separatorTree_.getRawPtr();
    else
      return NULL;
  }
  
  /*! \brief Get reference to separator column block.
   */
  inline index_t &NumSeparatorBlocks()
  {
    return separatorColBlocks_; 
  }

  /*! \brief returns 0 if permutation is valid, negative if invalid.
   */
  int validatePerm()
  {
    size_t n = getPermutationSize();

    if(!havePerm_) {
      return -1;
    }

    std::vector<int> count(getPermutationSize());
    for (size_t i = 0; i < n; i++) {
      count[i] = 0;
    }

    for (size_t i = 0; i < n; i++) {
      if ( (perm_[i] < 0) || (perm_[i] >= static_cast<index_t>(n)) )
        return -1;
      else
        count[perm_[i]]++;
    }

    // Each index should occur exactly once (count==1)
    for (size_t i = 0; i < n; i++) {
      if (count[i] != 1){
        return -2;
      }
    }

    return 0;
  }

protected:
  // Ordering solution consists of permutation vector(s).
  // Either perm or invperm should be computed by the algorithm.
  size_t perm_size_;
  bool havePerm_;                    // has perm_ been computed yet?
  bool haveInverse_;                 // has invperm_ been computed yet?
  bool haveSeparatorRange_;          // has sepRange_ been computed yet?
  bool haveSeparatorTree_;           // has sepTree_ been computed yet?
  ArrayRCP<index_t> perm_;           // zero-based local permutation
  ArrayRCP<index_t> invperm_;        // inverse of permutation above
  ArrayRCP<index_t> separatorRange_; // range iterator for separator tree
  ArrayRCP<index_t> separatorTree_;  // separator tree
  index_t separatorColBlocks_;       // number of column blocks in separator
};

template <typename lno_t>
class LocalOrderingSolution : public OrderingSolution<lno_t>
{
public:
  LocalOrderingSolution(lno_t perm_size) : OrderingSolution<lno_t>(perm_size) {}
};

template <typename gno_t>
class GlobalOrderingSolution : public OrderingSolution<gno_t>
{
public:
  GlobalOrderingSolution(gno_t perm_size) : OrderingSolution<gno_t>(perm_size) {}
};

}

#endif
