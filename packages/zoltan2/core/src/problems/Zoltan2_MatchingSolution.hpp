#if 0
// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_MatchingSolution.hpp
    \brief Defines the MatchingSolution class.
*/

#ifndef _ZOLTAN2_MATCHINGSOLUTION_HPP_
#define _ZOLTAN2_MATCHINGSOLUTION_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Solution.hpp>

namespace Zoltan2 {

/*! \brief The class containing coloring solution.

    Template parameters:
    \li \c adapter    input adapter

The coloring solution contains an array of colors, one per id.
Matchs are represented as int (sufficient for any reasonable use case). 
A special value, currently 0, is used for vertices that have not been colored.

*/

template <typename Adapter>
  class MatchingSolution : public Solution
{
private: 
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::lno_t lno_t;

public:

  /*! \brief Constructor allocates memory for the solution.
   */
  MatchingSolution(
    size_t length // This should be equal to nlids. TODO: Optional?
  )
  {
    HELLO;
    length_ = length;
    colors_  = ArrayRCP<int>(length_);
  }

  //////////////////////////////////////////////
  // Accessor functions, allowing algorithms to get ptrs to solution memory.
  // Algorithms can then load the memory.
  // Non-RCP versions are provided for applications to use.

  /*! \brief Get (local) size of color array.
   */
  inline size_t getMatchsSize() {return length_;}

  /*! \brief Get (local) color array by RCP.
   */
  inline ArrayRCP<int>  &getMatchsRCP()  {return colors_;}

  /*! \brief Get (local) color array by raw pointer (no RCP).
   */
  inline int * getMatchs()  {return &(*colors_);}

  /*! \brief Get local number of colors.
   *  This is computed from the coloring each time, as this is cheap.
   */
  int getNumColors()  
  { 
    int maxColor = 0;
    for (size_t i=0; i<length_; i++){
      if (colors_[i] > maxColor)
        maxColor = colors_[i];
    }
    return maxColor;
  } 

  /*! \brief Get global number of colors.
   */
  //int getGlobalNumColors(); // TODO
 
protected:
  // Matching solution consists of permutation vector(s).
  size_t length_;
  ArrayRCP<int> colors_;   // zero-based local color array
  //int numMatchs_;        // Number of colors (local on this proc)
  //int numMatchsGlobal_;  // For future distributed coloring
};

}

#endif
#endif
