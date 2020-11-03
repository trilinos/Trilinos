#if 0
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
