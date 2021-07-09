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
#ifndef _ZOLTAN2_SORT_HPP_
#define _ZOLTAN2_SORT_HPP_

#include <Teuchos_Array.hpp>
#include <algorithm>


////////////////////////////////////////////////////////////////////////
//! \file Zoltan2_Sort.hpp
//! \brief Sort vector of pairs (key, value) by value. 
//! \brief This class is needed so we also get the sorted keys (indices).
  
// TODO: This is a generic utility class; should move this source file.
// We could perhaps use Sort2 from Tpetra, but that uses a custom sort not std::sort

namespace Zoltan2{

namespace Details {

  // These helper functions can't be class members, so put in Details namespace.
  template<class KeyType, class ValueType>
  bool SortInc (const std::pair<KeyType, ValueType>& a, const std::pair<KeyType, ValueType>& b)
  {
    return a.second < b.second;
  }

  template<class KeyType, class ValueType>
  bool SortDec (const std::pair<KeyType, ValueType>& a, const std::pair<KeyType, ValueType>& b)
  {
    return a.second > b.second;
  }

} // namespace Details

template <typename key_t, typename value_t>
class SortPairs
{
  public:
    SortPairs()
    {
    }

  public:
    void sort(Teuchos::Array<std::pair<key_t,value_t> > &listofPairs, bool inc=true)
    {
      // Sort in increasing (default) or decreasing order of value
      if (inc)
        std::stable_sort (listofPairs.begin(), listofPairs.end(), Details::SortInc<key_t, value_t>);
      else
        std::stable_sort (listofPairs.begin(), listofPairs.end(), Details::SortDec<key_t, value_t>);
    }

};
}
#endif
