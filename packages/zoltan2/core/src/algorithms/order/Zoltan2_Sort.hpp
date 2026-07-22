// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
