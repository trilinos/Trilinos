// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_UTIL_HPP_
#define _ZOLTAN2_UTIL_HPP_

/*! \file Zoltan2_Util.hpp
*/

#include <Zoltan2_IdentifierTraits.hpp>

namespace Z2
{

template <typename T1, typename T2>
inline T2 Hash(T1 &key, T2 n)
{
  return T2(IdentifierTraits<T1>(key) % n);
}


}                   // namespace Z2

#endif
