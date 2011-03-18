// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_HASH_HPP_
#define _ZOLTAN2_HASH_HPP_

/*! \file Zoltan2_Hash.hpp
*/

#include <Teuchos_HashUtils.hpp>

namespace Z2
{

template <typename T1, typename T2>
T2 Hash(T1 &key, T2 n)
{
  int ikey = Teuchos::hashCode(key);

  return T2(ikey % n);
}

}                   // namespace Z2

#endif
