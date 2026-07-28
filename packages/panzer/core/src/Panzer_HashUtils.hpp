// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// *******************************************************************
// This file contains a copy of hash support code from boost that
// didn't make it into the stl.  We only needed two lines code so
// copied it here. Below is boost copyright.
// *******************************************************************

// Copyright 2005-2014 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  Based on Peter Dimov's proposal
//  http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2005/n1756.pdf
//  issue 6.18. 
//
//  This also contains public domain code from MurmurHash. From the
//  MurmurHash header:

// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.

// ******************************************************************* 
// ******************************************************************* 

#ifndef PANZER_HASH_UTILS_HPP
#define PANZER_HASH_UTILS_HPP

namespace panzer {

  /** \brief Mixes the hash of v into seed, in place. Used to build a combined
    * hash from multiple values.
    *
    * \tparam T the type of the value being hashed; must be usable with std::hash.
    * \param[in,out] seed the running hash to mix v's hash into.
    * \param[in] v the value whose hash should be mixed in.
    */
  template <class T>
  inline void hash_combine(std::size_t& seed, const T& v)
  {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
  }

  /// \brief A hash functor for std::pair, combining the hashes of both elements via hash_combine(). Useful as the Hash template argument of an unordered container keyed on std::pair.
  struct pair_hash
  {
    template<typename T1, typename T2>
    std::size_t operator()(const std::pair<T1,T2>& v) const
    {
      std::size_t seed = 0;
      panzer::hash_combine(seed, v.first);
      panzer::hash_combine(seed, v.second);
      return seed;
    }
  };

}

namespace std
{
/// \brief Specializes std::hash for std::pair, combining the hashes of both elements via panzer::hash_combine(). Allows std::pair to be used directly as a key in std::unordered_map/std::unordered_set.
template <typename T1, typename T2>
struct hash<std::pair<T1,T2> >
{
  std::size_t operator()(const std::pair<T1,T2>& v) const
  {
    std::size_t seed = 0;
    panzer::hash_combine(seed, v.first);
    panzer::hash_combine(seed, v.second);
    return seed;
  }
};
}

#endif
