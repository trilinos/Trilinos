/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#ifndef TPETRA_DETAILS_HASH_HPP
#define TPETRA_DETAILS_HASH_HPP

#include "Tpetra_ConfigDefs.hpp"
#ifdef TPETRA_USE_MURMUR_HASH
#  include <Kokkos_Functional.hpp> // hash function used by Kokkos::UnorderedMap
#endif // TPETRA_USE_MURMUR_HASH
#include <type_traits> // make_signed

namespace Tpetra {
namespace Details {

namespace Impl {

//! Implementation of getRecommendedSize (see below) for int result_type.
int getRecommendedSizeInt (const int size);

} // namespace Impl

/// \brief The hash function for FixedHashTable.
///
/// \tparam KeyType The type of hash function inputs.  Corresponds to
///   Tpetra's GlobalOrdinal type.
/// \tparam DeviceType Kokkos::Device specialization.
/// \tparam OffsetType Type of each entry of the "buckets" (offsets)
///   array \c ptr_ in FixedHashTable.
/// \tparam ResultType Type of the return value of the hash function.
template<class KeyType,
         class DeviceType,
         class OffsetType = typename std::make_signed<typename Kokkos::View<KeyType*, DeviceType>::size_type>::type,
         class ResultType = int>
struct Hash {
  /// \brief Type of the hash function's input.
  ///
  /// This is the same typedef as found in std::hash.
  typedef KeyType argument_type;

  /// \brief Type of the return value of the hash function.
  ///
  /// This is the same typedef as found in std::hash.
  typedef ResultType result_type;

  //! Type of offsets into the hash table's array of (key,value) pairs.
  typedef OffsetType offset_type;

  /// \brief The hash function.
  ///
  /// \param key [in] The key to hash.
  /// \param size [in] Number of "buckets" in the hash table.
  ///
  /// It is legal to call this function in a Kokkos parallel kernel.
  static KOKKOS_INLINE_FUNCTION result_type
  hashFunc (const argument_type& /*key*/, const offset_type& /*size*/) {
    static_assert (! std::is_same<result_type, int>::value,
                   "Not yet implemented for ResultType != int");
  }

  /// \brief Number of "buckets" that the constructor of
  ///   FixedHashTable should allocate.
  ///
  /// \param size [in] Number of entries in the table.
  ///   (FixedHashTable fixes the number of entries in its
  ///   constructor.)
  ///
  /// This is included in Hash, because a correct and efficient
  /// implementation is a function of result_type.  The point is to
  /// factor out all of the non-generic stuff into this Hash struct,
  /// which can be specialized for the types of interest.
  static result_type getRecommendedSize (const offset_type /*size*/) {
    static_assert (! std::is_same<result_type, int>::value,
                   "Not yet implemented for ResultType != int");
  }
};

/// \brief Specialization for ResultType = int.
///
/// \tparam KeyType The type of hash function inputs.  Corresponds to
///   Tpetra's GlobalOrdinal type.
/// \tparam DeviceType Kokkos::Device specialization.
/// \tparam OffsetType Type of each entry of the "buckets" (offsets)
///   array \c ptr_ in FixedHashTable.
///
/// This hash function currently always returns \c int, no matter the
/// value type.  This is unfortunate, because it limits Tpetra to have
/// <tt>LocalOrdinal = int</tt> (or smaller).  I would like to fix
/// this at some point.  Tpetra provides this partial specialization
/// as a hook for fixing this later.
template<class KeyType, class DeviceType, class OffsetType>
struct Hash<KeyType, DeviceType, OffsetType, int> {
  /// \brief Type of the hash function's input.
  ///
  /// This is the same typedef as found in std::hash.
  typedef KeyType argument_type;

  /// \brief Type of the return value of the hash function.
  ///
  /// This is the same typedef as found in std::hash.
  typedef int result_type;

  //! Type of offsets into the hash table's array of (key,value) pairs.
  typedef OffsetType offset_type;

  /// \brief The hash function.
  ///
  /// \param key [in] The key to hash.
  /// \param size [in] Number of "buckets" in the hash table.
  ///
  /// It is legal to call this function in a Kokkos parallel kernel.
  static KOKKOS_INLINE_FUNCTION result_type
  hashFunc (const argument_type& key, const offset_type& size)
  {
#ifdef TPETRA_USE_MURMUR_HASH
    Kokkos::pod_hash<argument_type> hash;
    const uint32_t k = hash (key);
    return static_cast<result_type> (k % size);
#else
    // We are using Epetra's hash function by default, as we have
    // observed that it is much faster than the Murmur hash
    // function. However, this is not a good hash function for general
    // sets of keys.  For our typical use case, this is good.  Use
    // Murmur hash if the maps are sparse.
    const unsigned int seed = (2654435761U);
    const int intkey = (int) ((key & 0x000000007fffffffLL) +
                              ((key & 0x7fffffff80000000LL) >> 31));
    return static_cast<result_type> ((seed ^ intkey) % static_cast<int> (size));
#endif
  }

  /// \brief Number of "buckets" that the constructor of
  ///   FixedHashTable should allocate.
  ///
  /// \param size [in] Number of entries in the table.
  ///   (FixedHashTable fixes the number of entries in its
  ///   constructor.)
  static result_type getRecommendedSize (const offset_type size)
  {
    return Impl::getRecommendedSizeInt (static_cast<int> (size));
  }
};

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_HASH_HPP
