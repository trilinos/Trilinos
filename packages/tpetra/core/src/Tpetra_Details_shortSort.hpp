// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_SHORTSORT_HPP
#define TPETRA_DETAILS_SHORTSORT_HPP

/// \file Tpetra_Details_shortSort.hpp
/// \brief Declaration and definition of functions for sorting "short"
///   arrays of keys and corresponding values.
///
/// "Short" means "having few entries."  This matters because optimal
/// (in terms of fewest comparisons) sorting algorithms are different
/// for short arrays.

#include "TpetraCore_config.h"
#include "Kokkos_Macros.hpp"
#include <type_traits>

namespace Tpetra {
namespace Details {

// Make sure that the macro defined below wasn't defined somewhere else.
#ifdef TPETRA_DETAILS_SWAP_KEYSANDVALUES
#  error "The TPETRA_DETAILS_SWAP_KEYSANDVALUES macro is already defined."
#endif // TPETRA_DETAILS_SWAP_KEYSANDVALUES

/// \brief Macro that swaps the i and j entries of keys and values, if
///   keys[i] > keys[j] (i.e., if keys[i] and keys[j] are out of order).
///
/// \warning This macro is an implementation detail of
///   shortSortKeysAndValues_2, shortSortKeysAndValues_3,
///   shortSortKeysAndValues_4, and shortSortKeysAndValues_8 (see
///   below).  Do not rely on this macro anywhere in your code.  The
///   macro must be called inside those functions, because it assumes
///   that certain symbols exist (see list below).
///
/// \param i [in] Integer index.  i <= j.
/// \param j [in] Integer index.  i <= j.
///
/// Captured symbols:
/// <ul>
/// <li> \c keys: Array of keys; each entry has type \c KeyType </li>
/// <li> \c values: Array of values; each entry has type \c ValueType </li>
/// <li> \c KeyType: Type of each entry of \c keys </li>
/// <li> \c ValueType: Type of each entry of \c values </li>
/// </ul>
#define TPETRA_DETAILS_SWAP_KEYSANDVALUES( i, j )  \
  if (keys[i] > keys[j]) { \
    const KeyType tmpKey (keys[i]); \
    keys[i] = keys[j]; \
    keys[j] = tmpKey; \
    const ValueType tmpVal (values[i]); \
    values[i] = values[j]; \
    values[j] = tmpVal; \
  }

// Make sure that the macro defined below wasn't defined somewhere else.
#ifdef TPETRA_DETAILS_SWAP_KEYS
#  error "The TPETRA_DETAILS_SWAP_KEYS macro is already defined."
#endif // TPETRA_DETAILS_SWAP_KEYSANDVALUES

/// \brief Macro that swaps the i and j entries of keys, if keys[i] >
///   keys[j] (i.e., if keys[i] and keys[j] are out of order).
///
/// \warning This macro is an implementation detail of
///   shortSortKeys_2, shortSortKeys_3, shortSortKeys_4, and
///   shortSortKeys_8 (see below).  Do not rely on this macro anywhere
///   in your code.  The macro must be called inside those functions,
///   because it assumes that certain symbols exist (see list below).
///
/// \param i [in] Integer index.  i <= j.
/// \param j [in] Integer index.  i <= j.
///
/// Captured symbols:
/// <ul>
/// <li> \c keys: Array of keys; each entry has type \c KeyType </li>
/// <li> \c KeyType: Type of each entry of \c keys </li>
/// </ul>
#define TPETRA_DETAILS_SWAP_KEYS( i, j )  \
  if (keys[i] > keys[j]) { \
    const KeyType tmpKey (keys[i]); \
    keys[i] = keys[j]; \
    keys[j] = tmpKey; \
  }

/// \brief Sort keys and values jointly, by keys, for arrays of length 2.
///
/// \tparam KeyType Greater-than comparable, copy constructible,
///   assignable
/// \tparam ValueType Copy constructible, assignable
///
/// \param keys [in/out] Length 2 array of keys.  This function
///   sorts this \c keys array, and applies the same permutation to
///   the \c values array.
/// \param values [in/out] Length 2 array of values.
template<class KeyType, class ValueType>
KOKKOS_FUNCTION void
shortSortKeysAndValues_2 (KeyType keys[2], ValueType values[2])
{
  // Since this function takes a constant number of entries, I use a
  // sorting network here.  For 2 entries, the sorting network is
  // nearly trivial.
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(0, 1);
}

/// \brief Sort length-2 array of keys.
///
/// \tparam KeyType Greater-than comparable, copy constructible,
///   assignable
///
/// \param keys [in/out] Length-2 array of keys to sort.
template<class KeyType>
KOKKOS_FUNCTION void
shortSortKeys_2 (KeyType keys[2])
{
  // Since this function takes a constant number of entries, I use a
  // sorting network here.  For 2 entries, the sorting network is
  // nearly trivial.
  TPETRA_DETAILS_SWAP_KEYS(0, 1);
}

/// \brief Sort keys and values jointly, by keys, for arrays of length 3.
///
/// \tparam KeyType Greater-than comparable, copy constructible,
///   assignable
/// \tparam ValueType Copy constructible, assignable
///
/// \param keys [in/out] Length 3 array of keys.  This function
///   sorts this \c keys array, and applies the same permutation to
///   the \c values array.
/// \param values [in/out] Length 3 array of values.
template<class KeyType, class ValueType>
KOKKOS_FUNCTION void
shortSortKeysAndValues_3 (KeyType keys[3], ValueType values[3])
{
  // Since this function takes a constant number of entries, I use a
  // sorting network here.  To make the network, I used the generator
  // at
  //
  // http://pages.ripco.net/~jgamble/nw.html
  //
  // "Best" algorithm in this case is the Bose-Nelson Algorithm.
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(1, 2);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(0, 2);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(0, 1);
}

/// \brief Sort length-3 array of keys.
///
/// \tparam KeyType Greater-than comparable, copy constructible,
///   assignable
///
/// \param keys [in/out] Length-3 array of keys to sort.
template<class KeyType>
KOKKOS_FUNCTION void
shortSortKeys_3 (KeyType keys[3])
{
  // Since this function takes a constant number of entries, I use a
  // sorting network here.  To make the network, I used the generator
  // at
  //
  // http://pages.ripco.net/~jgamble/nw.html
  //
  // "Best" algorithm in this case is the Bose-Nelson Algorithm.
  TPETRA_DETAILS_SWAP_KEYS(1, 2);
  TPETRA_DETAILS_SWAP_KEYS(0, 2);
  TPETRA_DETAILS_SWAP_KEYS(0, 1);
}

/// \brief Sort keys and values jointly, by keys, for arrays of length 4.
///
/// \tparam KeyType Greater-than comparable, copy constructible,
///   assignable
/// \tparam ValueType Copy constructible, assignable
///
/// \param keys [in/out] Length 4 array of keys.  This function
///   sorts this \c keys array, and applies the same permutation to
///   the \c values array.
/// \param values [in/out] Length 4 array of values.
template<class KeyType, class ValueType>
KOKKOS_FUNCTION void
shortSortKeysAndValues_4 (KeyType keys[4], ValueType values[4])
{
  // Since this function takes a constant number of entries, I use a
  // sorting network here.  To make the network, I used the generator
  // at
  //
  // http://pages.ripco.net/~jgamble/nw.html
  //
  // "Best" algorithm in this case is the Bose-Nelson Algorithm.
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(0, 1);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(2, 3);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(0, 2);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(1, 3);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(1, 2);
}

/// \brief Sort length-4 array of keys.
///
/// \tparam KeyType Greater-than comparable, copy constructible,
///   assignable
///
/// \param keys [in/out] Length-4 array of keys to sort.
template<class KeyType>
KOKKOS_FUNCTION void
shortSortKeys_4 (KeyType keys[4])
{
  // Since this function takes a constant number of entries, I use a
  // sorting network here.  To make the network, I used the generator
  // at
  //
  // http://pages.ripco.net/~jgamble/nw.html
  //
  // "Best" algorithm in this case is the Bose-Nelson Algorithm.
  TPETRA_DETAILS_SWAP_KEYS(0, 1);
  TPETRA_DETAILS_SWAP_KEYS(2, 3);
  TPETRA_DETAILS_SWAP_KEYS(0, 2);
  TPETRA_DETAILS_SWAP_KEYS(1, 3);
  TPETRA_DETAILS_SWAP_KEYS(1, 2);
}

/// \brief Sort keys and values jointly, by keys, for arrays of length 8.
///
/// \tparam KeyType Greater-than comparable, copy constructible,
///   assignable
/// \tparam ValueType Copy constructible, assignable
///
/// \param keys [in/out] Length 8 array of keys.  This function
///   sorts this \c keys array, and applies the same permutation to
///   the \c values array.
/// \param values [in/out] Length 8 array of values.
template<class KeyType, class ValueType>
KOKKOS_FUNCTION void
shortSortKeysAndValues_8 (KeyType keys[8], ValueType values[8])
{
  // Since this function takes a constant number of entries, I use a
  // sorting network here.  To make the network, I used the generator
  // at
  //
  // http://pages.ripco.net/~jgamble/nw.html
  //
  // "Best" algorithm in this case is the Bose-Nelson Algorithm.
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(0, 1);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(2, 3);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(0, 2);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(1, 3);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(1, 2);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(4, 5);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(6, 7);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(4, 6);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(5, 7);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(5, 6);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(0, 4);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(1, 5);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(1, 4);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(2, 6);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(3, 7);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(3, 6);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(2, 4);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(3, 5);
  TPETRA_DETAILS_SWAP_KEYSANDVALUES(3, 4);
}

/// \brief Sort length-8 array of keys.
///
/// \tparam KeyType Greater-than comparable, copy constructible,
///   assignable
///
/// \param keys [in/out] Length-8 array of keys to sort.
template<class KeyType>
KOKKOS_FUNCTION void
shortSortKeys_8 (KeyType keys[8])
{
  // Since this function takes a constant number of entries, I use a
  // sorting network here.  To make the network, I used the generator
  // at
  //
  // http://pages.ripco.net/~jgamble/nw.html
  //
  // "Best" algorithm in this case is the Bose-Nelson Algorithm.
  TPETRA_DETAILS_SWAP_KEYS(0, 1);
  TPETRA_DETAILS_SWAP_KEYS(2, 3);
  TPETRA_DETAILS_SWAP_KEYS(0, 2);
  TPETRA_DETAILS_SWAP_KEYS(1, 3);
  TPETRA_DETAILS_SWAP_KEYS(1, 2);
  TPETRA_DETAILS_SWAP_KEYS(4, 5);
  TPETRA_DETAILS_SWAP_KEYS(6, 7);
  TPETRA_DETAILS_SWAP_KEYS(4, 6);
  TPETRA_DETAILS_SWAP_KEYS(5, 7);
  TPETRA_DETAILS_SWAP_KEYS(5, 6);
  TPETRA_DETAILS_SWAP_KEYS(0, 4);
  TPETRA_DETAILS_SWAP_KEYS(1, 5);
  TPETRA_DETAILS_SWAP_KEYS(1, 4);
  TPETRA_DETAILS_SWAP_KEYS(2, 6);
  TPETRA_DETAILS_SWAP_KEYS(3, 7);
  TPETRA_DETAILS_SWAP_KEYS(3, 6);
  TPETRA_DETAILS_SWAP_KEYS(2, 4);
  TPETRA_DETAILS_SWAP_KEYS(3, 5);
  TPETRA_DETAILS_SWAP_KEYS(3, 4);
}

/// \brief Shellsort (yes, it's one word) the input array \c keys, and
///   apply the resulting permutation to the input array \c values.
///
/// mfh 28 Nov 2016, 17 Dec 2016: I adapted this function from
/// sh_sort2 in Tpetra_Util.hpp (in this directory).
template<class KeyType, class ValueType, class IndexType>
KOKKOS_FUNCTION void
shellSortKeysAndValues (KeyType keys[],
                        ValueType values[],
                        const IndexType n)
{
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a signed integer type.");
  static_assert (std::is_signed<IndexType>::value,
                 "IndexType must be a signed integer type.  "
                 "This implementation does a count-down loop, "
                 "and may thus loop forever "
                 "if one attempts to use it with unsigned IndexType.");
  constexpr IndexType ZERO = 0;
  IndexType midpoint = n / static_cast<IndexType> (2);

  while (midpoint > ZERO) {
    // Avoid names like "max" in case they collide with macros.
    const IndexType theMax = n - midpoint;
    for (IndexType j = 0; j < theMax; ++j) {
      // IndexType is signed, so it's legit to compare >= 0.
      for (IndexType k = j; k >= 0; k -= midpoint) {
        if (keys[k + midpoint] >= keys[k]) {
          break;
        }
        const KeyType tmpKey = keys[k + midpoint];
        keys[k + midpoint] = keys[k];
        keys[k] = tmpKey;
        const ValueType tmpVal = values[k + midpoint];
        values[k + midpoint] = values[k];
        values[k] = tmpVal;
      }
    }
    midpoint = midpoint / 2;
  }
}

/// \brief Shellsort (yes, it's one word) the input array \c keys.
///
/// \param keys [in/out] Input array of keys to sort.
/// \param n [in] Length of the input array \c keys.
template<class KeyType, class IndexType>
KOKKOS_FUNCTION void
shellSortKeys (KeyType keys[], const IndexType n)
{
  static_assert (std::is_integral<IndexType>::value,
                 "IndexType must be a signed integer type.");
  static_assert (std::is_signed<IndexType>::value,
                 "IndexType must be a signed integer type.  "
                 "This implementation does a count-down loop, "
                 "and may thus loop forever "
                 "if one attempts to use it with unsigned IndexType.");
  constexpr IndexType ZERO = 0;
  IndexType midpoint = n / static_cast<IndexType> (2);

  while (midpoint > ZERO) {
    // Avoid names like "max" in case they collide with macros.
    const IndexType theMax = n - midpoint;
    for (int j = 0; j < theMax; ++j) {
      // IndexType must be signed, so it's legit to compare >= 0.
      for (int k = j; k >= 0; k -= midpoint) {
        if (keys[k + midpoint] >= keys[k]) {
          break;
        }
        const KeyType tmpKey = keys[k + midpoint];
        keys[k + midpoint] = keys[k];
        keys[k] = tmpKey;
      }
    }
    midpoint = midpoint / 2;
  }
}

template<typename KeyType, typename ValueType, typename IndexType>
KOKKOS_FUNCTION void
shellSortKeysAndValues2(KeyType* keys, ValueType* values, IndexType n)
{
  IndexType ngaps = 10;
  //Use Ciura's gap choices
  IndexType gaps[] = {3548, 1577, 701, 301, 132, 57, 23, 10, 4, 1};
  for(IndexType gapIndex = 0; gapIndex < ngaps; gapIndex++)
  {
    auto gap = gaps[gapIndex];
    if(n < gap)
      continue;
    //insertion sort the array {keys[0*gap], keys[1*gap], keys[2*gap], ...}
    for(IndexType gapOffset = 0; gapOffset < gap; gapOffset++)
    {
      for(IndexType i = gap + gapOffset; i < n; i += gap)
      {
        //avoid extra swaps: scan for the final position of keys[i]
        if(keys[i - gap] > keys[i])
        {
          IndexType finalPos = i - gap;
          while(finalPos - gap >= 0 && keys[finalPos - gap] > keys[i])
          {
            finalPos -= gap;
          }
          //save keys/values [i], then shift up all keys/values between finalPos and i-gap (inclusive)
          auto tempKey = keys[i];
          auto tempVal = values[i];
          for(IndexType j = i - gap; j >= finalPos; j -= gap)
          {
            keys[j + gap] = keys[j];
            values[j + gap] = values[j];
          }
          keys[finalPos] = tempKey;
          values[finalPos] = tempVal;
        }
      }
    }
  }
#undef SHELL_SWAP
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_SHORTSORT_HPP
