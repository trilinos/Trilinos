// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Details_shortSort.hpp"
#include <iterator>
#include <utility> // std::swap

namespace { // (anonymous)

using Tpetra::Details::shortSortKeys_2;
using Tpetra::Details::shortSortKeys_3;
using Tpetra::Details::shortSortKeys_4;
using Tpetra::Details::shortSortKeys_8;
using std::endl;

// Shell sort the input array 'keys' (length n).
//
// mfh 17 Dec 2016: I adapted this function from
// shellSortKeysAndValues in Sort.cpp (in this directory).
template<class KeyType>
void
shellSortKeys (KeyType keys[],
               const int n)
{
  const int ZERO = 0;
  int midpoint = n / 2;

  while (midpoint > ZERO) {
    // Avoid names like "max" in case they collide with macros.
    const int theMax = n - midpoint;
    for (int j = 0; j < theMax; ++j) {
      // int is signed, so it's legit to compare >= 0.
      for (int k = j; k >= 0; k -= midpoint) {
        if (keys[k + midpoint] >= keys[k]) {
          break;
        }
        std::swap (keys[k + midpoint], keys[k]);
      }
    }
    midpoint = midpoint / 2;
  }
}

// Is the given range of keys sorted?
template<class KeyType>
bool
isSorted (const KeyType keys[], const int len)
{
  if (len <= 1) {
    return true;
  }
  else {
    for (int k = 1; k < len; ++k) {
      if (keys[k] < keys[k-1]) { // only use less-than compare
        return false;
      }
    }
    return true;
  }
}

// Copy 'in' into 'out'.
template<class EntryType>
void
copyArray (EntryType out[], const EntryType in[], const int len)
{
  for (int k = 0; k < len; ++k) {
    out[k] = in[k];
  }
}

// Are all the corresponding entries of the two arrays 'x' and 'y' equal?
template<class EntryType>
bool
isEqual (const EntryType x[], const EntryType y[], const size_t len)
{
  for (size_t k = 0; k < len; ++k) {
    if (x[k] != y[k]) {
      return false;
    }
  }
  return true;
}

// Convert run-time choice of array length, into compile-time choice.
template<class KeyType, const int arrayLength>
struct ShortSortKeys {
  // Call the appropriate shortSortKeys_* function, if it exists for
  // the given arrayLength.  If it does not exist, throw
  // std::logic_error.
  static void call (KeyType /* keys */ []) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Not implemented for arrayLength = "
       << arrayLength << ".");
  }
};

template<class KeyType>
struct ShortSortKeys<KeyType, 2> {
  static void call (KeyType keys[2]) {
    Tpetra::Details::shortSortKeys_2 (keys);
  }
};

template<class KeyType>
struct ShortSortKeys<KeyType, 3> {
  static void call (KeyType keys[3]) {
    Tpetra::Details::shortSortKeys_3 (keys);
  }
};

template<class KeyType>
struct ShortSortKeys<KeyType, 4> {
  static void call (KeyType keys[4]) {
    Tpetra::Details::shortSortKeys_4 (keys);
  }
};

template<class KeyType>
struct ShortSortKeys<KeyType, 8> {
  static void call (KeyType keys[8]) {
    Tpetra::Details::shortSortKeys_8 (keys);
  }
};

// \brief Test shortSortKeys_${arrayLength} for fixed \c KeyType type,
//   and a fixed \c arrayLength value.
//
// \tparam KeyType Type of each entry of the \c keys array.
// \tparam arrayLength Length of the keys array.
//
// \param success [out] Whether the test is successful.
// \param out [out] Output stream.
// \param keys [in] Input array of keys for the sort.  The test must
//   make a deep copy of this before sorting it.
template<class KeyType, const int arrayLength>
void
test_fixedTypes_fixedArrayLength (bool& success,
                                  Teuchos::FancyOStream& out,
                                  const KeyType keys[arrayLength])
{
  out << "Test shortSortKeys_" << arrayLength << endl;
  Teuchos::OSTab tab1 (out);

  // Copy the input, since we'll need the original input arrays later.
  KeyType keysCopy[arrayLength];
  copyArray (keysCopy, keys, arrayLength);

  // Call the sort function.
  typedef ShortSortKeys<KeyType, arrayLength> impl_type;
  TEST_NOTHROW( impl_type::call (keysCopy) );

  // Make sure the keys got sorted.
  const bool sorted = isSorted (keysCopy, arrayLength);
  TEST_ASSERT( sorted );

  // Compare against the above shell sort implementation.
  KeyType keysCopy2[arrayLength];
  copyArray (keysCopy2, keys, arrayLength);
  shellSortKeys (keysCopy2, arrayLength);
  const bool equalKeys = isEqual (keysCopy, keysCopy2, arrayLength);
  TEST_ASSERT( equalKeys );

  // Compare against Tpetra::Details::shellSortKeys.
  KeyType keysCopy3[arrayLength];
  copyArray (keysCopy3, keys, arrayLength);
  ::Tpetra::Details::shellSortKeys (keysCopy3, arrayLength);
  const bool equalKeys3 = isEqual (keysCopy2, keysCopy3, arrayLength);
  TEST_ASSERT( equalKeys3 );
}

// \brief Test shortSortKeys_${arrayLength} for fixed \c KeyType type,
//   and for all implemented \c arrayLength values.
//
// \tparam KeyType Type of each entry of the \c keys array.
//
// \param keyTypeName [in] Human-readable name of \c KeyType.
template<class KeyType>
void
test_fixedTypes (bool& success,
                 Teuchos::FancyOStream& out,
                 const char keyTypeName[])
{
  out << "Test shortSortKeys_*" << " for KeyType=" << keyTypeName << endl;
  Teuchos::OSTab tab1 (out);

  {
    constexpr int arrayLength = 2;
    out << "Test arrayLength=" << arrayLength << endl;
    Teuchos::OSTab tab2 (out);

    {
      KeyType keys[arrayLength] = {0, 1};
      out << "Test in-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 0};
      out << "Test out-of-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 1};
      out << "Test duplicate keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
  }

  {
    constexpr int arrayLength = 3;
    out << "Test arrayLength=" << arrayLength << endl;
    Teuchos::OSTab tab2 (out);

    {
      KeyType keys[arrayLength] = {1, 3, 5};
      out << "Test in-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 1};
      out << "Test reverse-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {3, 5, 1};
      out << "Test out-of-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3};
      out << "Test in-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 3, 3};
      out << "Test in-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {3, 3, 1};
      out << "Test out-of-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 3};
      out << "Test out-of-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 5};
      out << "Test out-of-order keys with separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
  }

  {
    constexpr int arrayLength = 4;
    out << "Test arrayLength=" << arrayLength << endl;
    Teuchos::OSTab tab2 (out);

    {
      KeyType keys[arrayLength] = {1, 3, 5, 7};
      out << "Test in-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {7, 5, 3, 1};
      out << "Test reverse-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {3, 7, 5, 1};
      out << "Test out-of-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 1, 1};
      out << "Test in-order keys, all duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3, 5};
      out << "Test in-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 3, 5, 5};
      out << "Test in-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {5, 5, 3, 1};
      out << "Test out-of-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 1, 1};
      out << "Test out-of-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 1, 5};
      out << "Test out-of-order keys with separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3, 3};
      out << "Test in-order keys with multiple duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {3, 3, 1, 1};
      out << "Test out-of-order keys with multiple duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {3, 1, 3, 1};
      out << "Test out-of-order keys with multiple, separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
  }

  {
    constexpr int arrayLength = 8;
    out << "Test arrayLength=" << arrayLength << endl;
    Teuchos::OSTab tab2 (out);

    {
      KeyType keys[arrayLength] = {1, 3, 5, 7, 9, 11, 13, 15};
      out << "Test in-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {15, 13, 11, 9, 7, 5, 3, 1};
      out << "Test reverse-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {15, 11, 3, 7, 5, 1, 13, 9};
      out << "Test out-of-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 1, 1, 1, 1, 1, 1};
      out << "Test in-order keys, all duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3, 5, 7, 9, 11, 13};
      out << "Test in-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 3, 5, 7, 9, 11, 13, 13};
      out << "Test in-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {13, 13, 11, 9, 7, 5, 3, 1};
      out << "Test out-of-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {13, 11, 9, 7, 5, 3, 1, 1};
      out << "Test out-of-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {5, 13, 11, 9, 7, 3, 1, 5};
      out << "Test out-of-order keys with separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3, 3, 5, 5, 7, 7};
      out << "Test in-order keys with multiple duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {7, 7, 5, 5, 3, 3, 1, 1};
      out << "Test out-of-order keys with multiple duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
    {
      KeyType keys[arrayLength] = {3, 1, 3, 7, 5, 1, 9, 11};
      out << "Test out-of-order keys with multiple, separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      test_fixedTypes_fixedArrayLength<KeyType, arrayLength> (success, out, keys);
    }
  }
}

TEUCHOS_UNIT_TEST( Sort, shortSortKeys )
{
  out << "Test shortSortKeys_*" << endl;
  Teuchos::OSTab tab1 (out);

  test_fixedTypes<short> (success, out, "short");
  test_fixedTypes<int> (success, out, "int");
  test_fixedTypes<long> (success, out, "long");
  test_fixedTypes<long long> (success, out, "long long");

  test_fixedTypes<unsigned short> (success, out, "unsigned short");
  test_fixedTypes<unsigned int> (success, out, "unsigned int");
  test_fixedTypes<unsigned long> (success, out, "unsigned long");
  test_fixedTypes<unsigned long long> (success, out, "unsigned long long");
}

} // namespace (anonymous)

