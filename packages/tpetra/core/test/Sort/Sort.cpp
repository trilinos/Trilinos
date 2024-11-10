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
#include "Tpetra_Details_radixSort.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <iterator>
#include <utility> // std::swap
//numeric_limits<unsigned short> - not provided by Teuchos OrdinalTraits):
#include <limits>

namespace { // (anonymous)

using Tpetra::Details::shortSortKeysAndValues_2;
using Tpetra::Details::shortSortKeysAndValues_3;
using Tpetra::Details::shortSortKeysAndValues_4;
using Tpetra::Details::shortSortKeysAndValues_8;
using std::endl;

// Shell sort the input array 'keys' (length n), and apply the
// resulting permutation to the input array 'values'.
//
// mfh 28 Nov 2016: I adapted this function from sh_sort2 in
// Tpetra_Util.hpp (in this directory).
template<class KeyType, class ValueType>
void
shellSortKeysAndValues (KeyType keys[],
                        ValueType values[],
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
        std::swap (values[k + midpoint], values[k]);
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
template<class KeyType, class ValueType, const int arrayLength>
struct ShortSortKeysAndValues {
  // Call the appropriate shortSortKeysAndValues_* function, if it
  // exists for the given arrayLength.  If it does not exist, throw
  // std::logic_error.
  static void call (KeyType /* keys */ [], ValueType /* values */ []) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Not implemented for arrayLength = "
       << arrayLength << ".");
  }
};

template<class KeyType, class ValueType>
struct ShortSortKeysAndValues<KeyType, ValueType, 2> {
  static void call (KeyType keys[2], ValueType values[2]) {
    Tpetra::Details::shortSortKeysAndValues_2 (keys, values);
  }
};

template<class KeyType, class ValueType>
struct ShortSortKeysAndValues<KeyType, ValueType, 3> {
  static void call (KeyType keys[3], ValueType values[3]) {
    Tpetra::Details::shortSortKeysAndValues_3 (keys, values);
  }
};

template<class KeyType, class ValueType>
struct ShortSortKeysAndValues<KeyType, ValueType, 4> {
  static void call (KeyType keys[4], ValueType values[4]) {
    Tpetra::Details::shortSortKeysAndValues_4 (keys, values);
  }
};

template<class KeyType, class ValueType>
struct ShortSortKeysAndValues<KeyType, ValueType, 8> {
  static void call (KeyType keys[8], ValueType values[8]) {
    Tpetra::Details::shortSortKeysAndValues_8 (keys, values);
  }
};

// \brief Test shortSortKeysAndValues_${arrayLength} for fixed
//   \c KeyType and \c ValueType types, and a fixed \c arrayLength
//   value.
//
// \tparam KeyType Type of each entry of the \c keys array.
// \tparam ValueType Type of each entry of the \c values array.
// \tparam arrayLength Length of the keys and values arrays.
//
// \param success [out] Whether the test is successful.
// \param out [out] Output stream.
// \param keys [in] Input array of keys for the sort.  The test must
//   make a deep copy of this before sorting it.
// \param values [in] Input array of values (corresponding to the
//   above keys) for the sort.  The test must make a deep copy of this
//   before sorting it.
// \param keysMayContainDuplicates [in] Whether the input \c keys
//   array may contain duplicate entries.  This matters because if our
//   shortSortKeysAndValues_* implementations are not stable sorts,
//   then they do not promise to preserve the original ordering of the
//   entries in \c values.
// \param expectStableSort [in] Whether we expect our
//   shortSortKeysAndValues_* implementations to be stable sorts.
template<class KeyType, class ValueType, const int arrayLength>
void
test_fixedTypes_fixedArrayLength (bool& success,
                                  Teuchos::FancyOStream& out,
                                  const KeyType keys[arrayLength],
                                  const ValueType values[arrayLength],
                                  const bool keysMayContainDuplicates = true,
                                  const bool expectStableSort = false)
{
  out << "Test shortSortKeysAndValues_" << arrayLength << endl;
  Teuchos::OSTab tab1 (out);

  // Copy the input, since we'll need the original input arrays later.
  KeyType keysCopy[arrayLength];
  ValueType valuesCopy[arrayLength];
  copyArray (keysCopy, keys, arrayLength);
  copyArray (valuesCopy, values, arrayLength);

  // Call the sort function.
  typedef ShortSortKeysAndValues<KeyType, ValueType, arrayLength> impl_type;
  TEST_NOTHROW( impl_type::call (keysCopy, valuesCopy) );

  // Make sure the keys got sorted.
  const bool sorted = isSorted (keysCopy, arrayLength);
  TEST_ASSERT( sorted );

  // Compare against the above shell sort implementation.
  KeyType keysCopy2[arrayLength];
  ValueType valuesCopy2[arrayLength];
  copyArray (keysCopy2, keys, arrayLength);
  copyArray (valuesCopy2, values, arrayLength);
  shellSortKeysAndValues (keysCopy2, valuesCopy2, arrayLength);
  const bool equalKeys = isEqual (keysCopy, keysCopy2, arrayLength);
  TEST_ASSERT( equalKeys );

  if (! keysMayContainDuplicates || expectStableSort) {
    const bool equalValues = isEqual (valuesCopy, valuesCopy2, arrayLength);
    TEST_ASSERT( equalValues );
  }

  // Compare against Tpetra::Details::shellSortKeysAndValues.
  KeyType keysCopy3[arrayLength];
  ValueType valuesCopy3[arrayLength];
  copyArray (keysCopy3, keys, arrayLength);
  copyArray (valuesCopy3, values, arrayLength);
  ::Tpetra::Details::shellSortKeysAndValues (keysCopy3, valuesCopy3, arrayLength);
  const bool equalKeys3 = isEqual (keysCopy2, keysCopy3, arrayLength);
  TEST_ASSERT( equalKeys3 );

  if (! keysMayContainDuplicates || expectStableSort) {
    const bool equalValues3 = isEqual (valuesCopy2, valuesCopy3, arrayLength);
    TEST_ASSERT( equalValues3 );
  }

  // Compare against Tpetra::Details::radixSortKeysAndValues.
  // LSB-first radix sort is always stable and can always handle duplicate keys
  KeyType keysCopy4[arrayLength];
  KeyType keysCopy4Aux[arrayLength];
  ValueType valuesCopy4[arrayLength];
  ValueType valuesCopy4Aux[arrayLength];
  copyArray (keysCopy4, keys, arrayLength);
  copyArray (valuesCopy4, values, arrayLength);
  ::Tpetra::Details::radixSortKeysAndValues<KeyType, ValueType, size_t>(keysCopy4, keysCopy4Aux, valuesCopy4, valuesCopy4Aux, arrayLength, std::numeric_limits<KeyType>::max());
  const bool equalKeys4 = isEqual (keysCopy2, keysCopy4, arrayLength);
  TEST_ASSERT( equalKeys4 );
  const bool equalValues4 = isEqual (valuesCopy2, valuesCopy4, arrayLength);
  TEST_ASSERT( equalValues4 );
}

// Fill the given array with distinct values.
// Kokkos::ArithTraits must have a specialization for ValueType.
template<class ValueType>
void
fillValues (ValueType values[],
            const int arrayLength)
{
  typedef Kokkos::ArithTraits<ValueType> KAT;
  const ValueType ONE = KAT::one ();

  if (arrayLength >= 1) {
    values[0] = ONE;
    for (int k = 1; k < arrayLength; ++k) {
      values[k] = values[k-1] + ONE;
    }
  }
}

// \brief Test shortSortKeysAndValues_${arrayLength} for fixed
//   \c KeyType and \c ValueType types, and for all implemented
//   \c arrayLength values.
//
// \tparam KeyType Type of each entry of the \c keys array.
// \tparam ValueType Type of each entry of the \c values array.
//
// \param keyTypeName [in] Human-readable name of \c KeyType.
// \param valueTypeName [in] Human-readable name of \c ValueType.
template<class KeyType, class ValueType>
void
test_fixedTypes (bool& success,
                 Teuchos::FancyOStream& out,
                 const char keyTypeName[],
                 const char valueTypeName[])
{
  out << "Test shortSortKeysAndValues_*" << " for KeyType="
      << keyTypeName << " and ValueType=" << valueTypeName << endl;
  Teuchos::OSTab tab1 (out);

  const bool expectStableSort = false;

  {
    constexpr int arrayLength = 2;
    out << "Test arrayLength=" << arrayLength << endl;
    Teuchos::OSTab tab2 (out);

    ValueType values[arrayLength];
    fillValues (values, arrayLength);
    {
      KeyType keys[arrayLength] = {0, 1};
      out << "Test in-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 0};
      out << "Test out-of-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 1};
      out << "Test duplicate keys (enforce stable sort)" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
  }

  {
    constexpr int arrayLength = 3;
    out << "Test arrayLength=" << arrayLength << endl;
    Teuchos::OSTab tab2 (out);

    ValueType values[arrayLength];
    fillValues (values, arrayLength);
    {
      KeyType keys[arrayLength] = {1, 3, 5};
      out << "Test in-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 1};
      out << "Test reverse-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {3, 5, 1};
      out << "Test out-of-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3};
      out << "Test in-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 3, 3};
      out << "Test in-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {3, 3, 1};
      out << "Test out-of-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 3};
      out << "Test out-of-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 5};
      out << "Test out-of-order keys with separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
  }

  {
    constexpr int arrayLength = 4;
    out << "Test arrayLength=" << arrayLength << endl;
    Teuchos::OSTab tab2 (out);

    ValueType values[arrayLength];
    fillValues (values, arrayLength);
    {
      KeyType keys[arrayLength] = {1, 3, 5, 7};
      out << "Test in-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {7, 5, 3, 1};
      out << "Test reverse-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {3, 7, 5, 1};
      out << "Test out-of-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 1, 1};
      out << "Test in-order keys, all duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3, 5};
      out << "Test in-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 3, 5, 5};
      out << "Test in-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {5, 5, 3, 1};
      out << "Test out-of-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 1, 1};
      out << "Test out-of-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {5, 3, 1, 5};
      out << "Test out-of-order keys with separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3, 3};
      out << "Test in-order keys with multiple duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {3, 3, 1, 1};
      out << "Test out-of-order keys with multiple duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {3, 1, 3, 1};
      out << "Test out-of-order keys with multiple, separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
  }

  {
    constexpr int arrayLength = 8;
    out << "Test arrayLength=" << arrayLength << endl;
    Teuchos::OSTab tab2 (out);

    ValueType values[arrayLength];
    fillValues (values, arrayLength);
    {
      KeyType keys[arrayLength] = {1, 3, 5, 7, 9, 11, 13, 15};
      out << "Test in-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {15, 13, 11, 9, 7, 5, 3, 1};
      out << "Test reverse-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {15, 11, 3, 7, 5, 1, 13, 9};
      out << "Test out-of-order keys" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = false;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 1, 1, 1, 1, 1, 1};
      out << "Test in-order keys, all duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3, 5, 7, 9, 11, 13};
      out << "Test in-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 3, 5, 7, 9, 11, 13, 13};
      out << "Test in-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {13, 13, 11, 9, 7, 5, 3, 1};
      out << "Test out-of-order keys with duplicates at beginning" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {13, 11, 9, 7, 5, 3, 1, 1};
      out << "Test out-of-order keys with duplicates at end" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {5, 13, 11, 9, 7, 3, 1, 5};
      out << "Test out-of-order keys with separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {1, 1, 3, 3, 5, 5, 7, 7};
      out << "Test in-order keys with multiple duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {7, 7, 5, 5, 3, 3, 1, 1};
      out << "Test out-of-order keys with multiple duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
    {
      KeyType keys[arrayLength] = {3, 1, 3, 7, 5, 1, 9, 11};
      out << "Test out-of-order keys with multiple, separated duplicates" << endl;
      Teuchos::OSTab tab3 (out);
      const bool keysMayContainDuplicates = true;
      test_fixedTypes_fixedArrayLength<KeyType, ValueType, arrayLength> (success, out, keys, values, keysMayContainDuplicates, expectStableSort);
    }
  }
}

TEUCHOS_UNIT_TEST( Sort, shortSortKeysAndValues )
{
  out << "Test shortSortKeysAndValues_*" << endl;
  Teuchos::OSTab tab1 (out);

  test_fixedTypes<short, double> (success, out, "short", "double");
  test_fixedTypes<int, double> (success, out, "int", "double");
  test_fixedTypes<long, double> (success, out, "long", "double");
  test_fixedTypes<long long, double> (success, out, "long long", "double");

  // Experimental design suggests that we don't have to test all
  // possibilities, as long as we get coverage of all independent
  // variables.  Furthermore, ValueType matters very little to the
  // sort routines.  I mainly want to test that they compile for both
  // real and complex ValueType.

  // test_fixedTypes<unsigned short, double> (success, out, "unsigned short", "double");
  // test_fixedTypes<unsigned int, double> (success, out, "unsigned int", "double");
  // test_fixedTypes<unsigned long, double> (success, out, "unsigned long", "double");
  // test_fixedTypes<unsigned long long, double> (success, out, "unsigned long long", "double");

  // test_fixedTypes<short, Kokkos::complex<double> > (success, out, "short", "Kokkos::complex<double>");
  // test_fixedTypes<int, Kokkos::complex<double> > (success, out, "int", "Kokkos::complex<double>");
  // test_fixedTypes<long, Kokkos::complex<double> > (success, out, "long", "Kokkos::complex<double>");
  // test_fixedTypes<long long, Kokkos::complex<double> > (success, out, "long long", "Kokkos::complex<double>");

  test_fixedTypes<unsigned short, Kokkos::complex<double> > (success, out, "unsigned short", "Kokkos::complex<double>");
  test_fixedTypes<unsigned int, Kokkos::complex<double> > (success, out, "unsigned int", "Kokkos::complex<double>");
  test_fixedTypes<unsigned long, Kokkos::complex<double> > (success, out, "unsigned long", "Kokkos::complex<double>");
  test_fixedTypes<unsigned long long, Kokkos::complex<double> > (success, out, "unsigned long long", "Kokkos::complex<double>");

  // test_fixedTypes<short, float> (success, out, "short", "float");
  // test_fixedTypes<int, float> (success, out, "int", "float");
  // test_fixedTypes<long, float> (success, out, "long", "float");
  // test_fixedTypes<long long, float> (success, out, "long long", "float");

  // test_fixedTypes<unsigned short, float> (success, out, "unsigned short", "float");
  // test_fixedTypes<unsigned int, float> (success, out, "unsigned int", "float");
  // test_fixedTypes<unsigned long, float> (success, out, "unsigned long", "float");
  // test_fixedTypes<unsigned long long, float> (success, out, "unsigned long long", "float");

  // test_fixedTypes<short, Kokkos::complex<float> > (success, out, "short", "Kokkos::complex<float>");
  // test_fixedTypes<int, Kokkos::complex<float> > (success, out, "int", "Kokkos::complex<float>");
  // test_fixedTypes<long, Kokkos::complex<float> > (success, out, "long", "Kokkos::complex<float>");
  // test_fixedTypes<long long, Kokkos::complex<float> > (success, out, "long long", "Kokkos::complex<float>");

  // test_fixedTypes<unsigned short, Kokkos::complex<float> > (success, out, "unsigned short", "Kokkos::complex<float>");
  // test_fixedTypes<unsigned int, Kokkos::complex<float> > (success, out, "unsigned int", "Kokkos::complex<float>");
  // test_fixedTypes<unsigned long, Kokkos::complex<float> > (success, out, "unsigned long", "Kokkos::complex<float>");
  // test_fixedTypes<unsigned long long, Kokkos::complex<float> > (success, out, "unsigned long long", "Kokkos::complex<float>");
}

} // namespace (anonymous)

