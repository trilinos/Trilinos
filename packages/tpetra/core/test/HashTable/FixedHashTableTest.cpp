// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Details_FixedHashTable.hpp>
#include <Kokkos_Core.hpp>
#include "Teuchos_OrdinalTraits.hpp"
#include "KokkosCompat_View.hpp"
#include <cstdlib> // atexit
#include <memory>

namespace { // (anonymous)

  // atexit() callback that finalizes the given Kokkos execution
  // space, if it is currently initialized.
  template<class ExecSpace>
  void finalizeExecSpace () {
    if (Kokkos::is_initialized ()) {
      Kokkos::finalize ();
    }
  }

  // Take care of Kokkos execution space initialization automatically.
  // This ensures that each execution space gets initialized and
  // finalized at most once, in that order, over all the tests in this
  // file.  Kokkos requires this, just like MPI does for MPI_Init and
  // MPI_Finalize.
  template<class ExecSpace>
  struct InitExecSpace {
    InitExecSpace () {
      if (! Kokkos::is_initialized ()) {
        Kokkos::initialize ();
      }
      // How should we respond if atexit() fails to register our hook?
      // That means that the Kokkos execution space won't get
      // finalized at the end of the program.  Kokkos already prints
      // out a message in that case.  Thus, we don't have to do
      // anything here.  It's not Kokkos' fault if atexit() ran out of
      // space for all the hooks.
      if (! registeredExitHook_) {
        (void) atexit (finalizeExecSpace<ExecSpace>);
        registeredExitHook_ = true;
      }
    }

    bool isInitialized () {
      return Kokkos::is_initialized ();
    }

    static bool registeredExitHook_;
  };

#ifdef KOKKOS_ENABLE_CUDA
  template<>
  struct InitExecSpace<Kokkos::Cuda> {
    typedef Kokkos::Cuda ExecSpace;

    InitExecSpace () {
      // Make sure that HostSpace's execution space is initialized
      // first.  Otherwise, Cuda::initialize() throws an exception.
      InitExecSpace<Kokkos::HostSpace::execution_space> init2;

      if (! Kokkos::is_initialized ()) {
        Kokkos::initialize ();
      }
      // How should we respond if atexit() fails to register our hook?
      // That means that the Kokkos execution space won't get
      // finalized at the end of the program.  Kokkos already prints
      // out a message in that case.  Thus, we don't have to do
      // anything here.  It's not Kokkos' fault if atexit() ran out of
      // space for all the hooks.
      if (! registeredExitHook_) {
        (void) atexit (finalizeExecSpace<ExecSpace>);
        registeredExitHook_ = true;
      }
    }

    bool isInitialized () {
      return Kokkos::is_initialized ();
    }

    static bool registeredExitHook_;
  };
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_SERIAL
  template<> bool InitExecSpace<Kokkos::Serial>::registeredExitHook_ = false;
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_OPENMP
  template<> bool InitExecSpace<Kokkos::OpenMP>::registeredExitHook_ = false;
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_THREADS
  template<> bool InitExecSpace<Kokkos::Threads>::registeredExitHook_ = false;
#endif // KOKKOS_ENABLE_THREADS

#ifdef KOKKOS_ENABLE_CUDA
  bool InitExecSpace<Kokkos::Cuda>::registeredExitHook_ = false;
#endif // KOKKOS_ENABLE_CUDA


  template<class KeyType, class ValueType, class DeviceType>
  struct TestFixedHashTable {
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    // Fix the layout, so that it's the same on all devices.  This
    // makes testing FixedHashTable's templated copy constructor
    // easier.
    typedef Kokkos::View<const KeyType*, Kokkos::LayoutLeft, DeviceType> keys_type;
    typedef Kokkos::View<const ValueType*, Kokkos::LayoutLeft, DeviceType> vals_type;
    typedef typename keys_type::size_type size_type;

    // For duplicate keys with different values, lookups
    // (FixedHashTable::get) could get different values.  This means
    // that we should only test whether the given keys are in the
    // table; we shouldn't test their values.  testValues = false only
    // tests whether FixedHashTable::get incorrectly returns
    // Teuchos::OrdinalTraits<ValueType>::invalid() for a key that
    // should be in the table; testValues = true also tests whether
    // the returned value matches the expected value.
    static bool
    testKeys (std::ostream& out,
              const table_type& table,
              const keys_type& keys,
              const vals_type& vals,
              const bool testValues = true)
    {
      using std::endl;
      const size_type numKeys = keys.extent (0);

      out << "Test " << numKeys << " key" << (numKeys != 1 ? "s" : "") << ":  ";
      TEUCHOS_TEST_FOR_EXCEPTION
        (numKeys != vals.extent (0), std::logic_error,
         "keys and vals are not the same length!  keys.extent(0) = "
         << numKeys << " != vals.extent(0) = " << vals.extent (0)
         << ".");

      size_type badCount = 0;
      Kokkos::parallel_reduce("testKeys",
        Kokkos::RangePolicy<typename DeviceType::execution_space>(0, numKeys),
        KOKKOS_LAMBDA(const size_t &i, size_type &myBadCount) {
          const KeyType key = keys(i);
          const ValueType expectedVal = vals(i);
          const ValueType actualVal = table.get (key);
          bool curBad = false;
  
          if (actualVal == Tpetra::Details::OrdinalTraits<ValueType>::invalid ()) {
            curBad = true;
            printf("get(key=%lld) = invalid, should have been %d\n",
                    (long long) key, expectedVal);
          }
          else if (testValues && actualVal != expectedVal) {
            curBad = true;
            printf("get(key=%lld) = %d, should have been %d\n",
                    (long long) key, actualVal, expectedVal);
          }
  
          if (curBad) {
            ++myBadCount;
          }
        }, badCount);

      if (badCount == 0) {
        out << "PASSED" << endl;
        return true;
      }
      else {
        out << "FAILED (" << badCount << " out of " << numKeys << ")" << endl;
        return false;
      }
    }
  };

  // Test that an empty table really is empty.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable, Empty, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename Kokkos::View<const KeyType*, Kokkos::LayoutLeft, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test empty table" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }

    const size_type numKeys = 0;
    Kokkos::View<KeyType*, Kokkos::LayoutLeft, DeviceType> keys ("keys", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);

    // Pick something other than 0, just to make sure that it works.
    const ValueType startingValue = 1;

    Teuchos::ArrayView<const KeyType> keys_av = Kokkos::Compat::getArrayView (keys_h);
    out << "Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, startingValue)) );
    if (table.is_null ()) {
      return; // stop the test now to prevent dereferencing null
    }

    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    if (! success) {
      out << "table.hasDuplicateKeys() raised an exception; no sense in "
        "continuing." << endl;
      return;
    }
    TEST_EQUALITY_CONST( duplicateKeys, false );

    const size_type numTestKeys = 4;
    Kokkos::View<KeyType*, Kokkos::LayoutLeft, DeviceType> testKeys ("testKeys", numTestKeys);
    auto testKeys_h = Kokkos::create_mirror_view (testKeys);
    testKeys_h(0) = static_cast<KeyType> (0);
    testKeys_h(1) = static_cast<KeyType> (1);
    testKeys_h(2) = static_cast<KeyType> (-1);
    testKeys_h(3) = static_cast<KeyType> (42);
    Kokkos::deep_copy (testKeys, testKeys_h);

    size_type badCount = 0;
    table_type tt = *table;
    Kokkos::parallel_reduce("testEmpty",
      Kokkos::RangePolicy<typename DeviceType::execution_space>(0, numTestKeys),
      KOKKOS_LAMBDA(const size_t &i, size_type &myBadCount) {
        if (tt.get(testKeys(i)) != Tpetra::Details::OrdinalTraits<ValueType>::invalid()) {
          myBadCount++;
        }
      }, badCount);

    TEST_EQUALITY( badCount, 0);
  }

  // Test contiguous keys, with the constructor that takes a
  // Teuchos::ArrayView of keys (hence the "_T" in the test name) and
  // a single starting value.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_T, ContigKeysStartingValue, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename Kokkos::View<const KeyType*, Kokkos::LayoutLeft, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test contiguous keys, with constructor that takes a single "
      "starting value and produces contiguous values" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }

    const size_type numKeys = 10;
    Kokkos::View<KeyType*, Kokkos::LayoutLeft, DeviceType> keys ("keys", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);
    // Start with some number other than 0, just to make sure that
    // it works.
    for (size_type i = 0; i < numKeys; ++i) {
      keys_h(i) = static_cast<KeyType> (i + 20);
    }
    Kokkos::deep_copy (keys, keys_h);

    // Pick something other than 0, just to make sure that it works.
    const ValueType startingValue = 1;

    // The hash table doesn't need this; we use it only for testing.
    Kokkos::View<ValueType*, Kokkos::LayoutLeft, DeviceType> vals ("vals", numKeys);
    auto vals_h = Kokkos::create_mirror_view (vals);
    for (size_type i = 0; i < numKeys; ++i) {
      vals_h(i) = static_cast<ValueType> (i) + startingValue;
    }
    Kokkos::deep_copy (vals, vals_h);

    Teuchos::ArrayView<const KeyType> keys_av (keys_h.data (), numKeys);
    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, startingValue)) );
    if (table.is_null ()) {
      return; // stop the test now to prevent dereferencing null
    }

    TEST_EQUALITY( keys_h(0), table->minKey () );
    TEST_EQUALITY( keys_h(numKeys-1), table->maxKey () );
    TEST_EQUALITY( startingValue, table->minVal () );
    TEST_EQUALITY( static_cast<ValueType> (startingValue + numKeys - 1), table->maxVal () );
    TEST_EQUALITY( static_cast<size_t> (numKeys), static_cast<size_t> (table->numPairs ()) );

    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    TEST_EQUALITY_CONST( duplicateKeys, false );

    {
      Teuchos::OSTab tab1 (out);
      success = TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table, keys, vals);
      TEST_EQUALITY_CONST( success, true );
    }
  }

  // Test contiguous keys, with the constructor that takes a
  // Kokkos::View of keys (hence the "_K" in the test name) and a
  // single starting value.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_K, ContigKeysStartingValue, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename table_type::keys_type keys_type;
    typedef typename keys_type::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test contiguous keys (as a Kokkos::View), with constructor that "
      "takes a single starting value and produces contiguous values" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }

    const size_type numKeys = 10;
    typename keys_type::non_const_type keys ("keys", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);
    // Start with some number other than 0, just to make sure that
    // it works.
    for (size_type i = 0; i < numKeys; ++i) {
      keys_h(i) = static_cast<KeyType> (i + 20);
    }
    Kokkos::deep_copy (keys, keys_h);

    // Pick something other than 0, just to make sure that it works.
    const ValueType startingValue = 1;

    // The hash table doesn't need this; we use it only for testing.
    Kokkos::View<ValueType*, Kokkos::LayoutLeft, DeviceType> vals ("vals", numKeys);
    auto vals_h = Kokkos::create_mirror_view (vals);
    for (size_type i = 0; i < numKeys; ++i) {
      vals_h(i) = static_cast<ValueType> (i) + startingValue;
    }
    Kokkos::deep_copy (vals, vals_h);

    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys, startingValue)) );
    if (table.is_null ()) {
      return; // stop the test now to prevent dereferencing null
    }

    TEST_EQUALITY( keys_h(0), table->minKey () );
    TEST_EQUALITY( keys_h(numKeys-1), table->maxKey () );
    TEST_EQUALITY( startingValue, table->minVal () );
    TEST_EQUALITY( static_cast<ValueType> (startingValue + numKeys - 1), table->maxVal () );
    TEST_EQUALITY( static_cast<size_t> (numKeys), static_cast<size_t> (table->numPairs ()) );

    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    TEST_EQUALITY_CONST( duplicateKeys, false );

    {
      Teuchos::OSTab tab1 (out);
      success = TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table, keys, vals);
      TEST_EQUALITY_CONST( success, true );
    }
  }

  // Test noncontiguous keys, with the constructor that takes a
  // Teuchos::ArrayView of keys (hence the "_T" in the test name) and
  // a single starting value.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_T, NoncontigKeysStartingValue, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename Kokkos::View<const KeyType*, Kokkos::LayoutLeft, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test noncontiguous keys, with constructor that takes a single "
      "starting value and produces contiguous values" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }

    const size_type numKeys = 5;
    Kokkos::View<KeyType*, Kokkos::LayoutLeft, DeviceType> keys ("keys", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);
    keys_h(0) = static_cast<KeyType> (10);
    keys_h(1) = static_cast<KeyType> (8);
    keys_h(2) = static_cast<KeyType> (12);
    keys_h(3) = static_cast<KeyType> (17);
    keys_h(4) = static_cast<KeyType> (7);
    Kokkos::deep_copy (keys, keys_h);

    // Pick something other than 0, just to make sure that it works.
    const ValueType startingValue = 1;

    // The hash table doesn't need this; we use it only for testing.
    Kokkos::View<ValueType*, Kokkos::LayoutLeft, DeviceType> vals ("vals", numKeys);
    auto vals_h = Kokkos::create_mirror_view (vals);
    for (size_type i = 0; i < numKeys; ++i) {
      vals_h(i) = static_cast<ValueType> (i) + startingValue;
    }
    Kokkos::deep_copy (vals, vals_h);

    Teuchos::ArrayView<const KeyType> keys_av (keys_h.data (), numKeys);
    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, startingValue)) );
    if (table.is_null ()) {
      return; // stop the test now to prevent dereferencing null
    }

    TEST_EQUALITY( keys_h(4), table->minKey () );
    TEST_EQUALITY( keys_h(3), table->maxKey () );
    TEST_EQUALITY( startingValue, table->minVal () );
    TEST_EQUALITY( static_cast<ValueType> (startingValue + numKeys - 1), table->maxVal () );
    TEST_EQUALITY( static_cast<size_t> (numKeys), static_cast<size_t> (table->numPairs ()) );

    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    TEST_EQUALITY_CONST( duplicateKeys, false );

    {
      Teuchos::OSTab tab1 (out);
      success = TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table, keys, vals);
      TEST_EQUALITY_CONST( success, true );
    }
  }

  // Test noncontiguous keys, with the constructor that takes a
  // Kokkos::View of keys (hence the "_K" in the test name) and
  // a single starting value.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_K, NoncontigKeysStartingValue, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename table_type::keys_type keys_type;
    typedef typename keys_type::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test noncontiguous keys (given as a Kokkos::View), with "
      "constructor that takes a single starting value and produces "
      "contiguous values" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }

    const size_type numKeys = 5;
    typename keys_type::non_const_type keys ("keys", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);
    keys_h(0) = static_cast<KeyType> (10);
    keys_h(1) = static_cast<KeyType> (8);
    keys_h(2) = static_cast<KeyType> (12);
    keys_h(3) = static_cast<KeyType> (17);
    keys_h(4) = static_cast<KeyType> (7);
    Kokkos::deep_copy (keys, keys_h);

    // Pick something other than 0, just to make sure that it works.
    const ValueType startingValue = 1;

    // The hash table doesn't need this; we use it only for testing.
    Kokkos::View<ValueType*, Kokkos::LayoutLeft, DeviceType> vals ("vals", numKeys);
    auto vals_h = Kokkos::create_mirror_view (vals);
    for (size_type i = 0; i < numKeys; ++i) {
      vals_h(i) = static_cast<ValueType> (i) + startingValue;
    }
    Kokkos::deep_copy (vals, vals_h);

    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys, startingValue)) );
    if (table.is_null ()) {
      return; // stop the test now to prevent dereferencing null
    }

    TEST_EQUALITY( keys_h(4), table->minKey () );
    TEST_EQUALITY( keys_h(3), table->maxKey () );
    TEST_EQUALITY( startingValue, table->minVal () );
    TEST_EQUALITY( static_cast<ValueType> (startingValue + numKeys - 1), table->maxVal () );
    TEST_EQUALITY( static_cast<size_t> (numKeys), static_cast<size_t> (table->numPairs ()) );

    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    TEST_EQUALITY_CONST( duplicateKeys, false );

    {
      Teuchos::OSTab tab1 (out);
      success = TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table, keys, vals);
      TEST_EQUALITY_CONST( success, true );
    }
  }

  // Test noncontiguous keys, with the constructor that takes a
  // Teuchos::ArrayView of keys and a Teuchos::ArrayView of values.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_T, NoncontigKeysAndVals, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename Kokkos::View<const KeyType*, Kokkos::LayoutLeft, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test noncontiguous keys, with constructor that takes "
      "lists of keys and values" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }

    const size_type numKeys = 5;
    Kokkos::View<KeyType*, Kokkos::LayoutLeft, DeviceType> keys ("keys", numKeys);
    Kokkos::View<ValueType*, Kokkos::LayoutLeft, DeviceType> vals ("vals", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);
    keys_h(0) = static_cast<KeyType> (10);
    keys_h(1) = static_cast<KeyType> (8);
    keys_h(2) = static_cast<KeyType> (12);
    keys_h(3) = static_cast<KeyType> (17);
    keys_h(4) = static_cast<KeyType> (7);
    Kokkos::deep_copy (keys, keys_h);

    // I've chosen the min and max values to occur in different
    // positions than the min and max keys, for maximum generality.
    auto vals_h = Kokkos::create_mirror_view (vals);
    vals_h(0) = static_cast<ValueType> (600);
    vals_h(0) = static_cast<ValueType> (200);
    vals_h(0) = static_cast<ValueType> (500);
    vals_h(0) = static_cast<ValueType> (300);
    vals_h(0) = static_cast<ValueType> (400);
    Kokkos::deep_copy (vals, vals_h);

    Teuchos::ArrayView<const KeyType> keys_av (keys_h.data (), numKeys);
    Teuchos::ArrayView<const ValueType> vals_av (vals_h.data (), numKeys);
    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, vals_av)) );
    if (table.is_null ()) {
      return; // stop the test now to prevent dereferencing null
    }

    TEST_EQUALITY( keys_h(4), table->minKey () );
    TEST_EQUALITY( keys_h(3), table->maxKey () );
    TEST_EQUALITY( vals_h(0), table->minVal () );
    TEST_EQUALITY( vals_h(1), table->maxVal () );
    TEST_EQUALITY( static_cast<size_t> (numKeys), static_cast<size_t> (table->numPairs ()) );

    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    TEST_EQUALITY_CONST( duplicateKeys, false );

    {
      Teuchos::OSTab tab1 (out);
      success = TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table, keys, vals);
      TEST_EQUALITY_CONST( success, true );
    }
  }

  // Test whether FixedHashTable successfully detects duplicate keys.
  // Use the constructor that takes a Teuchos::ArrayView of keys and a
  // single starting value.  (It shouldn't matter which constructor we
  // use, because detection is separate from construction.  However,
  // it would be wise to test both cases.)  The above tests exercise
  // hasDuplicateKeys() for tables that do NOT have duplicate keys,
  // including empty and nonempty tables.  Thus, it suffices here to
  // test a nonempty table with duplicate keys.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_T, DuplicateKeys, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename Kokkos::View<const KeyType*, Kokkos::LayoutLeft, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test detection of duplicate keys" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }

    const size_type numKeys = 6;
    Kokkos::View<KeyType*, Kokkos::LayoutLeft, DeviceType> keys ("keys", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);
    keys_h(0) = static_cast<KeyType> (100);
    keys_h(1) = static_cast<KeyType> (3);
    keys_h(2) = static_cast<KeyType> (50);
    keys_h(3) = static_cast<KeyType> (120);
    keys_h(4) = static_cast<KeyType> (0);
    keys_h(5) = static_cast<KeyType> (100);
    Kokkos::deep_copy (keys, keys_h);

    // Pick something other than 0, just to make sure that it works.
    const ValueType startingValue = 1;

    // The hash table doesn't need this; we use it only for testing.
    Kokkos::View<ValueType*, Kokkos::LayoutLeft, DeviceType> vals ("vals", numKeys);
    auto vals_h = Kokkos::create_mirror_view (vals);
    for (size_type i = 0; i < numKeys; ++i) {
      vals_h(i) = static_cast<ValueType> (i) + startingValue;
    }
    Kokkos::deep_copy (vals, vals_h);

    Teuchos::ArrayView<const KeyType> keys_av (keys_h.data (), numKeys);
    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, startingValue)) );
    if (table.is_null ()) {
      return; // stop the test now to prevent dereferencing null
    }

    // We still require that the table correctly report the min and
    // max keys and values, even if there were duplicate keys.
    TEST_EQUALITY( keys_h(4), table->minKey () );
    TEST_EQUALITY( keys_h(3), table->maxKey () );
    TEST_EQUALITY( startingValue, table->minVal () );
    TEST_EQUALITY( static_cast<ValueType> (startingValue + numKeys - 1), table->maxVal () );
    // The table is supposed to count duplicate keys separately.
    TEST_EQUALITY( static_cast<size_t> (numKeys), static_cast<size_t> (table->numPairs ()) );

    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    // This table actually has duplicate keys.
    TEST_EQUALITY_CONST( duplicateKeys, true );

    // Testing for duplicates should not affect the min and max keys
    // or values.
    TEST_EQUALITY( keys_h(4), table->minKey () );
    TEST_EQUALITY( keys_h(3), table->maxKey () );
    TEST_EQUALITY( startingValue, table->minVal () );
    TEST_EQUALITY( static_cast<ValueType> (startingValue + numKeys - 1), table->maxVal () );
    // The table is supposed to count duplicate keys separately.
    // Asking if the table has duplicate keys must NOT merge those
    // keys.
    TEST_EQUALITY( static_cast<size_t> (numKeys), static_cast<size_t> (table->numPairs ()) );

    // Furthermore, asking for the min and max keys should not affect
    // whether the table reports that it has duplicate keys.
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    // This table actually has duplicate keys.
    TEST_EQUALITY_CONST( duplicateKeys, true );

    // For duplicate keys with different values, lookups
    // (FixedHashTable::get) could get different values.  This means
    // that we should only test whether the given keys are in the
    // table; we shouldn't test their values.
    const bool testValues = false;
    {
      Teuchos::OSTab tab1 (out);
      success =
        TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table,
                                                                      keys, vals,
                                                                      testValues);
      TEST_EQUALITY_CONST( success, true );
    }
  }


  template<class KeyType, class ValueType, class OutDeviceType, class InDeviceType>
  struct TestCopyCtor {
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, OutDeviceType> out_table_type;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, InDeviceType> in_table_type;
    typedef Kokkos::View<const KeyType*, Kokkos::LayoutLeft, InDeviceType> keys_type;
    typedef Kokkos::View<const ValueType*, Kokkos::LayoutLeft, InDeviceType> vals_type;

    static void
    test (std::ostream& out,
          bool& success,
          in_table_type& inTable,
          const keys_type& keys,
          const vals_type& vals,
          const std::string& outDeviceName,
          const std::string& inDeviceName,
          const bool testValues = true)
    {
      using std::endl;

      out << "Test FixedHashTable copy constructor from " << inDeviceName
           << " to " << outDeviceName << endl;
      Teuchos::OSTab tab1 (out);

      // Make sure that the input device's execution space has been
      // initialized.
      TEST_EQUALITY_CONST( Kokkos::is_initialized (), true );
      if (! Kokkos::is_initialized ()) {
        return; // avoid crashes if initialization failed
      }

      // Initialize the output device's execution space, if necessary.
      InitExecSpace<typename OutDeviceType::execution_space> init;
      // This avoids warnings for 'init' being unused.
      TEST_EQUALITY_CONST( init.isInitialized (), true );
      if (! init.isInitialized ()) {
        return; // avoid crashes if initialization failed
      }

      Teuchos::RCP<out_table_type> outTable;
      TEST_NOTHROW( outTable = Teuchos::rcp (new out_table_type (inTable)) );
      if (outTable.is_null ()) {
        return;
      }

      TEST_EQUALITY( inTable.minKey (), outTable->minKey () );
      TEST_EQUALITY( inTable.maxKey (), outTable->maxKey () );
      TEST_EQUALITY( inTable.minVal (), outTable->minVal () );
      TEST_EQUALITY( inTable.maxVal (), outTable->maxVal () );
      TEST_EQUALITY( inTable.numPairs (), outTable->numPairs () );

      // Make sure the new table has duplicate keys if and only if the
      // original table has duplicate keys.
      const bool originalHasDuplicateKeys = inTable.hasDuplicateKeys ();
      const bool newTableHasDuplicateKeys = outTable->hasDuplicateKeys ();
      TEST_EQUALITY( originalHasDuplicateKeys, newTableHasDuplicateKeys );

      // Make sure that computing whether the table has duplicate keys
      // doesn't affect the min or max keys.
      TEST_EQUALITY( inTable.minKey (), outTable->minKey () );
      TEST_EQUALITY( inTable.maxKey (), outTable->maxKey () );
      TEST_EQUALITY( inTable.minVal (), outTable->minVal () );
      TEST_EQUALITY( inTable.maxVal (), outTable->maxVal () );
      TEST_EQUALITY( inTable.numPairs (), outTable->numPairs () );

      Kokkos::View<KeyType*, Kokkos::LayoutLeft, OutDeviceType> keys_out ("keys", keys.extent (0));
      Kokkos::deep_copy (keys_out, keys);
      Kokkos::View<ValueType*, Kokkos::LayoutLeft, OutDeviceType> vals_out ("vals", vals.extent (0));
      Kokkos::deep_copy (vals_out, vals);

      success = TestFixedHashTable<KeyType, ValueType, OutDeviceType>::testKeys (out, *outTable, keys_out, vals_out, testValues);
      TEST_EQUALITY_CONST( success, true );
    }
  };

  // Test FixedHashTable's templated copy constructor, with a list of
  // keys that does not have duplicates.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_T, CopyCtorNoDupKeys, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename Kokkos::View<const KeyType*, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test FixedHashTable's templated copy constructor" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }
    out << "Initialized execution space " << typeid (execution_space).name ()
        << endl;

    // Create the same set of key,value pairs as in the previous test.
    const size_type numKeys = 5;
    Kokkos::View<KeyType*, DeviceType> keys ("keys", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);
    keys_h(0) = static_cast<KeyType> (10);
    keys_h(1) = static_cast<KeyType> (8);
    keys_h(2) = static_cast<KeyType> (12);
    keys_h(3) = static_cast<KeyType> (17);
    keys_h(4) = static_cast<KeyType> (7);
    Kokkos::deep_copy (keys, keys_h);

    // Pick something other than 0, just to make sure that it works.
    const ValueType startingValue = 1;

    // The hash table doesn't need this; we use it only for testing.
    Kokkos::View<ValueType*, DeviceType> vals ("vals", numKeys);
    auto vals_h = Kokkos::create_mirror_view (vals);
    for (size_type i = 0; i < numKeys; ++i) {
      vals_h(i) = static_cast<ValueType> (i) + startingValue;
    }
    Kokkos::deep_copy (vals, vals_h);

    Teuchos::ArrayView<const KeyType> keys_av (keys_h.data (), numKeys);
    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, startingValue)) );
    if (table.is_null ()) {
      return; // fail fast to avoid dereferencing null
    }
    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    // This table does NOT have duplicate keys.
    TEST_EQUALITY_CONST( duplicateKeys, false );

    {
      Teuchos::OSTab tab1 (out);
      success = TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table, keys, vals);
      TEST_EQUALITY_CONST( success, true );
    }

    // Print a warning if only one execution space is enabled in
    // Kokkos, because this means that we can't test the templated
    // copy constructor.
    bool testedAtLeastOnce = false;

#ifdef KOKKOS_ENABLE_SERIAL
    if (! std::is_same<execution_space, Kokkos::Serial>::value) {
      out << "Testing copy constructor to Serial" << endl;
      // The test initializes the output device's execution space if necessary.
      TestCopyCtor<KeyType, ValueType, Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
        DeviceType>::test (out, success, *table, keys, vals,
                           "Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>",
                           typeid(DeviceType).name ());
      testedAtLeastOnce = true;
    }
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_OPENMP
    if (! std::is_same<execution_space, Kokkos::OpenMP>::value) {
      out << "Testing copy constructor to OpenMP" << endl;
      TestCopyCtor<KeyType, ValueType, Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
        DeviceType>::test (out, success, *table, keys, vals,
                           "Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>",
                           typeid(DeviceType).name ());
      testedAtLeastOnce = true;
    }
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_THREADS
    if (! std::is_same<execution_space, Kokkos::Threads>::value) {
      out << "Testing copy constructor to Threads" << endl;
      TestCopyCtor<KeyType, ValueType, Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
        DeviceType>::test (out, success, *table, keys, vals,
                           "Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>",
                           typeid(DeviceType).name ());
      testedAtLeastOnce = true;
    }
#endif // KOKKOS_ENABLE_THREADS

#ifdef KOKKOS_ENABLE_CUDA
    {
      if (! std::is_same<typename DeviceType::memory_space, Kokkos::CudaSpace>::value) {
        out << "Testing copy constructor to CudaSpace" << endl;
        TestCopyCtor<KeyType, ValueType, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
          DeviceType>::test (out, success, *table, keys, vals,
                             "Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>",
                             typeid(DeviceType).name ());
        testedAtLeastOnce = true;
      }
    }
#endif // KOKKOS_ENABLE_CUDA
    if (! testedAtLeastOnce) {
      out << "*** WARNING: Did not actually test FixedHashTable's templated "
        "copy constructor, since only one Kokkos execution space is enabled!"
        " ***" << endl;
    }

    out << "Done testing copy constructor through different devices.  "
      "Now try invoking the FixedHashTable's destructor" << endl;

    TEST_NOTHROW( table = Teuchos::null );

    out << "Got to the end!" << endl;
  }

  // Test FixedHashTable's templated copy constructor, with a list of
  // keys that DOES have duplicates.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_T, CopyCtorDupKeys, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType, DeviceType> table_type;
    typedef typename Kokkos::View<const KeyType*, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test FixedHashTable's templated copy constructor, "
      "with a list of keys that DOES have duplicates" << endl;
    Teuchos::OSTab tab0 (out);

    InitExecSpace<execution_space> init;
    // This avoids warnings for 'init' being unused.
    TEST_EQUALITY_CONST( init.isInitialized (), true );
    if (! init.isInitialized ()) {
      return; // avoid crashes if initialization failed
    }
    out << "Initialized execution space " << typeid (execution_space).name ()
        << endl;

    // Create the same set of key,value pairs as in the previous test.
    const size_type numKeys = 6;
    Kokkos::View<KeyType*, DeviceType> keys ("keys", numKeys);
    auto keys_h = Kokkos::create_mirror_view (keys);
    keys_h(0) = static_cast<KeyType> (10);
    keys_h(1) = static_cast<KeyType> (8);
    keys_h(2) = static_cast<KeyType> (12);
    keys_h(3) = static_cast<KeyType> (17);
    keys_h(4) = static_cast<KeyType> (17);
    keys_h(5) = static_cast<KeyType> (7);
    Kokkos::deep_copy (keys, keys_h);

    // Pick something other than 0, just to make sure that it works.
    const ValueType startingValue = 1;

    // The hash table doesn't need this; we use it only for testing.
    Kokkos::View<ValueType*, DeviceType> vals ("vals", numKeys);
    auto vals_h = Kokkos::create_mirror_view (vals);
    for (size_type i = 0; i < numKeys; ++i) {
      vals_h(i) = static_cast<ValueType> (i) + startingValue;
    }
    Kokkos::deep_copy (vals, vals_h);

    Teuchos::ArrayView<const KeyType> keys_av (keys_h.data (), numKeys);
    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, startingValue)) );
    if (table.is_null ()) {
      return; // fail fast to avoid dereferencing null
    }
    bool duplicateKeys = false;
    TEST_NOTHROW( duplicateKeys = table->hasDuplicateKeys () );
    // This table DOES have duplicate keys.
    TEST_EQUALITY_CONST( duplicateKeys, true );

    // For tables with duplicate keys, we can't promise which value we
    // actually get when we look up one of those keys.
    const bool testValues = false;
    {
      Teuchos::OSTab tab1 (out);
      const bool lclSuccess =
        TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table,
                                                                      keys, vals,
                                                                      testValues);
      TEST_EQUALITY_CONST( lclSuccess, true );
      success = success && lclSuccess;
    }

    // Print a warning if only one execution space is enabled in
    // Kokkos, because this means that we can't test the templated
    // copy constructor.
    bool testedAtLeastOnce = false;

#ifdef KOKKOS_ENABLE_SERIAL
    if (! std::is_same<execution_space, Kokkos::Serial>::value) {
      out << "Testing copy constructor to Serial" << endl;
      // The test initializes the output device's execution space if necessary.
      TestCopyCtor<KeyType, ValueType, Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>,
        DeviceType>::test (out, success, *table, keys, vals,
                           "Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>",
                           typeid (DeviceType).name (), testValues);
      testedAtLeastOnce = true;
    }
#endif // KOKKOS_ENABLE_SERIAL

#ifdef KOKKOS_ENABLE_OPENMP
    if (! std::is_same<execution_space, Kokkos::OpenMP>::value) {
      out << "Testing copy constructor to OpenMP" << endl;
      TestCopyCtor<KeyType, ValueType, Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>,
        DeviceType>::test (out, success, *table, keys, vals,
                           "Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>",
                           typeid (DeviceType).name (), testValues);
      testedAtLeastOnce = true;
    }
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_THREADS
    if (! std::is_same<execution_space, Kokkos::Threads>::value) {
      out << "Testing copy constructor to Threads" << endl;
      TestCopyCtor<KeyType, ValueType, Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>,
        DeviceType>::test (out, success, *table, keys, vals,
                           "Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>",
                           typeid (DeviceType).name (), testValues);
      testedAtLeastOnce = true;
    }
#endif // KOKKOS_ENABLE_THREADS

#ifdef KOKKOS_ENABLE_CUDA
    {
      if (! std::is_same<typename DeviceType::memory_space, Kokkos::CudaSpace>::value) {
        out << "Testing copy constructor to CudaSpace" << endl;
        TestCopyCtor<KeyType, ValueType, Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>,
          DeviceType>::test (out, success, *table, keys, vals,
                             "Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>",
                             typeid(DeviceType).name (), testValues);
        testedAtLeastOnce = true;
      }
    }
#endif // KOKKOS_ENABLE_CUDA
    if (! testedAtLeastOnce) {
      out << "*** WARNING: Did not actually test FixedHashTable's templated "
        "copy constructor, since only one Kokkos execution space is enabled!"
        " ***" << endl;
    }

    out << "Done testing copy constructor through different devices.  "
      "Now try invoking the FixedHashTable's destructor" << endl;

    TEST_NOTHROW( table = Teuchos::null );

    out << "Got to the end!" << endl;
  }

//! @test Check that a default-initialized @ref FixedHashTable claims 0 pairs. This test was proposed in issue #13034.
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_T, ZeroPairsOnDefaultInit, ValueType, KeyType, DeviceType)
  {
    using fixed_hash_table_t = Tpetra::Details::FixedHashTable<KeyType,ValueType,DeviceType>;
    
    fixed_hash_table_t table {};
    
    TEST_EQUALITY_CONST(table.numPairs(), 0);
  }
  
  //
  // Instantiations of the templated unit test(s) above.
  //

  // Magic that makes Tpetra tests work.  (Not _quite_ magic; it
  // relates to macros not liking arguments with commas or spaces in
  // them.  It defines some typedefs to avoid this.)
  TPETRA_ETI_MANGLING_TYPEDEFS()

  // Macro that instantiates all unit tests, templated on all three
  // template parameters.  We use this macro below to instantiate for
  // only the Kokkos execution spaces (DEVICE) that are actually
  // enabled.
#define UNIT_TEST_GROUP_3( LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable, Empty, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_T, ContigKeysStartingValue, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_K, ContigKeysStartingValue, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_T, NoncontigKeysStartingValue, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_K, NoncontigKeysStartingValue, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_T, NoncontigKeysAndVals, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_T, DuplicateKeys, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_T, CopyCtorNoDupKeys, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_T, CopyCtorDupKeys, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_T, ZeroPairsOnDefaultInit, LO, GO, DEVICE )

  // The typedefs below are there because macros don't like arguments
  // with commas in them.

#ifdef KOKKOS_ENABLE_SERIAL
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> serial_device_type;

#define UNIT_TEST_GROUP_SERIAL( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, serial_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_SERIAL )

#else
#  define UNIT_TEST_GROUP_SERIAL( LO, GO )
#endif // KOKKOS_ENABLE_SERIAL


#ifdef KOKKOS_ENABLE_OPENMP
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> openmp_device_type;

#define UNIT_TEST_GROUP_OPENMP( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, openmp_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_OPENMP )

#else
#  define UNIT_TEST_GROUP_OPENMP( LO, GO )
#endif // KOKKOS_ENABLE_OPENMP


#ifdef KOKKOS_ENABLE_THREADS
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> threads_device_type;

#define UNIT_TEST_GROUP_PTHREAD( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, threads_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_PTHREAD )

#else
#  define UNIT_TEST_GROUP_PTHREAD( LO, GO )
#endif // KOKKOS_ENABLE_THREADS


#ifdef KOKKOS_ENABLE_CUDA
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> cuda_device_type;

#define UNIT_TEST_GROUP_CUDA( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, cuda_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_CUDA )

#else
#  define UNIT_TEST_GROUP_CUDA( LO, GO )
#endif // KOKKOS_ENABLE_CUDA

} // namespace (anonymous)
