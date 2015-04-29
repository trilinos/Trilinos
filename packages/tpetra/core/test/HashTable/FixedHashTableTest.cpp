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

#include <Teuchos_UnitTestHarness.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Details_FixedHashTable.hpp>
#include <Kokkos_Core.hpp>

namespace { // (anonymous)

  template<class KeyType, class ValueType, class DeviceType>
  struct TestFixedHashTable {
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType> table_type;
    typedef typename Kokkos::View<const KeyType*, DeviceType>::size_type size_type;

    static bool
    testKeys (std::ostream& out,
              const table_type& table,
              const Kokkos::View<const KeyType*, DeviceType>& keys,
              const Kokkos::View<const ValueType*, DeviceType>& vals)
    {
      using std::endl;
      const size_type numKeys = keys.size ();

      out << "Test " << numKeys << "key" << (numKeys != 1 ? "s" : "") << ":  ";
      TEUCHOS_TEST_FOR_EXCEPTION
        (numKeys != vals.dimension_0 (), std::logic_error,
         "keys and vals are not the same length!  keys.dimension_0() = "
         << numKeys << " != vals.dimension_0() = " << vals.dimension_0 ()
         << ".");

      auto keys_h = Kokkos::create_mirror_view (keys);
      Kokkos::deep_copy (keys_h, keys);
      auto vals_h = Kokkos::create_mirror_view (vals);
      Kokkos::deep_copy (vals_h, vals);
      size_type badCount = 0;
      for (size_type i = 0; i < keys.size (); ++i) {
        const KeyType key = keys_h(i);
        const ValueType expectedVal = vals_h(i);
        const ValueType actualVal = table.get (key);
        if (actualVal != expectedVal) {
          ++badCount;
        }
      }

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

  // Test contiguous keys, with the constructor that takes a
  // Teuchos::ArrayView of keys and a single starting value.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_ArrayView, ContigKeysStartingValue, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType> table_type;
    typedef typename Kokkos::View<const KeyType*, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test contiguous keys, with constructor that takes a single "
      "starting value and produces contiguous values" << endl;
    Teuchos::OSTab tab0 (out);

    bool mustFinalize = false;
    if (! execution_space::is_initialized ()) {
      execution_space::initialize ();
      mustFinalize = true;
    }

    const size_type numKeys = 10;
    Kokkos::View<KeyType*, DeviceType> keys ("keys", numKeys);
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
    Kokkos::View<ValueType*, DeviceType> vals ("vals", numKeys);
    auto vals_h = Kokkos::create_mirror_view (vals);
    for (size_type i = 0; i < numKeys; ++i) {
      vals_h(i) = static_cast<ValueType> (i) + startingValue;
    }
    Kokkos::deep_copy (vals, vals_h);

    Teuchos::ArrayView<const KeyType> keys_av (keys_h.ptr_on_device (), numKeys);
    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, startingValue)) );
    {
      Teuchos::OSTab tab1 (out);
      success = TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table, keys, vals);
      TEST_EQUALITY_CONST( success, true );
    }

    if (mustFinalize) {
      execution_space::finalize ();
    }
  }

  // Test noncontiguous keys, with the constructor that takes a
  // Teuchos::ArrayView of keys and a single starting value.
  //
  // ValueType and KeyType are "backwards" because they correspond to
  // LO resp. GO.  (LO, GO) is the natural order for Tpetra's test
  // macros.
  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(FixedHashTable_ArrayView, NoncontigKeysStartingValue, ValueType, KeyType, DeviceType)
  {
    using std::endl;
    typedef Tpetra::Details::FixedHashTable<KeyType, ValueType> table_type;
    typedef typename Kokkos::View<const KeyType*, DeviceType>::size_type size_type;
    typedef typename DeviceType::execution_space execution_space;

    out << "Test noncontiguous keys, with constructor that takes a single "
      "starting value and produces contiguous values" << endl;
    Teuchos::OSTab tab0 (out);

    bool mustFinalize = false;
    if (! execution_space::is_initialized ()) {
      execution_space::initialize ();
      mustFinalize = true;
    }

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

    Teuchos::ArrayView<const KeyType> keys_av (keys_h.ptr_on_device (), numKeys);
    out << " Create table" << endl;

    Teuchos::RCP<table_type> table;
    TEST_NOTHROW( table = Teuchos::rcp (new table_type (keys_av, startingValue)) );
    {
      Teuchos::OSTab tab1 (out);
      success = TestFixedHashTable<KeyType, ValueType, DeviceType>::testKeys (out, *table, keys, vals);
      TEST_EQUALITY_CONST( success, true );
    }

    if (mustFinalize) {
      execution_space::finalize ();
    }
  }

  //
  // Instantiations of the templated unit test(s) above.
  //

  // Magic that makes Tpetra tests work.  (Not _quite_ magic; it
  // relates to macros not liking arguments with commas or spaces in
  // them.  It defines some typedefs to avoid this.)
  TPETRA_ETI_MANGLING_TYPEDEFS()

  // Set of all unit tests, templated on all three template parameters.
#define UNIT_TEST_GROUP_3( LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_ArrayView, ContigKeysStartingValue, LO, GO, DEVICE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( FixedHashTable_ArrayView, NoncontigKeysStartingValue, LO, GO, DEVICE )

  // The typedefs below are there because macros don't like arguments
  // with commas in them.

#ifdef KOKKOS_HAVE_SERIAL
  typedef Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace> serial_device_type;

#define UNIT_TEST_GROUP_SERIAL( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, serial_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_SERIAL )

#else
#  define UNIT_TEST_GROUP_SERIAL( LO, GO )
#endif // KOKKOS_HAVE_SERIAL


#ifdef KOKKOS_HAVE_OPENMP
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> openmp_device_type;

#define UNIT_TEST_GROUP_OPENMP( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, openmp_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_OPENMP )

#else
#  define UNIT_TEST_GROUP_OPENMP( LO, GO )
#endif // KOKKOS_HAVE_OPENMP


#ifdef KOKKOS_HAVE_PTHREAD
  typedef Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace> threads_device_type;

#define UNIT_TEST_GROUP_PTHREAD( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, threads_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_PTHREAD )

#else
#  define UNIT_TEST_GROUP_PTHREAD( LO, GO )
#endif // KOKKOS_HAVE_PTHREAD


#ifdef KOKKOS_HAVE_CUDA
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> cuda_device_type;

#define UNIT_TEST_GROUP_CUDA( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, cuda_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_CUDA )

#else
#  define UNIT_TEST_GROUP_CUDA( LO, GO )
#endif // KOKKOS_HAVE_CUDA


#ifdef KOKKOS_HAVE_CUDA
  typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace> cuda_uvm_device_type;

#define UNIT_TEST_GROUP_CUDA_UVM( LO, GO ) \
  UNIT_TEST_GROUP_3( LO, GO, cuda_uvm_device_type )

  TPETRA_INSTANTIATE_LG( UNIT_TEST_GROUP_CUDA_UVM )

#else
#  define UNIT_TEST_GROUP_CUDA_UVM( LO, GO )
#endif // KOKKOS_HAVE_CUDA

} // namespace (anonymous)
