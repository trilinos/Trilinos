// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>
#include <type_traits>

namespace TpetraTest {
  // CUDA 7.5 doesn't always like functors in anonymous namespaces,
  // so we put this one in a named namespace.
  template<class OffsetType, class CountType, class DeviceType>
  class SumOfCounts {
  public:
    SumOfCounts (const Kokkos::View<const CountType*, DeviceType>& counts) :
      counts_ (counts) {}
    KOKKOS_INLINE_FUNCTION void
    operator () (const CountType& k, OffsetType& inout) const {
      inout += counts_(k);
    }
  private:
    Kokkos::View<const CountType*, DeviceType> counts_;
  };

  // Compute sum of counts.  We use this for error checking.
  template<class OffsetType, class CountType, class DeviceType>
  OffsetType sumOfCounts (const Kokkos::View<const CountType*, DeviceType>& counts) {
    typedef typename Kokkos::View<const CountType*, DeviceType>::size_type size_type;
    typedef Kokkos::RangePolicy<typename DeviceType::execution_space, size_type> range_type;
    typedef SumOfCounts<OffsetType, CountType, DeviceType> functor_type;

    OffsetType total = 0;
    range_type range (0, counts.extent (0));
    Kokkos::parallel_reduce (range, functor_type (counts), total);
    return total;
  }
} // namespace TpetraTest


namespace { // (anonymous)

  using std::endl;

  template<class ExecutionSpace>
  struct ExecSpaceName {
    static const char* name () { return ExecutionSpace().name(); }
  };

  template<class MemorySpace>
  struct MemorySpaceName {
    static const char* name () { return MemorySpace().name(); }
  };

  template<class DeviceType>
  std::string deviceName ()
  {
    return std::string ("Kokkos::Device<") +
      DeviceType::execution_space::name() +
      std::string (", ") +
      DeviceType::memory_space::name() +
      std::string (" >");
  }

  template<class OffsetType, class CountType, class DeviceType>
  void
  testComputeOffsetsTmpl (bool& success,
                          Teuchos::FancyOStream& originalOutputStream,
                          const char offsetTypeName[],
                          const char countTypeName[],
                          const bool debug)
  {
    static_assert (std::is_integral<OffsetType>::value,
                   "OffsetType must be a built-in integer type.");
    static_assert (std::is_integral<CountType>::value,
                   "CountType must be a built-in integer type.");

    // In debug mode, print output right away.  In release mode, only
    // print at the end of this function if something fails.
    Teuchos::RCP<Teuchos::FancyOStream> outPtr;
    Teuchos::RCP<std::ostringstream> releaseOutputStream;
    if (debug) {
      outPtr = Teuchos::rcpFromRef (originalOutputStream);
    }
    else {
      releaseOutputStream = Teuchos::rcp (new std::ostringstream ());
      outPtr = Teuchos::getFancyOStream (releaseOutputStream);
    }
    // The Teuchos unit test macros assume 'out' exists in their scope.
    Teuchos::FancyOStream& out = *outPtr;

    Teuchos::OSTab tab0 (out);
    out << "Test OffsetType = " << offsetTypeName
        << ", CountType = " << countTypeName << endl;
    Teuchos::OSTab tab1 (out);

    const CountType numCounts = 10;

    // Set up counts array.  Fill it with entries, all of which are
    // different, but whose running sums we can easily calculate.
    Kokkos::View<CountType*, DeviceType> counts ("counts", numCounts);
    auto counts_h = Kokkos::create_mirror_view (counts);
    for (CountType k = 0; k < numCounts; ++k) {
      counts_h(k) = k + 1;
    }
    Kokkos::deep_copy (counts, counts_h);

    const OffsetType ZERO = 0;
    const OffsetType ONE = 1;
    const OffsetType TWO = 2;

    // Make sure that our sum formula is correct.
    const OffsetType expectedTotal = TpetraTest::sumOfCounts<OffsetType, CountType, DeviceType> (counts);
    TEST_EQUALITY( expectedTotal, (numCounts*(numCounts+ONE)) / TWO );

    Kokkos::View<OffsetType*, DeviceType> offsets ("offsets", numCounts+1);
    // The initial contents shouldn't matter, so fill offsets with
    // an "invalid" flag value (-1 for signed types).
    Kokkos::deep_copy (offsets, Tpetra::Details::OrdinalTraits<OffsetType>::invalid ());

    using ::Tpetra::Details::computeOffsetsFromCounts;

    if (std::is_same<OffsetType, CountType>::value) {
      out << "Test the case where counts and offsets alias one another" << endl;
      Teuchos::OSTab tab2 (out);
      using Kokkos::subview;
      typedef Kokkos::pair<size_t, size_t> range_type;
      auto counts_in = subview (offsets, range_type (0, counts.extent (0)));
      Kokkos::deep_copy (counts_in, counts);

      OffsetType computedTotal = 0;
      TEST_NOTHROW( computedTotal = computeOffsetsFromCounts (offsets, counts_in) );
      TEST_EQUALITY( expectedTotal, computedTotal );

      auto offsets_h = Kokkos::create_mirror_view (offsets);
      Kokkos::deep_copy (offsets_h, offsets);

      TEST_EQUALITY( offsets_h(0), ZERO );
      for (CountType k = 0; k < numCounts; ++k) {
        // Test result against sequential computation
        TEST_EQUALITY( offsets_h(k+1), offsets_h(k) + counts_h(k) );
        // Test against closed-form formula for partial sums
        TEST_EQUALITY( offsets_h(k+1), ((k + ONE)*(k + TWO)) / TWO );
        // Another sanity check
        TEST_EQUALITY( static_cast<CountType> (offsets_h(k+1) - offsets_h(k)), counts_h(k) );
      }

      if (! success) {
        out << "Test FAILED; returning early" << endl;
        return;
      }
    }

    out << "Test the case where counts and offsets do not alias one another, "
      "but live in the same memory space" << endl;
    {
      Teuchos::OSTab tab2 (out);
      OffsetType computedTotal = 0;
      TEST_NOTHROW( computedTotal = computeOffsetsFromCounts (offsets, counts) );
      TEST_EQUALITY( expectedTotal, computedTotal );

      // Make sure that computeOffsetsFromCounts didn't change counts.
      {
        const OffsetType total = TpetraTest::sumOfCounts<OffsetType, CountType, DeviceType> (counts);
        TEST_EQUALITY( total, expectedTotal );
      }

      auto offsets_h = Kokkos::create_mirror_view (offsets);
      Kokkos::deep_copy (offsets_h, offsets);

      TEST_EQUALITY( offsets_h(0), ZERO );
      for (CountType k = 0; k < numCounts; ++k) {
        // Test result against sequential computation
        TEST_EQUALITY( offsets_h(k+1), offsets_h(k) + counts_h(k) );
        // Test against closed-form formula for partial sums
        TEST_EQUALITY( offsets_h(k+1), ((k + ONE)*(k + TWO)) / TWO );
        // Another sanity check
        TEST_EQUALITY( static_cast<CountType> (offsets_h(k+1) - offsets_h(k)), counts_h(k) );
      }
    }

    // Now test the case where counts lives in host memory, and
    // offsets in device memory.  We only need to test this if device
    // != host.
    if (! std::is_same<typename DeviceType::memory_space, Kokkos::HostSpace>::value) {
      out << "Test the case where counts lives in host memory and "
        "offsets in device memory" << endl;
      Teuchos::OSTab tab2 (out);

      // Kokkos guarantees that if an execution space is initialized,
      // its HostMirror's execution space is also initialized.
      typedef typename Kokkos::View<CountType*, DeviceType>::HostMirror::execution_space host_execution_space;
      // If DeviceType::memory_space is CudaUVMSpace, HostMirror's
      // memory space may also be CudaUVMSpace (since this is
      // accessible from host).  We really want the count array to
      // live in host memory for this test.
      typedef Kokkos::HostSpace host_memory_space;
      typedef Kokkos::Device<host_execution_space, host_memory_space> host_device_type;

      Kokkos::View<CountType*, host_device_type> counts_host ("counts", numCounts);
      for (CountType k = 0; k < numCounts; ++k) {
        counts_host(k) = k + 1;
      }

      OffsetType computedTotal = 0;
      TEST_NOTHROW( computedTotal = computeOffsetsFromCounts (offsets, counts_host) );
      TEST_EQUALITY( expectedTotal, computedTotal );

      // Make sure that computeOffsetsFromCounts didn't change counts_host.
      {
        const OffsetType total = TpetraTest::sumOfCounts<OffsetType, CountType, host_device_type> (counts_host);
        TEST_EQUALITY( total, expectedTotal );
      }

      auto offsets_h = Kokkos::create_mirror_view (offsets);
      Kokkos::deep_copy (offsets_h, offsets);

      TEST_EQUALITY( offsets_h(0), ZERO );
      for (CountType k = 0; k < numCounts; ++k) {
        // Test result against sequential computation
        TEST_EQUALITY( offsets_h(k+1), offsets_h(k) + counts_host(k) );
        // Test against closed-form formula for partial sums
        TEST_EQUALITY( offsets_h(k+1), ((k + ONE)*(k + TWO)) / TWO );
        // Another sanity check
        TEST_EQUALITY( static_cast<CountType> (offsets_h(k+1) - offsets_h(k)), counts_host(k) );
      }
    }

    // In release mode, only print at the end of this function if
    // something fails.
    if (! debug && ! success) {
      originalOutputStream << releaseOutputStream->str ();
    }
  }


  template<class DeviceType>
  void
  testComputeOffsets (bool& success,
                      Teuchos::FancyOStream& out,
                      const bool debug)
  {
    Teuchos::OSTab tab0 (out);
    out << "Test DeviceType = " << deviceName<DeviceType> () << endl;

    // OffsetType must be able to hold sums of CountType values.
    // Thus, OffsetType can't be smaller than CountType, and
    // OffsetType can't be signed if CountType is unsigned.
    //
    // sizeof(long) == sizeof(int) on some platforms (e.g., Windows),
    // so we must not test (e.g.,) CountType = unsigned int and
    // OffsetType = long.
    //
    // The typical non-GPU case is OffsetType = size_t, CountType =
    // size_t or int.  We don't use size_t explicitly here because
    // this is a typedef with platform-dependent type.

    {
      const char countTypeName[] = "int";
      testComputeOffsetsTmpl<int, int, DeviceType> (success, out,
                                                    "int",
                                                    countTypeName, debug);
      if (! success) {
        out << "Test with " << "int, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<long, int, DeviceType> (success, out,
                                                     "long",
                                                     countTypeName, debug);
      if (! success) {
        out << "Test with " << "long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<long long, int, DeviceType> (success, out,
                                                          "long long",
                                                          countTypeName, debug);
      if (! success) {
        out << "Test with " << "long long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<unsigned int, int, DeviceType> (success, out,
                                                             "unsigned int",
                                                             countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned int, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<unsigned long, int, DeviceType> (success, out,
                                                              "unsigned long",
                                                              countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<unsigned long long, int, DeviceType> (success, out,
                                                                   "unsigned long long",
                                                                   countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
    }
    {
      const char countTypeName[] = "unsigned int";
      testComputeOffsetsTmpl<unsigned int, unsigned int, DeviceType> (success, out,
                                                                      "unsigned int",
                                                                      countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned int, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<unsigned long, unsigned int, DeviceType> (success, out,
                                                                       "unsigned long",
                                                                       countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<unsigned long long, unsigned int, DeviceType> (success, out,
                                                                            "unsigned long long",
                                                                            countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
    }
    {
      const char countTypeName[] = "long";
      testComputeOffsetsTmpl<long, long, DeviceType> (success, out,
                                                      "long",
                                                      countTypeName, debug);
      if (! success) {
        out << "Test with " << "long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<long long, long, DeviceType> (success, out,
                                                           "long long",
                                                           countTypeName, debug);
      if (! success) {
        out << "Test with " << "long long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<unsigned long, long, DeviceType> (success, out,
                                                               "unsigned long",
                                                               countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<unsigned long long, long, DeviceType> (success, out,
                                                                    "unsigned long long",
                                                                    countTypeName, debug);
    }
    {
      const char countTypeName[] = "unsigned long";
      testComputeOffsetsTmpl<unsigned long, unsigned long, DeviceType> (success, out,
                                                                        "unsigned long",
                                                                        countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      // We can't test OffsetType = long long here, in case
      // sizeof(long) == sizeof(long long).
      testComputeOffsetsTmpl<unsigned long long, unsigned long, DeviceType> (success, out,
                                                                             "unsigned long long",
                                                                             countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
    }
    {
      const char countTypeName[] = "long long";
      testComputeOffsetsTmpl<long long, long long, DeviceType> (success, out,
                                                                "long long",
                                                                countTypeName, debug);
      if (! success) {
        out << "Test with " << "long long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
      testComputeOffsetsTmpl<unsigned long long, long long, DeviceType> (success, out,
                                                                         "unsigned long long",
                                                                         countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
    }
    {
      const char countTypeName[] = "unsigned long long";
      testComputeOffsetsTmpl<unsigned long long, unsigned long long, DeviceType> (success, out,
                                                                                  "unsigned long long",
                                                                                  countTypeName, debug);
      if (! success) {
        out << "Test with " << "unsigned long long, " << countTypeName << ", "
            << deviceName<DeviceType> () << " failed; returning early" << endl;
        return;
      }
    }
  }

  void
  runTests (bool& success,
            Teuchos::FancyOStream& out,
            const bool debug)
  {
    // There is no DefaultMemorySpace typedef in Kokkos.
    typedef Kokkos::Device<Kokkos::DefaultExecutionSpace,
      Kokkos::DefaultExecutionSpace::memory_space> device_type;

    // Test the default execution space, and the execution space of
    // its host mirror.  Kokkos::initialize always initializes both.
    testComputeOffsets<device_type> (success, out, debug);
    typedef Kokkos::View<int*, device_type>::HostMirror::device_type host_device_type;
    // The host mirror may be the same; don't run the test again in
    // that case.  It's still correct to compile the test again in
    // that case.
    if (! std::is_same<device_type, host_device_type>::value) {
      testComputeOffsets<host_device_type> (success, out, debug);

      // Host mirror execution space of Cuda may be Serial (or OpenMP
      // or Threads), but host mirror memory space of CudaUVMSpace
      // (which may or may not be the default memory space of Cuda) is
      // CudaUVMSpace, not HostSpace.  Thus, we need to test HostSpace
      // too.
      if (! std::is_same<host_device_type::memory_space, Kokkos::HostSpace>::value) {
        typedef Kokkos::Device<host_device_type::execution_space, Kokkos::HostSpace> cur_device_type;
        testComputeOffsets<cur_device_type> (success, out, debug);
      }
    }

#ifdef KOKKOS_ENABLE_CUDA
    {
      // Make sure that we test both without and with UVM.
      // We only have to test once for each case.
      using Kokkos::Cuda;
      using Kokkos::CudaSpace;
      using Kokkos::CudaUVMSpace;
      using mem_space = typename device_type::memory_space;
      if (! std::is_same<mem_space, CudaSpace>::value) {
        using cur_device_type = Kokkos::Device<Cuda, CudaSpace>;
        testComputeOffsets<cur_device_type> (success, out, debug);
      }
      if (! std::is_same<mem_space, CudaUVMSpace>::value) {
        using cur_device_type = Kokkos::Device<Cuda, CudaUVMSpace>;
        testComputeOffsets<cur_device_type> (success, out, debug);
      }
    }
#endif // KOKKOS_ENABLE_CUDA
  }
} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  const bool throwExceptionsOnParseError = true;
  // Teuchos doesn't know about Kokkos' command-line options, so tell
  // it to ignore options it doesn't recognize.
  const bool reportErrorOnUnrecognizedOption = false;
  Teuchos::CommandLineProcessor clp (throwExceptionsOnParseError,
                                     reportErrorOnUnrecognizedOption);
  bool debug = false;
  clp.setOption ("debug", "release", &debug, "Set debug mode.  In debug mode, "
                 "print all output.  In release mode, most output only prints "
                 "if the test fails.  Debug mode is useful if the test crashes "
                 "before it detects failure.");
  const auto clpResult = clp.parse (argc, argv);
  if (clpResult == Teuchos::CommandLineProcessor::PARSE_ERROR) {
    std::cout << "Failed to parse command-line arguments!" << endl
              << "End Result: TEST FAILED" << endl;
  }

  bool success = true;
  {
    Kokkos::ScopeGuard kokkosScope (argc, argv);
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    RCP<FancyOStream> outPtr =
      getFancyOStream (rcpFromRef (debug ? std::cerr : std::cout));
    runTests (success, *outPtr, debug);
  }

  // The Teuchos unit test framework needs to see this to figure out
  // whether the test passed.
  if (success) {
    std::cout << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}
