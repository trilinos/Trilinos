// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>

namespace { // (anonymous)

  using std::endl;

  template<class OffsetType, class CountType>
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
    const CountType count = 3;

    Kokkos::View<OffsetType*> offsets ("offsets", numCounts+1);
    // The initial contents shouldn't matter, so fill offsets with
    // an "invalid" flag value (-1 for signed types).
    Kokkos::deep_copy (offsets, Tpetra::Details::OrdinalTraits<OffsetType>::invalid ());

    using ::Tpetra::Details::computeOffsetsFromConstantCount;
    TEST_NOTHROW( computeOffsetsFromConstantCount (offsets, count) );

    auto offsets_h = Kokkos::create_mirror_view (offsets);
    Kokkos::deep_copy (offsets_h, offsets);

    TEST_EQUALITY( offsets_h(0), static_cast<OffsetType> (0) );
    for (CountType k = 0; k < numCounts; ++k) {
      // Test result against sequential computation
      TEST_EQUALITY( offsets_h(k+1), offsets_h(k) + count );
      // Test against closed-form formula for partial sums
      TEST_EQUALITY( offsets_h(k+1),
                     static_cast<OffsetType> (count) *
                     static_cast<OffsetType> (k+1) );
      // Another sanity check
      TEST_EQUALITY( static_cast<CountType> (offsets_h(k+1) - offsets_h(k)), count );
    }

    // In release mode, only print at the end of this function if
    // something fails.
    if (! debug && ! success) {
      originalOutputStream << releaseOutputStream->str ();
    }
  }

  void
  testComputeOffsets (bool& success,
                      Teuchos::FancyOStream& out,
                      const bool debug)
  {
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
      testComputeOffsetsTmpl<int, int> (success, out,
                                        "int",
                                        countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<long, int> (success, out,
                                         "long",
                                         countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<long long, int> (success, out,
                                              "long long",
                                              countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<unsigned int, int> (success, out,
                                                 "unsigned int",
                                                 countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<unsigned long, int> (success, out,
                                                  "unsigned long",
                                                  countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<unsigned long long, int> (success, out,
                                                       "unsigned long long",
                                                       countTypeName, debug);
      if (! success) {
        return;
      }
    }
    {
      const char countTypeName[] = "unsigned int";
      testComputeOffsetsTmpl<unsigned int, unsigned int> (success, out,
                                                          "unsigned int",
                                                          countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<unsigned long, unsigned int> (success, out,
                                                           "unsigned long",
                                                           countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<unsigned long long, unsigned int> (success, out,
                                                                "unsigned long long",
                                                                countTypeName, debug);
      if (! success) {
        return;
      }
    }
    {
      const char countTypeName[] = "long";
      testComputeOffsetsTmpl<long, long> (success, out,
                                          "long",
                                          countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<long long, long> (success, out,
                                               "long long",
                                               countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<unsigned long, long> (success, out,
                                                   "unsigned long",
                                                   countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<unsigned long long, long> (success, out,
                                                        "unsigned long long",
                                                        countTypeName, debug);
    }
    {
      const char countTypeName[] = "unsigned long";
      testComputeOffsetsTmpl<unsigned long, unsigned long> (success, out,
                                                            "unsigned long",
                                                            countTypeName, debug);
      if (! success) {
        return;
      }
      // We can't test OffsetType = long long here, in case
      // sizeof(long) == sizeof(long long).
      testComputeOffsetsTmpl<unsigned long long, unsigned long> (success, out,
                                                                 "unsigned long long",
                                                                 countTypeName, debug);
      if (! success) {
        return;
      }
    }
    {
      const char countTypeName[] = "long long";
      testComputeOffsetsTmpl<long long, long long> (success, out,
                                                    "long long",
                                                    countTypeName, debug);
      if (! success) {
        return;
      }
      testComputeOffsetsTmpl<unsigned long long, long long> (success, out,
                                                             "unsigned long long",
                                                             countTypeName, debug);
      if (! success) {
        return;
      }
    }
    {
      const char countTypeName[] = "unsigned long long";
      testComputeOffsetsTmpl<unsigned long long, unsigned long long> (success, out,
                                                                      "unsigned long long",
                                                                      countTypeName, debug);
    }
  }

} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  Kokkos::initialize (argc, argv);

  const bool debug = true;
  Teuchos::RCP<Teuchos::FancyOStream> outPtr =
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (debug ? std::cerr : std::cout));
  bool success = true;
  testComputeOffsets (success, *outPtr, debug);

  Kokkos::finalize ();

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
