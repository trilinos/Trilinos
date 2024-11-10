// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_copyOffsets.hpp"
#include <limits>

namespace { // (anonymous)

  using Kokkos::View;
  using Kokkos::view_alloc;
  using Kokkos::WithoutInitializing;

  template<class OutputType, class InputType>
  void
  testNoOverflow (bool& success, Teuchos::FancyOStream& out)
  {
    const size_t len (6);
    View<InputType*> src (view_alloc ("src", WithoutInitializing), len);
    View<OutputType*> dst (view_alloc ("dst", WithoutInitializing), len);

    auto src_h = Kokkos::create_mirror_view (src);
    src_h[0] = InputType (0);
    src_h[1] = InputType (1);
    src_h[2] = InputType (3);
    src_h[3] = InputType (5);
    src_h[4] = InputType (7);
    src_h[5] = InputType (9);
    Kokkos::deep_copy (src, src_h);

    Tpetra::Details::copyOffsets (dst, src);

    auto dst_h = Kokkos::create_mirror_view (dst);
    Kokkos::deep_copy (dst_h, dst);

    TEUCHOS_ASSERT( dst_h[0] == OutputType (0) );
    TEUCHOS_ASSERT( dst_h[1] == OutputType (1) );
    TEUCHOS_ASSERT( dst_h[2] == OutputType (3) );
    TEUCHOS_ASSERT( dst_h[3] == OutputType (5) );
    TEUCHOS_ASSERT( dst_h[4] == OutputType (7) );
    TEUCHOS_ASSERT( dst_h[5] == OutputType (9) );
  }

  // "Convert" between two signed types that are the same

  TEUCHOS_UNIT_TEST( CopyOffsets, int_to_int )
  {
    testNoOverflow<int, int> (success, out);
  }

  // "Convert" between two unsigned types that are the same

  TEUCHOS_UNIT_TEST( CopyOffsets, unsigned_int_to_unsigned_int )
  {
    testNoOverflow<unsigned int, unsigned int> (success, out);
  }

  // Convert between smaller and larger types, both signed

  TEUCHOS_UNIT_TEST( CopyOffsets, int_to_long_long )
  {
    testNoOverflow<long long, int> (success, out);
  }

  TEUCHOS_UNIT_TEST( CopyOffsets, long_long_to_int )
  {
    testNoOverflow<int, long long> (success, out);
  }

  // Convert between smaller and larger types, both unsigned

  TEUCHOS_UNIT_TEST( CopyOffsets, unsigned_int_to_unsigned_long_long )
  {
    testNoOverflow<unsigned long long, unsigned int> (success, out);
  }

  TEUCHOS_UNIT_TEST( CopyOffsets, unsigned_long_long_to_unsigned_int )
  {
    testNoOverflow<unsigned int, unsigned long long> (success, out);
  }

  // Convert between smaller signed type and larger unsigned type

  TEUCHOS_UNIT_TEST( CopyOffsets, size_t_to_int )
  {
    testNoOverflow<int, size_t> (success, out);
  }

  TEUCHOS_UNIT_TEST( CopyOffsets, int_to_size_t )
  {
    testNoOverflow<size_t, int> (success, out);
  }

  // Convert between smaller unsigned type and larger signed type

  TEUCHOS_UNIT_TEST( CopyOffsets, unsigned_int_to_long_long )
  {
    testNoOverflow<long long, unsigned int> (success, out);
  }

  TEUCHOS_UNIT_TEST( CopyOffsets, long_long_to_unsigned_int )
  {
    testNoOverflow<unsigned int, long long> (success, out);
  }

  TEUCHOS_UNIT_TEST( CopyOffsets, TestOverflow )
  {
    using InputType = size_t;
    using OutputType = int;

    const bool debug = Tpetra::Details::Behavior::debug ();
    if (debug) {
      const size_t len (6);
      View<InputType*> src (view_alloc ("src", WithoutInitializing), len);
      View<OutputType*> dst (view_alloc ("dst", WithoutInitializing), len);

      auto src_h = Kokkos::create_mirror_view (src);
      src_h[0] = InputType (0);
      src_h[1] = InputType (1);
      src_h[2] = InputType (3);
      src_h[3] = InputType (5);
      src_h[4] = std::numeric_limits<InputType>::max (); // will overflow
      src_h[5] = InputType (9);
      Kokkos::deep_copy (src, src_h);

      TEST_THROW( Tpetra::Details::copyOffsets (dst, src),
                  std::runtime_error );
    }
  }

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraSession (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
