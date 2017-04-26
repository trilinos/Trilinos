/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <openmp/TestOpenMP.hpp>

namespace Test {

TEST_F( openmp, atomics )
{
  const int loop_count = 1e4;

  ASSERT_TRUE( ( TestAtomic::Loop< int, Kokkos::OpenMP >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< int, Kokkos::OpenMP >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< int, Kokkos::OpenMP >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< unsigned int, Kokkos::OpenMP >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< unsigned int, Kokkos::OpenMP >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< unsigned int, Kokkos::OpenMP >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< long int, Kokkos::OpenMP >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< long int, Kokkos::OpenMP >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< long int, Kokkos::OpenMP >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< unsigned long int, Kokkos::OpenMP >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< unsigned long int, Kokkos::OpenMP >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< unsigned long int, Kokkos::OpenMP >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< long long int, Kokkos::OpenMP >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< long long int, Kokkos::OpenMP >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< long long int, Kokkos::OpenMP >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< double, Kokkos::OpenMP >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< double, Kokkos::OpenMP >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< double, Kokkos::OpenMP >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< float, Kokkos::OpenMP >( 100, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< float, Kokkos::OpenMP >( 100, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< float, Kokkos::OpenMP >( 100, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< Kokkos::complex<double>, Kokkos::OpenMP >( 100, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< Kokkos::complex<double>, Kokkos::OpenMP >( 100, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< Kokkos::complex<double>, Kokkos::OpenMP >( 100, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< TestAtomic::SuperScalar<4>, Kokkos::OpenMP >( 100, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< TestAtomic::SuperScalar<4>, Kokkos::OpenMP >( 100, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< TestAtomic::SuperScalar<4>, Kokkos::OpenMP >( 100, 3 ) ) );
}

TEST_F( openmp, atomic_operations )
{
  const int start = 1; // Avoid zero for division.
  const int end = 11;

  for ( int i = start; i < end; ++i )
  {
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< int, Kokkos::OpenMP >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned int, Kokkos::OpenMP >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long int, Kokkos::OpenMP >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< unsigned long int, Kokkos::OpenMP >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 8 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 9 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 11 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestIntegralType< long long int, Kokkos::OpenMP >( start, end - i, 12 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< double, Kokkos::OpenMP >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< double, Kokkos::OpenMP >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< double, Kokkos::OpenMP >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< double, Kokkos::OpenMP >( start, end - i, 4 ) ) );

    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< float, Kokkos::OpenMP >( start, end - i, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< float, Kokkos::OpenMP >( start, end - i, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< float, Kokkos::OpenMP >( start, end - i, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicOperations::AtomicOperationsTestNonIntegralType< float, Kokkos::OpenMP >( start, end - i, 4 ) ) );
  }
}

TEST_F( openmp, atomic_views_integral )
{
  const long length = 1000000;
  {
    // Integral Types.
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestIntegralType< long, Kokkos::OpenMP >( length, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestIntegralType< long, Kokkos::OpenMP >( length, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestIntegralType< long, Kokkos::OpenMP >( length, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestIntegralType< long, Kokkos::OpenMP >( length, 4 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestIntegralType< long, Kokkos::OpenMP >( length, 5 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestIntegralType< long, Kokkos::OpenMP >( length, 6 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestIntegralType< long, Kokkos::OpenMP >( length, 7 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestIntegralType< long, Kokkos::OpenMP >( length, 8 ) ) );
  }
}

TEST_F( openmp, atomic_views_nonintegral )
{
  const long length = 1000000;
  {
    // Non-Integral Types.
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestNonIntegralType<  double, Kokkos::OpenMP >( length, 1 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestNonIntegralType<  double, Kokkos::OpenMP >( length, 2 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestNonIntegralType<  double, Kokkos::OpenMP >( length, 3 ) ) );
    ASSERT_TRUE( ( TestAtomicViews::AtomicViewsTestNonIntegralType<  double, Kokkos::OpenMP >( length, 4 ) ) );
  }
}

TEST_F( openmp, atomic_view_api )
{
  TestAtomicViews::TestAtomicViewAPI<int, Kokkos::OpenMP >();
}

} // namespace Test
