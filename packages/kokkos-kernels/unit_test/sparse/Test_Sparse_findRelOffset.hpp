/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

//#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Core.hpp"
#include <vector>
#include <iostream>
#include <gtest/gtest.h>
#include "KokkosSparse_findRelOffset.hpp"
#include "KokkosKernels_Utils.hpp"
#ifndef kokkos_complex_double
#define kokkos_complex_double Kokkos::complex<double>
#define kokkos_complex_float Kokkos::complex<float>
#endif

namespace Test{ // (anonymous)
  using std::endl;

  // Test findRelOffset with various array data types and for various cases.
  //
  // This takes the same arguments as if it were declared via the
  // TEUCHOS_UNIT_TEST macro.
  template <typename lno_t, typename DT>
  void generalTest (bool& success, std::ostream &out)
  {
    using KokkosSparse::findRelOffset;
    //typedef int lno_t;
    //typedef Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space> DT;
    typedef Kokkos::View<const lno_t*, DT> IVT;
    typedef Kokkos::View<lno_t*, DT> nIVT;

    //Teuchos::OSTab tab0 (out);
    out << "Test findRelOffset" << endl;
    //Teuchos::OSTab tab1 (out);

    out << "Test empty arrays" << endl;

    // Test with zero indices to search, using a raw array.
    {
      lno_t numEnt = 0;
      IVT nullView;
      const lno_t* indsToSearch = nullView.data();

      const lno_t indToFind = 42;
      for (lno_t hint = 0; hint < 3; ++hint) {
        // Length-zero array is trivially sorted, but try the unsorted
        // case just to make sure that branch of the code is right.
        lno_t offset =
          findRelOffset<lno_t, const lno_t* > (indsToSearch, numEnt,
                                         indToFind, hint, true);

        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array
        offset = findRelOffset<lno_t, const lno_t* > (indsToSearch, numEnt,
                                                indToFind, hint, false);
        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array
      }
    }

    out << "Test the sorted, nonempty array case" << endl;

    // Test the sorted case, with a raw array.
    {
      lno_t numEnt = 7;
      const lno_t indsToSearch[7] = {1, 1, 2, 3, 5, 8, 13};
      nIVT indsToSearch_view("indsToSearch", numEnt);
      typename nIVT::HostMirror h_indsToSearch_view = Kokkos::create_mirror_view (indsToSearch_view);
      for (int i = 0; i < numEnt; ++i) {
        //std::cout << "indsToSearch[i]:" << indsToSearch[i] << std::endl;
        h_indsToSearch_view(i) = indsToSearch[i];
      }
      Kokkos::deep_copy(indsToSearch_view, h_indsToSearch_view);
      Kokkos::fence();
      //KokkosKernels::Impl::kk_print_1Dview(indsToSearch_view);
      const bool isSorted = true;

      for (lno_t hint = 0; hint < 10; ++hint) {
        // Test an index that is not in the array.
        // This one is in [min, max].
        lno_t indNotThere = 4;

        lno_t offset =
          findRelOffset<lno_t, const lno_t* > (indsToSearch_view.data(), numEnt,
                                         indNotThere, hint, isSorted);

        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array

        // Test another index that is not in the array.
        // This one is _not_ in [min, max].
        indNotThere = 42;
        offset = findRelOffset<lno_t, const lno_t* > (indsToSearch_view.data(), numEnt,
                                                indNotThere, hint, isSorted);
        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array

        // Test all indices that are in the array.
        for (lno_t k = 0; k < numEnt; ++k) {

          const lno_t indToFind = indsToSearch[k]; // in the array
          offset = findRelOffset<lno_t, const lno_t* > (indsToSearch_view.data(), numEnt,
                                                  indToFind, hint, isSorted);
          if (indToFind == static_cast<lno_t> (1)) {
            // 1 is a duplicate in this example.  Treat it as a special
            // case.  We don't specify which instance of duplicates the
            // function must return, so either one is fine.

        	ASSERT_TRUE( (  offset == static_cast<lno_t> (0) || offset == static_cast<lno_t> (1) ));
        	/*
            TEST_ASSERT( offset == static_cast<LO> (0) ||
                         offset == static_cast<LO> (1) );
                         */
          }
          else {
        	EXPECT_TRUE( (offset == k ));
            //TEST_EQUALITY( offset, k );
          }
        }
      }
    }

    // Test the sorted case, with a Kokkos::View.
    {
      lno_t numEnt = 7;
      const lno_t indsToSearch[7] = {1, 1, 2, 3, 5, 8, 13};
      nIVT indsToSearch_view("indsToSearch", numEnt);

      typename nIVT::HostMirror h_indsToSearch_view = Kokkos::create_mirror_view (indsToSearch_view);
      for (int i = 0; i < numEnt; ++i) h_indsToSearch_view(i) = indsToSearch[i];
      Kokkos::deep_copy(indsToSearch_view, h_indsToSearch_view);
      Kokkos::fence();

      const bool isSorted = true;

      for (lno_t hint = 0; hint < 10; ++hint) {
        // Test an index that is not in the array.
        // This one is in [min, max].
        lno_t indNotThere = 4;
        lno_t offset = findRelOffset<lno_t, IVT> (indsToSearch_view, numEnt,
                                            indNotThere, hint, isSorted);
        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array

        // Test another index that is not in the array.
        // This one is _not_ in [min, max].
        indNotThere = 42;
        offset = findRelOffset<lno_t, IVT> (indsToSearch_view, numEnt,
                                         indNotThere, hint, isSorted);
        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array

        // Test all indices that are in the array.
        for (lno_t k = 0; k < numEnt; ++k) {
          const lno_t indToFind = indsToSearch[k]; // in the array
          offset = findRelOffset<lno_t, IVT> (indsToSearch_view, numEnt,
                                           indToFind, hint, isSorted);
          if (indToFind == static_cast<lno_t> (1)) {
            // 1 is a duplicate in this example.  Treat it as a special
            // case.  We don't specify which instance of duplicates the
            // function must return, so either one is fine.
          	ASSERT_TRUE( (  offset == static_cast<lno_t> (0) || offset == static_cast<lno_t> (1) ));

            //TEST_ASSERT( offset == static_cast<LO> (0) ||
            //             offset == static_cast<LO> (1) );
          }
          else {
          	EXPECT_TRUE( (offset == k ));
            //TEST_EQUALITY( offset, k );
          }
        }
      }
    }

    out << "Test the unsorted, nonempty array case" << endl;

    // Test the unsorted case, with a raw array.
    {
      lno_t numEnt = 7;
      const lno_t indsToSearch[7] = {8, 1, 13, 1, 3, 2, 5};
      const bool isSorted = false;

      nIVT indsToSearch_view("indsToSearch", numEnt);

      typename nIVT::HostMirror h_indsToSearch_view = Kokkos::create_mirror_view (indsToSearch_view);
      for (int i = 0; i < numEnt; ++i) h_indsToSearch_view(i) = indsToSearch[i];
      Kokkos::deep_copy(indsToSearch_view, h_indsToSearch_view);
      Kokkos::fence();


      for (lno_t hint = 0; hint < 10; ++hint) {
        // Test an index that is not in the array.
        // This one is in [min, max].
        lno_t indNotThere = 4;
        lno_t offset =
          findRelOffset<lno_t, const lno_t* > (indsToSearch_view.data(), numEnt,
                                         indNotThere, hint, isSorted);
        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array

        // Test another index that is not in the array.
        // This one is _not_ in [min, max].
        indNotThere = 42;
        offset = findRelOffset<lno_t, const lno_t* > (indsToSearch_view.data(), numEnt,
                                                indNotThere, hint, isSorted);
        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array

        // Test all indices that are in the array.
        for (lno_t k = 0; k < numEnt; ++k) {
          const lno_t indToFind = indsToSearch[k]; // in the array
          offset = findRelOffset<lno_t, const lno_t* > (indsToSearch_view.data(), numEnt,
                                                  indToFind, hint, isSorted);
          if (indToFind == static_cast<lno_t> (1)) {
            // 1 is a duplicate in this example.  Treat it as a special
            // case.  We don't specify which instance of duplicates the
            // function must return, so either one is fine.
           	ASSERT_TRUE( (  offset == static_cast<lno_t> (1) || offset == static_cast<lno_t> (3) ));
            //TEST_ASSERT( offset == static_cast<LO> (1) ||
            //             offset == static_cast<LO> (3) );
          }
          else {
        	EXPECT_TRUE( (offset == k ));
            //TEST_EQUALITY( offset, k );
          }
        }
      }
    }

    // Test the unsorted case, with a Kokkos::View.
    {
      lno_t numEnt = 7;
      lno_t indsToSearch[7];
      // This assumes UVM.
      indsToSearch[0] = 8;
      indsToSearch[1] = 1;
      indsToSearch[2] = 13;
      indsToSearch[3] = 1;
      indsToSearch[4] = 3;
      indsToSearch[5] = 2;
      indsToSearch[6] = 5;


      nIVT indsToSearch_view("indsToSearch", numEnt);

      typename nIVT::HostMirror h_indsToSearch_view = Kokkos::create_mirror_view (indsToSearch_view);
      for (int i = 0; i < numEnt; ++i) h_indsToSearch_view(i) = indsToSearch[i];
      Kokkos::deep_copy(indsToSearch_view, h_indsToSearch_view);
      Kokkos::fence();


      const bool isSorted = false;

      for (lno_t hint = 0; hint < 10; ++hint) {
        // Test an index that is not in the array.
        // This one is in [min, max].
        lno_t indNotThere = 4;
        lno_t offset = findRelOffset<lno_t, IVT> (indsToSearch_view, numEnt,
                                            indNotThere, hint, isSorted);
        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array

        // Test another index that is not in the array.
        // This one is _not_ in [min, max].
        indNotThere = 42;
        offset = findRelOffset<lno_t, IVT> (indsToSearch_view, numEnt,
                                         indNotThere, hint, isSorted);
        EXPECT_TRUE( (offset == numEnt ));
        //TEST_EQUALITY( offset, numEnt ); // not in the array

        // Test all indices that are in the array.
        for (lno_t k = 0; k < numEnt; ++k) {
          const lno_t indToFind = indsToSearch[k]; // in the array
          offset = findRelOffset<lno_t, IVT> (indsToSearch_view, numEnt,
                                           indToFind, hint, isSorted);
          if (indToFind == static_cast<lno_t> (1)) {
            // 1 is a duplicate in this example.  Treat it as a special
            // case.  We don't specify which instance of duplicates the
            // function must return, so either one is fine.
        	ASSERT_TRUE( (  offset == static_cast<lno_t> (1) || offset == static_cast<lno_t> (3) ));
        	/*
            TEST_ASSERT( offset == static_cast<LO> (1) ||
                         offset == static_cast<LO> (3) );
            */
          }
          else {
            EXPECT_TRUE( (offset == k ));

            //TEST_EQUALITY( offset, k );
          }
        }
      }
    }
  }


  // Test findRelOffset with a longer array.  This ensures that even
  // if findRelOffset optimizes for short arrays by using linear
  // search, we'll still get test coverage for longer arrays.
  //
  // This test doesn't need to exercise all the Kokkos device types.
  // Even if the aforementioned short-array optimization has different
  // constants for different Kokkos device types, a sufficiently long
  // array should exercise all cases.  Thus, this is not a templated
  // test, so we don't need to add it to the list of instantiations
  // for templated tests at the bottom of this file.
  //
  // This takes the same arguments as if it were declared via the
  // TEUCHOS_UNIT_TEST macro.
  template <typename lno_t, typename device_t>
  void testLongArray (bool& success, std::ostream &out)
  {
    using KokkosSparse::findRelOffset;
    //typedef long lno_t; // just for a change

    //Teuchos::OSTab tab0 (out);
    out << "Test findRelOffset with a long array" << endl;
    //Teuchos::OSTab tab1 (out);

    // Start with the array [0, 1, 2, ..., 2n], where the number of
    // entries N = 2n+1 for natural numbers n.  Permute every other
    // entry symmetrically about the middle entry (which exists
    // because the number of entries is odd).  For example, for n = 4:
    // [0 1 2 3 4 5 6 7 8] gets permuted to [8 1 6 3 4 5 2 7 0].  Use
    // this to test findRelOffset.  (We don't just reverse x, in case
    // implementations optimize for reverse contiguous order.)

    const lno_t n = 100;
    const lno_t N = 2*n + 1;
    //std::vector<lno_t> indsToSearch (N);

    typedef Kokkos::View<lno_t*, device_t> lno_view_t;
    lno_view_t indsToSearch("indsToSearch", N);
    typename lno_view_t::HostMirror h_indsToSearch = Kokkos::create_mirror_view (indsToSearch);


    for (lno_t k = 0; k < n; ++k) {
      h_indsToSearch[2*k] = 2*(n - k);
      h_indsToSearch[2*k + 1] = 2*k + 1;
    }
    Kokkos::deep_copy(indsToSearch, h_indsToSearch);
    // We don't need to test all possible hints, just two per search
    // value: the correct hint and some wrong hint.
    for (lno_t k = 0; k < N; ++k) {
      // We use std::vector<LO> in as the template parameter of
      // findRelOffset, because the function should work just fine
      // with anything that acts like a 1-D raw array.
      {
        const lno_t indToFind = indsToSearch[k];
        const lno_t expectedOffset = k;
        const lno_t correctHint = expectedOffset;
        // Add some number not 1 to make the "wrong hint," in case
        // there is a "search nearest" optimization (unlikely -- too
        // many branches).
        const lno_t wrongHint = expectedOffset + 7;

        const lno_t offset0 =
          findRelOffset<lno_t, /*std::vector<lno_t>*/ lno_view_t > (indsToSearch, N, indToFind,
                                               correctHint, false);

        EXPECT_TRUE( (offset0 == expectedOffset ));
        //TEST_EQUALITY( offset0, expectedOffset );
        const lno_t offset1 =
          findRelOffset<lno_t, /*std::vector<lno_t>*/ lno_view_t > (indsToSearch, N, indToFind,
                                               wrongHint, false);
        EXPECT_TRUE( (offset1 == expectedOffset ));
        //TEST_EQUALITY( offset1, expectedOffset );
      }
      {
        // This is the "index not in array" case.  We only need to
        // test one hint here, since all hints are wrong.
        const lno_t indToFind = N + 1; // not in the array
        const lno_t hint = 0;
        const lno_t offset0 =
          findRelOffset<lno_t, /*std::vector<lno_t>*/ lno_view_t > (indsToSearch, N, indToFind,
                                               hint, false);
        EXPECT_TRUE( (offset0 == N ));
        //TEST_EQUALITY( offset0, N );
      }
    }
  }
} // namespace (anonymous)

template<typename lno_t, typename device_t>
void test_findRelOffset()
{
  using namespace Test;

  class NullBuffer : public std::streambuf
  {
  public:
    int overflow(int c) { return c; }
  };
  NullBuffer null_buffer;
  //std::ostream &out = std::cout;
  std::ostream out(&null_buffer);
  out << "Test KokkosSparse::findRelOffset" << endl;

  bool success = true;
  //host test
  generalTest <lno_t, device_t>(success, out);
  EXPECT_TRUE( success);
  //host test
  testLongArray<lno_t, device_t> (success, out);
  EXPECT_TRUE( success);
}



#define EXECUTE_TEST(SCALAR, ORDINAL, OFFSET, DEVICE) \
TEST_F( TestCategory, sparse ## _ ##findRelOffset ## _ ## SCALAR ## _ ## ORDINAL ## _ ## OFFSET ## _ ## DEVICE ) { \
  test_findRelOffset<ORDINAL,DEVICE>(); \
}

#if (defined (KOKKOSKERNELS_INST_ORDINAL_INT) \
  ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int, int, TestExecSpace)
#endif

#if (defined (KOKKOSKERNELS_INST_ORDINAL_INT64_T) \
  ) || (!defined(KOKKOSKERNELS_ETI_ONLY) && !defined(KOKKOSKERNELS_IMPL_CHECK_ETI_CALLS))
 EXECUTE_TEST(double, int64_t, int, TestExecSpace)
#endif


