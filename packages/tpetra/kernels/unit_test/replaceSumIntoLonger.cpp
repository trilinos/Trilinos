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

#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Sparse_CrsMatrix.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <sstream>

namespace { // (anonymous)
  using std::endl;

  template<class CrsMatrixType, const int numEntToModify>
  class ModifyEntries {
  public:
    typedef typename CrsMatrixType::ordinal_type ordinal_type;
    typedef typename CrsMatrixType::value_type scalar_type;

    // The type of the reduction result.
    typedef ordinal_type value_type;

    ModifyEntries (const CrsMatrixType& A,
                   const bool replace,
                   const bool sorted,
                   const bool atomic) :
      A_ (A), replace_ (replace), sorted_ (sorted)
    {}

    KOKKOS_FUNCTION void
    operator () (const ordinal_type& lclRow, ordinal_type& numModified) const
    {
      typedef Kokkos::Details::ArithTraits<scalar_type> KAT;
      typedef typename KAT::mag_type mag_type;
      const scalar_type ONE = KAT::one ();

      const ordinal_type ncol =
        A_.numCols () < static_cast<ordinal_type> (numEntToModify) ?
        A_.numCols () :
        static_cast<ordinal_type> (numEntToModify);
      ordinal_type cols[numEntToModify];
      scalar_type vals[numEntToModify];

      // Indices A.numCols() - 1, A.numCols() - 2, ..., 0 always exist
      // in the row, given how we construct the matrix.  We put them
      // in reverse order here, just to check that the method works if
      // the input (not the matrix's row) is not sorted.

      for (ordinal_type k = 0; k < ncol; ++k) {
        const ordinal_type colToChange = A_.numCols () - k - static_cast<ordinal_type> (1);
        // Cast integers to mag_type first, since direct cast from
        // (e.g.,) int to Kokkos::complex (or std::complex) doesn't
        // work.
        const scalar_type curVal = ONE * static_cast<mag_type> (colToChange + 1);

        cols[k] = colToChange;
        // The expected result will always be -curVal, whether we're
        // doing replace or sumInto.  This lets us make sure that we
        // changed the right value.
        vals[k] = replace_ ? -curVal : -(curVal + curVal);

        // std::cout << "A.numCols(): " << A_.numCols ()
        //           << ", ncol: " << ncol
        //           << ", k: " << k
        //           << ", cols[k]: " << cols[k]
        //           << ", vals[k]: " << vals[k]
        //           << ", original value: " << A_.values[A_.graph.row_map[0] + k]
        //           << std::endl;
      }

      ordinal_type lclNumModified = 0;
      if (replace_) {
        lclNumModified = A_.replaceValues (lclRow, cols, ncol, vals, sorted_, atomic_);
      }
      else { // sumInto
        lclNumModified = A_.sumIntoValues (lclRow, cols, ncol, vals, sorted_, atomic_);
      }
      numModified += lclNumModified;
    }

  private:
    CrsMatrixType A_;
    bool replace_;
    bool sorted_;
    bool atomic_;
  };

  template<class CrsMatrixType, const int numEntToModify>
  void
  modifyEntries (bool& success,
                 Teuchos::FancyOStream& outRef, // see notes
                 const CrsMatrixType& A,
                 const bool replace,
                 const bool sorted,
                 const bool atomic,
                 const bool debug = false)
  {
    using Teuchos::RCP;
    typedef typename CrsMatrixType::device_type::execution_space execution_space;
    typedef typename CrsMatrixType::ordinal_type ordinal_type;
    typedef Kokkos::RangePolicy<execution_space, ordinal_type> policy_type;
    typedef ModifyEntries<CrsMatrixType, numEntToModify> functor_type;

    // If debug is false, we capture all output in an
    // std::ostringstream, and don't print it unless the test fails
    // inside here.  If debug is true, we print all output
    // immediately.
    //
    // Teuchos' unit test macros write to the Teuchos::FancyOStream&
    // 'out'.  They expect the variable to have that name.  That's why
    // the input argument to this function has a different name -- so
    // we can replace it here.
    RCP<std::ostringstream> dbgOutPtr;
    if (! debug) {
      dbgOutPtr = Teuchos::rcp (new std::ostringstream ());
    }
    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::rcpFromRef (outRef) :
      Teuchos::getFancyOStream (dbgOutPtr);
    Teuchos::FancyOStream& out = *outPtr;

    functor_type functor (A, replace, sorted, atomic);
    ordinal_type numModified = 0;
    policy_type range (0, A.numRows ());
    Kokkos::parallel_reduce (std::string ("ModifyEntries"), range, functor, numModified);

    const ordinal_type numEntShouldModify =
      static_cast<ordinal_type> (numEntToModify) <= A.numCols () ?
      static_cast<ordinal_type> (numEntToModify) :
      A.numCols ();
    TEST_EQUALITY( numModified, numEntShouldModify );

    if (! success && ! debug) {
      outRef << dbgOutPtr->str ();
    }
  }

  template<class CrsMatrixType, const int numEntToModify>
  void
  checkWhetherEntriesWereModified (bool& success,
                                   Teuchos::FancyOStream& outRef, // see notes
                                   const CrsMatrixType& A,
                                   const bool replace,
                                   const bool /* sorted */,
                                   const bool /* atomic */,
                                   const bool debug = false)
  {
    using Teuchos::RCP;
    typedef typename CrsMatrixType::value_type value_type;
    typedef typename CrsMatrixType::ordinal_type ordinal_type;
    typedef Kokkos::Details::ArithTraits<value_type> KAT;

    // If debug is false, we capture all output in an
    // std::ostringstream, and don't print it unless the test fails
    // inside here.  If debug is true, we print all output
    // immediately.
    //
    // Teuchos' unit test macros write to the Teuchos::FancyOStream&
    // 'out'.  They expect the variable to have that name.  That's why
    // the input argument to this function has a different name -- so
    // we can replace it here.
    RCP<std::ostringstream> dbgOutPtr;
    if (! debug) {
      dbgOutPtr = Teuchos::rcp (new std::ostringstream ());
    }
    RCP<Teuchos::FancyOStream> outPtr = debug ?
      Teuchos::rcpFromRef (outRef) :
      Teuchos::getFancyOStream (dbgOutPtr);
    Teuchos::FancyOStream& out = *outPtr;

    const value_type ONE = KAT::one ();
    const ordinal_type ncol =
      A.numCols () < static_cast<ordinal_type> (numEntToModify) ?
      A.numCols () :
      static_cast<ordinal_type> (numEntToModify);
    // modifyEntries changes entries with column indices N-1, N-2,
    // ..., max(N - numEntToModify, 0), where N = A.numCols().  Make
    // sure that the "lower bound" works for signed or unsigned
    // ordinal_type.
    const ordinal_type lowerBound = A.numCols () - ncol;

    Teuchos::OSTab tab0 (out);
    out << "check: "
        << "{numCols: " << A.numCols ()
        << ", numEntToModify: " << numEntToModify
        << ", ncol: " << ncol
        << ", lowerBound: " << lowerBound
        << "}" << endl;
    Teuchos::OSTab tab1 (out);

    auto val_h = Kokkos::create_mirror_view (A.values);
    Kokkos::deep_copy (val_h, A.values);
    auto ind_h = Kokkos::create_mirror_view (A.graph.entries);
    Kokkos::deep_copy (ind_h, A.graph.entries);

    const ordinal_type numRows = A.numRows ();
    TEST_EQUALITY( numRows, static_cast<ordinal_type> (1) );
    if (numRows != static_cast<ordinal_type> (1)) {
      return; // stop the test early
    }

    value_type curVal = ONE;
    for (ordinal_type k = 0; k < A.numCols (); ++k, curVal += ONE) {
      ordinal_type curCol;
      value_type expectedVal;

      // Cast integers to mag_type first before assigning to
      // value_type, since direct cast from (e.g.,) int to
      // Kokkos::complex (or std::complex) doesn't work.

      if (ind_h(k) < lowerBound) { // entry should not have been modified
        out << "ind_h(" << k << ") = " << ind_h(k) << "; entry should not have been modified" << endl;
        expectedVal = curVal;
      }
      else {
        out << "ind_h(" << k << ") = " << ind_h(k) << "; entry should have been modified" << endl;
        // The expected result for modified entries will always be
        // -curVal, whether we're doing replace or sumInto.  This lets
        // us make sure that we changed the right value.
        expectedVal = -curVal;
      }

      if (val_h(k) == expectedVal) {
        out << "CORRECT" << endl;
      }
      else {
        success = false;
        out << "ERROR: k: " << k << ", ind_h(k): " << ind_h(k) << ", "
            << "val_h(" << k << ") = " << val_h(k) << " != " << expectedVal
            << " (lowerBound = " << lowerBound << ")" << endl;
      }
    }

    if (! success && ! debug) {
      outRef << dbgOutPtr->str ();
    }
  }

  template<class CrsMatrixType, const int numEntriesToModify>
  void
  testOneCaseImpl (bool& success,
                   Teuchos::FancyOStream& out,
                   const CrsMatrixType& A,
                   const bool replace,
                   const bool sorted,
                   const bool atomic,
                   const bool debug = false)
  {
    typedef typename CrsMatrixType::value_type scalar_type;
    typedef typename CrsMatrixType::ordinal_type ordinal_type;

    if (A.numCols () >= static_cast<ordinal_type> (numEntriesToModify)) {
      Teuchos::OSTab tab0 (out);
      out << "numEntriesToModify: " << numEntriesToModify << endl;
      bool lclSuccess = true;
      modifyEntries<CrsMatrixType, numEntriesToModify> (lclSuccess, out, A, replace, sorted, atomic, debug);
      // If modifyEntries didn't work, no need to test further.
      if (lclSuccess) {
        checkWhetherEntriesWereModified<CrsMatrixType, numEntriesToModify> (lclSuccess, out, A, replace, sorted, atomic, debug);
        TEST_ASSERT( lclSuccess ); // this modifies 'success' and prints to 'out'
      }

      // Restore original values.
      auto val_h = Kokkos::create_mirror_view (A.values);
      const scalar_type ONE = Kokkos::Details::ArithTraits<scalar_type>::one ();
      scalar_type curVal = ONE;
      for (ordinal_type k = 0; k < A.numCols (); ++k, curVal += ONE) {
        val_h[k] = curVal;
      }
      Kokkos::deep_copy (A.values, val_h);
    }
  }


  template<class CrsMatrixType, const int numEntriesToModify>
  struct TestOneCase {
    static void
    test (bool& success,
          Teuchos::FancyOStream& out,
          const CrsMatrixType& A,
          const bool replace,
          const bool sorted,
          const bool atomic,
          const bool debug = false)
    {
      testOneCaseImpl<CrsMatrixType, numEntriesToModify> (success, out, A, replace, sorted, atomic, debug);
      if (! success) {
        return; // Don't bother continuing
      }
      // Yay template recursion!
      TestOneCase<CrsMatrixType, numEntriesToModify / 2 >::test (success, out, A, replace, sorted, atomic, debug);
    }
  };

  // Base case of template recursion for numEntriesToModify = 1.
  template<class CrsMatrixType>
  struct TestOneCase<CrsMatrixType, 1> {
    static void
    test (bool& success,
          Teuchos::FancyOStream& out,
          const CrsMatrixType& A,
          const bool replace,
          const bool sorted,
          const bool atomic,
          const bool debug = false)
    {
      constexpr int numEntriesToModify = 1;
      testOneCaseImpl<CrsMatrixType, numEntriesToModify> (success, out, A, replace, sorted, atomic, debug);
      // This is the base case, so don't recurse on numEntriesToModify.
    }
  };

  // We might never call this second base case of numEntriesToModify =
  // 0, but it ensures that the template recursion always terminates,
  // even if maxNumEntriesToModify (see below) is 0.
  template<class CrsMatrixType>
  struct TestOneCase<CrsMatrixType, 0> {
    static void
    test (bool& /* success */,
          Teuchos::FancyOStream& /* out */,
          const CrsMatrixType& /* A */,
          const bool /* replace */,
          const bool /* sorted */,
          const bool /* atomic */,
          const bool /* debug */ )
    {}
  };

  template<class CrsMatrixType>
  void
  testOneCase (bool& success,
               Teuchos::FancyOStream& out,
               const CrsMatrixType& A,
               const bool replace,
               const bool sorted,
               const bool atomic,
               const bool debug = false)
  {
    Teuchos::OSTab tab0 (out);
    out << "replace: " << (replace ? "true" : "false")
        << ", sorted: " << (sorted ? "true" : "false")
        << ", atomic: " << (atomic ? "true" : "false")
        << endl;
    Teuchos::OSTab tab1 (out);

    constexpr int maxNumEntriesToModify = 128;
    // Invoke template recursion.
    TestOneCase<CrsMatrixType, maxNumEntriesToModify>::test (success, out, A, replace, sorted, atomic, debug);
  }

  template<class CrsMatrixType>
  void
  testOneSize (bool& success,
               Teuchos::FancyOStream& out,
               const CrsMatrixType& A,
               const bool debug = false)
  {
    Teuchos::OSTab tab0 (out);
    out << "testOneSize: {numRows: " << A.numRows ()
        << ", numCols: " << A.numCols () << "}" << endl;

    for (int replaceInt = 0; replaceInt < 2; ++replaceInt) {
      const bool replace = replaceInt != 0;
      for (int sortedInt = 0; sortedInt < 2; ++sortedInt) {
        const bool sorted = sortedInt != 0;
        for (int atomicInt = 0; atomicInt < 2; ++atomicInt) {
          const bool atomic = atomicInt != 0;
          testOneCase (success, out, A, replace, sorted, atomic, debug);

          if (! success) {
            return; // Don't bother continuing
          }
        }
      }
    }
  }

  // Test KokkosSparse::CrsMatrix::{replace,sumInto}Values* with rows
  // of various lengths.  This ensures that we catch any issues with
  // special cases of these methods for short or long rows.
  template<class CrsMatrixType>
  void
  testAllSizes (bool& success,
                Teuchos::FancyOStream& out,
                const typename CrsMatrixType::size_type maxNumEnt,
                const bool debug = false)
  {
    typedef CrsMatrixType matrix_type;
    typedef typename matrix_type::value_type value_type;
    typedef typename matrix_type::ordinal_type ordinal_type;
    typedef typename matrix_type::size_type size_type;
    const value_type ONE = Kokkos::Details::ArithTraits<value_type>::one ();

    Teuchos::OSTab tab0 (out);
    out << "maxNumEnt: " << maxNumEnt << endl;
    Teuchos::OSTab tab1 (out);

    // This directory already has a test (replaceSumInto.cpp) for
    // matrices with more than one row.  Thus, for this test, we can
    // make a matrix with just one row, but vary the length of that
    // row.  This also lets us reuse storage.

    const ordinal_type numRows = 1;

    typename matrix_type::row_map_type::non_const_type ptr ("ptr", numRows+1);
    auto ptr_h = Kokkos::create_mirror_view (ptr);
    typename matrix_type::index_type::non_const_type ind_whole ("ind", maxNumEnt);
    auto ind_whole_h = Kokkos::create_mirror_view (ind_whole);
    typename matrix_type::values_type::non_const_type val_whole ("val", maxNumEnt);
    auto val_whole_h = Kokkos::create_mirror_view (val_whole);

    for (size_type numEnt = 1; numEnt <= maxNumEnt; numEnt *= static_cast<size_type> (2)) {
      const ordinal_type numCols = numEnt;

      out << "Test " << numRows << " x " << numCols << " matrix with " << numEnt
          << " entr" << (numEnt != static_cast<size_type> (1) ? "ies" : "y") << endl;

      ptr_h[0] = 0;
      ptr_h[1] = numEnt;
      Kokkos::deep_copy (ptr, ptr_h);

      const auto range = Kokkos::make_pair (static_cast<size_type> (0), numEnt);
      auto ind = Kokkos::subview (ind_whole, range);
      auto val = Kokkos::subview (val_whole, range);
      auto ind_h = Kokkos::subview (ind_whole_h, range);
      auto val_h = Kokkos::subview (val_whole_h, range);

      value_type curVal = ONE;
      for (size_type k = 0; k < numEnt; ++k, curVal += ONE) {
        ind_h[k] = static_cast<ordinal_type> (k); // sorted entries
        val_h[k] = curVal;
      }
      Kokkos::deep_copy (ind, ind_h);
      Kokkos::deep_copy (val, val_h);

      matrix_type A ("A", numRows, numCols, numEnt, val, ptr, ind);
      testOneSize (success, out, A, debug);
      if (! success) {
        return; // Don't bother continuing
      }
    }
  }

  // The first two arguments let us call Teuchos unit test macros
  // inside.  Those macros expect 'success' and 'out' to have exactly
  // those names.
  void
  generalTest (bool& success,
               Teuchos::FancyOStream& out,
               const bool debug = false)
  {
    typedef double SC;
    typedef int LO;
    typedef Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space> DT;
    typedef KokkosSparse::CrsMatrix<SC, LO, DT> matrix_type;

    Teuchos::OSTab tab0 (out);
    out << "Test KokkosSparse::CrsMatrix::{replace,sumInto}Values*" << endl;
    Teuchos::OSTab tab1 (out);

    const matrix_type::size_type maxNumEnt = 1024;
    testAllSizes<matrix_type> (success, out, maxNumEnt, debug);
  }

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using std::endl;

  Teuchos::RCP<Teuchos::FancyOStream> outPtr =
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  Teuchos::FancyOStream& out = *outPtr;

  out << "Call Kokkos::initialize" << endl;
  Kokkos::initialize (argc, argv);

  bool success = true;
  out << "Run test" << endl;
  const bool debug = false;
  generalTest (success, out, debug);

  out << "Call Kokkos::finalize" << endl;
  Kokkos::finalize ();

  if (success) {
    out << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}


