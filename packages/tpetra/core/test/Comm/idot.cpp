// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_idot.hpp"
#ifdef HAVE_TPETRACORE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_TPETRACORE_MPI
#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"

namespace { // (anonymous)

using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;

typedef Tpetra::Map<> map_type;
typedef Tpetra::MultiVector<> mv_type;
typedef Tpetra::Vector<> vec_type;
typedef vec_type::scalar_type SC;
typedef vec_type::mag_type mag_type;
typedef map_type::local_ordinal_type LO;
typedef map_type::global_ordinal_type GO;
typedef vec_type::device_type device_type;
typedef Teuchos::ScalarTraits<SC> STS;

/// \brief Test Tpetra::idot.
///
/// \param out [out] Output stream; valid (writeable) only on Process
///   0 in the input communicator.
/// \param comm [in] Communicator over which to do the test.
void
testIdot (bool& success,
          std::ostream& out,
          const RCP<const Comm<int> >& comm)
{
  const SC ZERO = STS::zero ();
  const SC ONE = STS::one ();
  const SC TWO = ONE + ONE;
  const SC THREE = TWO + ONE;
  const int numProcs = comm->getSize ();

  // lclSuccess: Local success status.  0 means a failure happened on
  //   the calling process.
  // gblSuccess [in/out] Global success status.  0 means a failure
  //   happened on some process in the input communicator.
  int lclSuccess = 1; // to be updated below
  int gblSuccess = 0; // output argument (see below)

  const LO lclNumRows = 10;
  const GO gblNumRows =
    static_cast<GO> (lclNumRows) * static_cast<GO> (numProcs);
  const GO indexBase = 0;
  RCP<const map_type> map =
    rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

  out << "Test Tpetra::Vector inputs and raw pointer output" << endl;
  {
    vec_type x (map);
    vec_type y (map);
    const SC valX = TWO;
    const SC valY = THREE;
    x.putScalar (valX);
    y.putScalar (valY);

    SC result = ZERO;
    out << "About to call idot" << endl;
    auto req = Tpetra::idot (&result, x, y);
    out << "Finished calling idot" << endl;
    req->wait ();
    out << "Finished wait" << endl;
    const SC N = static_cast<SC> (static_cast<mag_type> (gblNumRows));
    const SC expectedResult = N * valX * valY;

    TEST_EQUALITY( expectedResult, result );

    lclSuccess = success ? 1 : 0; // input argument
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    success = (gblSuccess != 0);
  }

  out << "Test Tpetra::Vector inputs and rank-0 Kokkos::View output" << endl;
  {
    vec_type x (map);
    vec_type y (map);
    const SC valX = TWO;
    const SC valY = THREE;
    x.putScalar (valX);
    y.putScalar (valY);

    Kokkos::View<vec_type::dot_type, vec_type::device_type> result ("result");
    auto result_h = Kokkos::create_mirror_view (result);
    result_h() = ZERO;
    Kokkos::deep_copy (result, result_h);

    auto req = Tpetra::idot (result, x, y);
    req->wait ();
    const SC N = static_cast<SC> (static_cast<mag_type> (gblNumRows));
    const SC expectedResult = N * valX * valY;

    Kokkos::deep_copy (result_h, result);

    TEST_EQUALITY( expectedResult, result_h() );

    lclSuccess = success ? 1 : 0; // input argument
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    success = (gblSuccess != 0);
  }

  out << "Test Tpetra::MultiVector inputs and rank-1 Kokkos::View output" << endl;
  {
    constexpr size_t numVecs = 3;
    mv_type x (map, numVecs);
    mv_type y (map, numVecs);
    const SC valX = TWO;
    const SC valY = THREE;
    x.putScalar (valX);
    y.putScalar (valY);

    Kokkos::View<SC*, device_type> results ("results[numVecs]", numVecs);
    auto results_h = Kokkos::create_mirror_view (results);
    for (size_t k = 0; k < numVecs; ++k) {
      results_h(k) = ZERO;
    }
    Kokkos::deep_copy (results, results_h);

    auto req = Tpetra::idot (results, x, y);
    req->wait ();
    Kokkos::deep_copy (results_h, results);
    const SC N = static_cast<SC> (static_cast<mag_type> (gblNumRows));
    const SC expectedResult = N * valX * valY;

    for (size_t k = 0; k < numVecs; ++k) {
      TEST_EQUALITY( expectedResult, results_h(k) );
    }

    lclSuccess = success ? 1 : 0; // input argument
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    success = (gblSuccess != 0);

    out << "Test special case of (multiple columns) dot (single column)" << endl;
    auto y_0 = y.getVector (0);
    req = Tpetra::idot (results, x, *y_0);
    req->wait ();
    Kokkos::deep_copy (results_h, results);
    for (size_t k = 0; k < numVecs; ++k) {
      TEST_EQUALITY( expectedResult, results_h(k) );
    }

    lclSuccess = success ? 1 : 0; // input argument
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    success = (gblSuccess != 0);

    out << "Test special case of (single column) dot (multiple columns)" << endl;
    auto x_0 = x.getVector (0);
    req = Tpetra::idot (results, *x_0, y);
    req->wait ();
    Kokkos::deep_copy (results_h, results);
    for (size_t k = 0; k < numVecs; ++k) {
      TEST_EQUALITY( expectedResult, results_h(k) );
    }

    lclSuccess = success ? 1 : 0; // input argument
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    success = (gblSuccess != 0);
  }

  out << "Test noncontiguous Tpetra::MultiVector inputs and raw pointer output" << endl;
  {
    out << " First, test contiguous Tpetra::MultiVector inputs "
      "and raw pointer output" << endl;
    constexpr size_t origNumVecs = 5;
    mv_type X (map, origNumVecs);
    mv_type Y (map, origNumVecs);

    SC valX = ONE;
    SC valY = TWO;
    for (size_t j = 0; j < origNumVecs; ++j, valX += ONE, valY += ONE) {
      X.getVectorNonConst (j)->putScalar (valX);
      Y.getVectorNonConst (j)->putScalar (valY);
    }

    SC origResults[origNumVecs];
    for (size_t k = 0; k < origNumVecs; ++k) {
      origResults[k] = ZERO;
    }
    auto req = Tpetra::idot (origResults, X, Y);
    req->wait ();

    // Print results all to a single string first, then to the output
    // stream, to make results more coherent across processes.
    {
      std::ostringstream os;
      os << "  Results: [";
      for (size_t j = 0; j < origNumVecs; ++j) {
        os << origResults[j];
        if (j + 1 < origNumVecs) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();
    }

    const SC N = static_cast<SC> (static_cast<mag_type> (gblNumRows));
    valX = ONE;
    valY = TWO;
    for (size_t j = 0; j < origNumVecs; ++j, valX += ONE, valY += ONE) {
      const SC expectedResult = N * valX * valY;
      TEST_EQUALITY( expectedResult, origResults[j] );
    }

    lclSuccess = success ? 1 : 0; // input argument
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    success = (gblSuccess != 0);

    out << " Now, test noncontiguous Tpetra::MultiVector inputs "
      "and raw pointer output" << endl;
    constexpr size_t newNumVecs = 2;
    Teuchos::Array<size_t> whichCols (newNumVecs);
    // For maximum generality, pick noncontiguous columns, and don't
    // use the zeroth column.
    whichCols[0] = 1;
    whichCols[1] = 3;
    RCP<mv_type> X_sub = X.subViewNonConst (whichCols ());
    RCP<mv_type> Y_sub = Y.subViewNonConst (whichCols ());

    SC newResults[newNumVecs];
    for (size_t k = 0; k < newNumVecs; ++k) {
      newResults[k] = ZERO;
    }
    req = Tpetra::idot (newResults, *X_sub, *Y_sub);
    req->wait ();

    // Print results all to a single string first, then to the output
    // stream, to make results more coherent across processes.
    {
      std::ostringstream os;
      os << "  Results: [";
      for (size_t j = 0; j < newNumVecs; ++j) {
        os << newResults[j];
        if (j + 1 < newNumVecs) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();
    }

    for (size_t k = 0; k < newNumVecs; ++k) {
      TEST_EQUALITY( origResults[whichCols[k]], newResults[k] );
    }

    lclSuccess = success ? 1 : 0; // input argument
    gblSuccess = 0; // output argument
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    success = (gblSuccess != 0);
  }
}


TEUCHOS_UNIT_TEST( idot, basic )
{
  out << "Testing Tpetra::idot" << endl;
  Teuchos::OSTab tab1 (out);

#ifdef HAVE_TPETRACORE_MPI
  RCP<const Comm<int> > comm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
#else
  RCP<const Comm<int> > comm = rcp (new Teuchos::SerialComm<int> ());
#endif // HAVE_TPETRACORE_MPI

  testIdot (success, std::cerr, comm);
  // Just to make sure that if the test fails, we still print
  // something that the unit test framework recognizes.
  TEST_ASSERT( success );
}

} // namespace (anonymous)
