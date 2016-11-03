/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

/// \file 762.cpp
/// \brief Regression test for Github Issue #762

#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Ifpack2_UnitTestHelpers.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include <sstream>

namespace { // (anonymous)

using Tpetra::global_size_t;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::cerr;
using std::endl;

typedef tif_utest::Node Node;

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Ifpack2, Issue762, Scalar, LO, GO)
{
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO,Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  // Whether to print to cerr, for more immediate output ('out'
  // doesn't actually print until the test returns, so it's not so
  // helpful if the test segfaults).
  const bool printToCerr = true;

  int lclSuccess = 1; // to be updated below
  int gblSuccess = 0; // output argument
  std::ostringstream errStrm;

  if (printToCerr) {
    cerr << "Ifpack2: Test #762" << endl;
  }
  else {
    out << "Ifpack2: Test #762" << endl;
  }
  // This is a scope guard, so it's better to do this unconditionally.
  Teuchos::OSTab tab1 (out);

  const global_size_t num_rows_per_proc = 5;
  RCP<const map_type> rowmap =
    tif_utest::create_tpetra_map<LO, GO, Node> (num_rows_per_proc);
  TEST_ASSERT( ! rowmap.is_null () );
  if (rowmap.is_null ()) {
    return; // that's the best we can do
  }
  RCP<const Teuchos::Comm<int> > comm = rowmap->getComm ();
  TEST_ASSERT( ! comm.is_null () );
  if (comm.is_null ()) {
    return; // that's the best we can do
  }
  if (comm->getSize () > 1) {
    if (printToCerr) {
      cerr << "This test may only be run in serial "
        "or with a single MPI process." << endl;
    }
    else {
      out << "This test may only be run in serial "
        "or with a single MPI process." << endl;
    }
    return;
  }

  //out << "Creating matrix" << endl;
  cerr << "Creating matrix" << endl;
  RCP<const crs_matrix_type> crsmatrix;
  try {
    crsmatrix = tif_utest::create_test_matrix2<Scalar,LO,GO,Node>(rowmap);
    if (printToCerr) {
      cerr << "create_test_matrix2 returned!" << endl;
    }
    // Don't print to 'out' here, since we only care about per-process
    // printing in this case.
  }
  catch (std::exception& e) {
    success = false;
    errStrm << "Process " << comm->getRank () << ": create_test_matrix2 threw "
      "an exception: " << e.what () << endl;
  }
  catch (...) {
    success = false;
    errStrm << "Process " << comm->getRank () << ": create_test_matrix2 threw "
      "an exception not a subclass of std::exception." << endl;
  }
  TEST_ASSERT( ! crsmatrix.is_null () );
  lclSuccess = (lclSuccess == 0 || ! success) ? 0 : 1;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    Tpetra::Details::gathervPrint (out, errStrm.str (), *comm);
    return;
  }

  if (printToCerr) {
    cerr << "Creating test problem" << endl;
  }
  else {
    out << "Creating test problem" << endl;
  }
  MV x (rowmap, 2);
  MV y (rowmap, 2);
  x.putScalar (STS::one ());

  {
    // First test alpha == 1 and beta == 0.
    const Scalar alpha = Teuchos::as<Scalar> (1.0);
    const Scalar beta = Teuchos::as<Scalar> (0.0);

    if (printToCerr) {
      cerr << "Creating copies of x and y" << endl;
    }
    else {
      out << "Creating copies of x and y" << endl;
    }
    MV x_copy = Tpetra::createCopy (x);
    MV y_copy = Tpetra::createCopy (y);

    if (printToCerr) {
      cerr << "Testing apply() for alpha = " << alpha
           << " and beta = " << beta << endl;
    }
    else {
      out << "Testing apply() for alpha = " << alpha
          << " and beta = " << beta << endl;
    }
    try {
      crsmatrix->apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
      if (printToCerr) {
        cerr << "apply (alpha = " << alpha << ", beta = " << beta << ") returned!"
             << endl;
      }
      // Don't print to 'out' here, since we only care about per-process
      // printing in this case.
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << comm->getRank () << ": CrsMatrix::apply threw "
        "an exception: " << e.what () << endl;
    }
    catch (...) {
      lclSuccess = 0;
      errStrm << "Process " << comm->getRank () << ": CrsMatrix::apply threw "
        "an exception not a subclass of std::exception." << endl;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      Tpetra::Details::gathervPrint (out, errStrm.str (), *comm);
      return;
    }
  }

  x.putScalar (STS::one ());
  y.putScalar (STS::zero ());
  {
    // Now test alpha != 1 and beta == 0.
    const Scalar alpha = Teuchos::as<Scalar> (2.0);
    const Scalar beta = Teuchos::as<Scalar> (0.0);

    if (printToCerr) {
      cerr << "Creating copies of x and y" << endl;
    }
    else {
      out << "Creating copies of x and y" << endl;
    }
    MV x_copy = Tpetra::createCopy (x);
    MV y_copy = Tpetra::createCopy (y);

    if (printToCerr) {
      cerr << "Testing apply() for alpha = " << alpha
           << " and beta = " << beta << endl;
    }
    else {
      out << "Testing apply() for alpha = " << alpha
          << " and beta = " << beta << endl;
    }
    try {
      // This (alpha = 2, beta = 0, 2 columns) was the #762 case.
      crsmatrix->apply (x_copy, y_copy, Teuchos::NO_TRANS, alpha, beta);
      cerr << "apply (alpha = " << alpha << ", beta = " << beta << ") returned!"
           << endl;
      // Don't print to 'out' here, since we only care about per-process
      // printing in this case.
    }
    catch (std::exception& e) {
      lclSuccess = 0;
      errStrm << "Process " << comm->getRank () << ": CrsMatrix::apply threw "
        "an exception: " << e.what () << endl;
    }
    catch (...) {
      lclSuccess = 0;
      errStrm << "Process " << comm->getRank () << ": CrsMatrix::apply threw "
        "an exception not a subclass of std::exception." << endl;
    }
    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY_CONST( gblSuccess, 1 );
    if (gblSuccess != 1) {
      Tpetra::Details::gathervPrint (out, errStrm.str (), *comm);
      return;
    }
    else {
      if (comm->getRank () == 0) {
        if (printToCerr) {
          cerr << "Yay, got through the test!" << endl;
        }
        else {
          out << "Yay, got through the test!" << endl;
        }
      }
    }
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO( SC, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Ifpack2, Issue762, SC, LO, GO )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types.

IFPACK2_INSTANTIATE_SLG( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

