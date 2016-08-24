/*
HEADER
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

/// \file Ifpack2_UnitTestSGSMT.cpp
/// \brief Unit test for multithreaded symmetric Gauss-Seidel

#include "Teuchos_UnitTestHarness.hpp"
#include "Ifpack2_UnitTestHelpers.hpp"
#include "Ifpack2_Relaxation.hpp"

#include <KokkosKernels_GaussSeidel.hpp>

namespace { // (anonymous)

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
typedef Tpetra::global_size_t GST;
typedef tif_utest::Node Node;
using std::endl;

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(SGS_MT, OneSweep, Scalar, LocalOrdinal, GlobalOrdinal)
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using std::endl;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO, Node> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> mv_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;

  out << "Make sure that one sweep of threaded SGS doesn't crash" << endl;

  RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  if (comm->getSize () == 1) {
    out << "The unit test's (MPI) communicator only contains one process."
        << endl << "This test only makes sense if the communicator contains "
        << "multiple processes." << endl << "I'll let the test pass trivially."
        << endl;
    return;
  }

  const LO lclNumRows = 100;
  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  RCP<const map_type> rowMap (new map_type (INVALID, lclNumRows, 0, comm));
  RCP<const crs_matrix_type> crsMatrix = tif_utest::create_test_matrix<Scalar,LO,GO,Node> (rowMap);

  // We don't really need prec to be a pointer here, but it's
  // convenient for putting the constructor call in a TEST_NOTHROW
  // macro.  Otherwise the declaration of prec would be in an inner
  // scope and we couldn't use it below.
  RCP<Ifpack2::Relaxation<row_matrix_type> > prec;
  {
    TEST_NOTHROW( prec = rcp (new Ifpack2::Relaxation<row_matrix_type> (crsMatrix)) );
  }

  ParameterList params;
  params.set ("relaxation: type", "MT Symmetric Gauss-Seidel");
  prec->setParameters (params);
  TEST_NOTHROW( prec->initialize () );
  TEST_NOTHROW( prec->compute () );

  mv_type X (rowMap, 2);
  X.putScalar (STS::zero ());
  mv_type B (rowMap, 2);
  B.putScalar (STS::one ());

  prec->apply (B, X);
}

// Check that multiple sweeps of Symmetric Gauss-Seidel (SGS) work correctly.
//
// The point is really to ensure that SGS does the Import (from the
// domain Map input X to the column Map version of X) in the right
// place, if an Import is needed.  Thus, this example uses a matrix
// with a nontrivial Import (i.e., where the domain and column Maps
// are not the same).
TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(SGS_MT, MultipleSweeps, Scalar, LO, GO)
{
  using Teuchos::ArrayView;
  using Teuchos::Array;
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using std::endl;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> crs_matrix_type;
  typedef Tpetra::Import<LO, GO, Node> import_type;
  typedef Tpetra::Map<LO, GO, Node> map_type;
  typedef Tpetra::Vector<Scalar, LO, GO, Node> vec_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO, Node> row_matrix_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef Teuchos::ScalarTraits<typename STS::magnitudeType> STM;

  Teuchos::OSTab tab0 (out);
  out << "Test multiple threaded SGS sweeps with nontrivial Import" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  if (numProcs == 1) {
    out << "The unit test's (MPI) communicator only contains one process."
        << endl << "This test only makes sense if the communicator contains "
        << "multiple processes." << endl << "I'll let the test pass trivially."
        << endl;
    return;
  }

  const GST gblNumRows = static_cast<GST> (numProcs * 2);
  const GO indexBase = 0;
  RCP<const map_type> rowMap (new map_type (gblNumRows, indexBase, comm));
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  RCP<crs_matrix_type> A =
    rcp (new crs_matrix_type (rowMap, 6, Tpetra::StaticProfile));

  {
    const size_t lclNumRows = rowMap->getNodeNumElements ();
    const Scalar ONE = STS::one ();
    const Scalar TWO = ONE + ONE;
    const Scalar FOUR = TWO + TWO;
    const Scalar EIGHT = FOUR + FOUR;
    const Scalar TWELVE = EIGHT + FOUR;

    if (myRank == 1) {
      std::cout << ONE << " " << TWO << " " <<  FOUR << " " << EIGHT << " " << TWELVE << std::endl;

    }
    Array<Scalar> vals (6);
    Array<GO> gblColInds (6);

    for (LO lclRow = 0; lclRow < static_cast<LO> (lclNumRows); ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const Scalar diagVal = TWELVE;
      size_t numEnt = 6;
      if (gblRow == rowMap->getMinAllGlobalIndex ()) {
        numEnt = 4;
        vals[0] = diagVal;
        vals[1] = ONE / TWO;
        vals[2] = -ONE;
        vals[3] = ONE;
        gblColInds[0] = rowMap->getMinAllGlobalIndex ();
        gblColInds[1] = gblColInds[0] + 1;
        gblColInds[2] = gblColInds[0] + 2;
        gblColInds[3] = gblColInds[0] + 3;

        //std::cout << "-row: " << lclRow << " vals:" << vals[0] << " " << vals[1] << " " << vals[2] << " " << vals[3] << std::endl;
      }
      else if (gblRow == rowMap->getMinAllGlobalIndex () + 1) {
        numEnt = 4;
        vals[0] = -ONE / TWO;
        vals[1] = diagVal / TWO;
        vals[2] = ONE;
        vals[3] = -ONE;
        gblColInds[0] = rowMap->getMinAllGlobalIndex ();
        gblColInds[1] = gblColInds[0] + 1;
        gblColInds[2] = gblColInds[0] + 2;
        gblColInds[3] = gblColInds[0] + 3;

        //std::cout << "--row: " << lclRow << " vals:" << vals[0] << " " << vals[1] << " " << vals[2] << " " << vals[3] << std::endl;
      }
      else if (gblRow == rowMap->getMaxAllGlobalIndex () - 1) {
        numEnt = 4;
        vals[0] = -ONE;
        vals[1] = ONE;
        vals[2] = diagVal;
        vals[3] = ONE / TWO;
        gblColInds[0] = rowMap->getMaxAllGlobalIndex () - 3;
        gblColInds[1] = rowMap->getMaxAllGlobalIndex () - 2;
        gblColInds[2] = rowMap->getMaxAllGlobalIndex () - 1;
        gblColInds[3] = rowMap->getMaxAllGlobalIndex ();


        //std::cout << "---row: " << lclRow << " gsblRow:" << gblRow << " vals:" << vals[0] << " " << vals[1] << " " << vals[2] << " " << vals[3] << std::endl;
      }
      else if (gblRow == rowMap->getMaxAllGlobalIndex ()) {
        numEnt = 4;
        vals[0] = ONE;
        vals[1] = -ONE;
        vals[2] = -ONE / TWO;
        vals[3] = diagVal / TWO;
        gblColInds[0] = rowMap->getMaxAllGlobalIndex () - 3;
        gblColInds[1] = rowMap->getMaxAllGlobalIndex () - 2;
        gblColInds[2] = rowMap->getMaxAllGlobalIndex () - 1;
        gblColInds[3] = rowMap->getMaxAllGlobalIndex ();
        //std::cout << "----row: " << lclRow << " vals:" << vals[0] << " " << vals[1] << " " << vals[2] << " " << vals[3] << std::endl;
      }
      else if (gblRow % 2 == static_cast<GO> (0)) {
        numEnt = 6;
        vals[0] = -ONE;
        vals[1] = ONE;
        vals[2] = diagVal;
        vals[3] = ONE / TWO;
        vals[4] = -ONE;
        vals[5] = ONE;
        gblColInds[0] = gblRow - 2;
        gblColInds[1] = gblRow - 1;
        gblColInds[2] = gblRow;
        gblColInds[3] = gblRow + 1;
        gblColInds[4] = gblRow + 2;
        gblColInds[5] = gblRow + 3;
        //std::cout << "-----row: " << lclRow << " vals:" << vals[0] << " " << vals[1] << " " << vals[2] << " " << vals[3] << " " << vals[4] << " " << vals[5] << std::endl;
      }
      else { // gblRow % 2 != 0
        numEnt = 6;
        vals[0] = ONE;
        vals[1] = -ONE;
        vals[2] = -ONE / TWO;
        vals[3] = diagVal;
        vals[4] = ONE;
        vals[5] = -ONE;
        gblColInds[0] = gblRow - 3;
        gblColInds[1] = gblRow - 2;
        gblColInds[2] = gblRow - 1;
        gblColInds[3] = gblRow;
        gblColInds[4] = gblRow + 1;
        gblColInds[5] = gblRow + 2;

        //std::cout << "------row: " << lclRow << " vals:" << vals[0] << " " << vals[1] << " " << vals[2] << " " << vals[3] << " " << vals[4] << " " << vals[5] << std::endl;
      }

      ArrayView<Scalar> valsView = vals (0, numEnt);
      ArrayView<GO> gblColIndsView = gblColInds (0, numEnt);
      A->insertGlobalValues (gblRow, gblColIndsView, valsView);
    }
  }
  A->fillComplete (domainMap, rangeMap);

  RCP<const map_type> gatherRowMap;
  {
    const size_t lclNumRows = (myRank == 0) ?
      static_cast<size_t> (gblNumRows) :
      static_cast<size_t> (0);
    gatherRowMap = rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
  }

  RCP<const map_type> gatherDomainMap = gatherRowMap;
  RCP<const map_type> gatherRangeMap = gatherRowMap;

  RCP<crs_matrix_type> A_gather = rcp (new crs_matrix_type (gatherRowMap, 6));
  import_type import (rowMap, gatherRowMap);
  A_gather->doImport (*A, import, Tpetra::INSERT);
  A_gather->fillComplete (gatherDomainMap, gatherRangeMap);

  vec_type X (domainMap);
  vec_type Y (rangeMap);
  vec_type X_gather (gatherDomainMap);
  vec_type Y_gather (gatherRangeMap);
  vec_type Y_diff (gatherRangeMap);

  // Test Symmetric Gauss-Seidel (SGS) with three sweeps.
  // Start by letting SGS set the starting solution to zero.
  ParameterList params;

  params.set ("relaxation: type", "MT Symmetric Gauss-Seidel");
  params.set ("relaxation: sweeps", 3);
  params.set ("relaxation: zero starting solution", true);

  Ifpack2::Relaxation<row_matrix_type> prec (A);
  prec.setParameters (params);
  TEST_NOTHROW( prec.initialize () );
  TEST_NOTHROW( prec.compute () );

  Ifpack2::Relaxation<row_matrix_type> gatherPrec (A_gather);
  gatherPrec.setParameters (params);
  TEST_NOTHROW( gatherPrec.initialize () );
  TEST_NOTHROW( gatherPrec.compute () );

  X.randomize ();
  X_gather.doImport (X, import, Tpetra::REPLACE);
  Y.randomize ();
  X_gather.doImport (X, import, Tpetra::REPLACE);

  prec.apply (X, Y);

  gatherPrec.apply (X_gather, Y_gather);
  Y_diff.doImport (Y, import, Tpetra::REPLACE);
  Y_diff.update (STS::one (), Y_gather, -STS::one ());

  typename STS::magnitudeType normInf = Y_diff.normInf ();
  TEST_EQUALITY(normInf, STM::zero ());

  out << "Repeat test without setting starting solution to zero" << endl;

  params.set ("relaxation: zero starting solution", false);

  prec.setParameters (params);
  TEST_NOTHROW( prec.initialize () );
  TEST_NOTHROW( prec.compute () );

  gatherPrec.setParameters (params);
  TEST_NOTHROW( gatherPrec.initialize () );
  TEST_NOTHROW( gatherPrec.compute () );

  X.randomize ();
  X_gather.doImport (X, import, Tpetra::REPLACE);
  Y.randomize ();
  X_gather.doImport (X, import, Tpetra::REPLACE);

  prec.apply (X, Y);
  gatherPrec.apply (X_gather, Y_gather);

  Y_diff.doImport (Y, import, Tpetra::REPLACE);
  Y_diff.update (STS::one (), Y_gather, -STS::one ());
  normInf = Y_diff.normInf ();
  TEST_EQUALITY( normInf, STM::zero () );
}


// Macro used inside the unit test below.  It tests for global error,
// and if so, prints each process' error message and quits the test
// early.
//
// 'out' only prints on Process 0.  It's really not OK for other
// processes to print to stdout, but it usually works and we need to
// do it for debugging.
#define IFPACK2RELAXATION_REPORT_GLOBAL_ERR( WHAT_STRING ) do { \
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess)); \
  TEST_EQUALITY_CONST( gblSuccess, 1 ); \
  if (gblSuccess != 1) { \
    out << WHAT_STRING << " FAILED on one or more processes!" << endl; \
    for (int p = 0; p < numProcs; ++p) { \
      if (myRank == p && lclSuccess != 1) { \
        std::cout << errStrm.str () << std::flush; \
      } \
      comm->barrier (); \
      comm->barrier (); \
      comm->barrier (); \
    } \
    std::cerr << "TEST FAILED; RETURNING EARLY" << endl; \
    return; \
  } \
} while (false)

#define UNIT_TEST_GROUP_SC_LO_GO( Scalar, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( SGS_MT, OneSweep, Scalar, LO, GO ) \
    TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( SGS_MT, MultipleSweeps, Scalar, LO, GO )

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// TODO (mfh 24 Aug 2016) Test complex Scalar types

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)


