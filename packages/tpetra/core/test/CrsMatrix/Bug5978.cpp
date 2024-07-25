// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This test is in response to bug documented in:
// https://software.sandia.gov/bugzilla/show_bug.cgi?id=5978
//
// Since the bugzilla server is not generally visible, a relevant comment from
// mhoemme describing the comment is the following:
//
// Ok, I think I get it.
//
// On both processes, the graph starts neither locally nor globally indexed. In
// your test, only Process 0 calls insertLocalValues(), since Process 1 has zero
// rows. The very first call to insertLocalValues() on Process 0 calls
// allocateValues(LocalIndices, GraphNotYetAllocated). This in turn calls
// allocateIndices(LocalIndices). This marks the graph as locally indexed, and
// as _not_ globally indexed:
//
// indicesAreLocal_  = (lg == LocalIndices);
// indicesAreGlobal_ = (lg == GlobalIndices);
//
// When your test calls fillComplete() on the matrix, indices have only been
// allocated (i.e., allocateIndices() has been called) only on Process 0, not on
// Process 1. Thus, in fillComplete, only Process 1 calls
// allocateValues(GlobalIndices, GraphNotYetAllocated). This sets the locally
// vs. globally indexed states to the opposite of what they are on Process 0.
//
// A comment above this line in fillComplete says "Allocate global, in case we
// do not have a column Map yet." The problem is that we _do_ have a column Map:
// your test gives it to the matrix's constructor.
//
// mhoemme indicates that the following commits provide fixes (but the commits
// are not in the Trilinos repository):
//
// - a0120d54b64fa0bb4718fd5aa4430c6ad9612512
// - 210b4176780172382181d05e815be92117d4324f
//

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <TpetraCore_ETIHelperMacros.h>
#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_TestingUtilities.hpp>

// @mhoemme: When compiled, this test gives a lot of warnings coming from the
// @mhoemme: Teuchos unit testing stuff.  I'm assuming because I don't return or call
// @mhoemme: any kind of assert?  Does that make sense?  I'd like to remove all the
// @mhoemme: warnings.

namespace { // anonymous

using std::endl;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::OrdinalTraits;
using Tpetra::Map;
using Tpetra::CrsMatrix;

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug5978, SC, LO, GO, NT)
{
  typedef Map<LO,GO,NT> MapType;
  typedef CrsMatrix<SC, LO, GO, NT> MatrixType;

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  const int numProc = comm->getSize ();
  const int myRank = comm->getRank ();

  if (numProc != 2) {
    out << "This test must be run with exactly 2 MPI processes, but you ran "
        << "it with " << numProc << " process" << (numProc != 1 ? "es" : "")
        << "." << endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "This test must be run with exactly "
      "2 MPI processes, but you ran it with " << numProc << " process"
      << (numProc != 1 ? "es" : "") << ".");
  }
  out << "Proc " << myRank << ": Running test" << endl;
  comm->barrier ();

  Array<GO> rowGIDs, colGIDs;
  if (myRank == 0) {
    rowGIDs.push_back (0);
    rowGIDs.push_back (1);
    colGIDs.push_back (0);
    colGIDs.push_back (1);
  }

  out << "Proc " << myRank << ": Creating row and column Maps" << endl;
  comm->barrier ();

  const GO INVALID = OrdinalTraits<GO>::invalid ();
  RCP<const MapType> rowMap =
    rcp (new MapType (INVALID, rowGIDs (), 0, comm));
  RCP<const MapType> colMap =
    rcp (new MapType (INVALID, colGIDs (), 0, comm));

  out << "Proc " << myRank << ": Creating matrix" << endl;
  comm->barrier ();

  ArrayRCP<size_t> count (rowGIDs.size ());

  // mfh 16 Jan 2014: The explicit cast will hopefully prevent gcc
  // 4.3.3 from confusing the non-template method assign(size_type,
  // const T&) from the template method assign(Iter, Iter).
  count.assign (rowGIDs.size (), static_cast<size_t> (1));
  MatrixType A (rowMap, colMap, count ());

  out << "Proc " << myRank << ": Filling matrix" << endl;
  comm->barrier ();

  Array<LO> column (1);
  Array<SC> entry (1, 7);
  for (int i = 0; i < int (rowGIDs.size ()); ++i) {
    column[0] = i;
    A.insertLocalValues (i, column (), entry ());
  }

  out << "Proc " << myRank << ": Calling fillComplete" << endl;
  comm->barrier ();

  A.fillComplete (colMap, rowMap);

  comm->barrier ();
  out << "Proc " << myRank << ": Done with test" << endl;
}

#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug5978, SC, LO, GO, NT)

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(UNIT_TEST_GROUP_SC_LO_GO_NO)

} // namespace (anonymous)
