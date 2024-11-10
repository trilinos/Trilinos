// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Teuchos_CommHelpers.hpp"

namespace { // (anonymous)

using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;

TEUCHOS_UNIT_TEST( CrsMatrix, sumIntoStaticProfileExtraSpace )
{
  using map_type = Tpetra::Map<>;
  using crs_matrix_type = Tpetra::CrsMatrix<>;
  using LO = map_type::local_ordinal_type;
  using GO = map_type::global_ordinal_type;
  using SC = Tpetra::Vector<>::scalar_type;
  int lclSuccess = 1;
  int gblSuccess = 0;

  out << "CrsMatrix: Test sumIntoLocalValues, "
    "not yet fill complete, extra space in each row" << endl;
  Teuchos::OSTab tab1 (out);

  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm ();

  // We need at least one nonzero local column index -- I'm presuming
  // that local column indices default-fill with zeros.

  const int numProcs = comm->getSize ();
  const LO lclNumInds = 5;
  const GO gblNumInds = GO (numProcs) * GO (lclNumInds);
  const GO indexBase = 0;

  RCP<const map_type> rowAndColMap
    (new map_type (gblNumInds, lclNumInds, indexBase, comm));

  constexpr size_t maxNumEntPerRow = 5;
  crs_matrix_type A (rowAndColMap, rowAndColMap, maxNumEntPerRow);

  for (LO lclRow = 0; lclRow < lclNumInds; ++lclRow) {
    const size_t numEnt = A.getNumEntriesInLocalRow (lclRow);
    TEST_ASSERT( numEnt == 0 );
  }

  lclSuccess = success ? 1 : 0;
  reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "Exiting early" << endl;
    return;
  }

  {
    out << "Insert local column indices 1,2 into each row" << endl;
    Teuchos::OSTab tab2 (out);

    const SC vals[2] = { 5.0, 6.0 };
    const LO inds[2] = { LO (1), LO (2) };
    const LO numToInsert = 2;

    for (LO lclRow = 0; lclRow < lclNumInds; ++lclRow) {
      A.insertLocalValues (lclRow, numToInsert, vals, inds);
      const size_t newNumEnt = A.getNumEntriesInLocalRow (lclRow);
      TEST_ASSERT( newNumEnt == 2 );

      typename crs_matrix_type::local_inds_host_view_type inds_av;
      typename crs_matrix_type::values_host_view_type vals_av;;
      A.getLocalRowView (lclRow, inds_av, vals_av);
      TEST_ASSERT( inds_av.size () == ptrdiff_t (2) );
      TEST_ASSERT( vals_av.size () == ptrdiff_t (2) );
      if (success) {
        TEST_ASSERT( inds_av[0] == inds[0] && inds_av[1] == inds[1] );
        TEST_ASSERT( vals_av[0] == vals[0] && vals_av[1] == vals[1] );
      }
    }
  }

  lclSuccess = success ? 1 : 0;
  reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "Exiting early" << endl;
    return;
  }

  {
    out << "Try to sumInto a local index that should not exist" << endl;
    Teuchos::OSTab tab2 (out);

    const SC vals[2] = { 20.0 };
    const LO inds[2] = { LO (0) }; // not in graph/matrix yet
    const LO numToInsert = 1;

    for (LO lclRow = 0; lclRow < lclNumInds; ++lclRow) {
      out << "lclRow=" << lclRow << ": " << endl;
      Teuchos::OSTab tab3 (out);

      const size_t oldNumEnt = A.getNumEntriesInLocalRow (lclRow);
      out << "Before: numEnt=" << oldNumEnt << endl;

      const LO numModified = A.sumIntoLocalValues (lclRow, numToInsert, vals, inds);
      out << "numModified=" << numModified << endl;
      TEST_ASSERT( numModified == 0 );

      const size_t newNumEnt = A.getNumEntriesInLocalRow (lclRow);
      out << "After: numEnt=" << newNumEnt << endl;
      TEST_ASSERT( newNumEnt == 2 );
    }
  }

  lclSuccess = success ? 1 : 0;
  reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY_CONST( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "Exiting early" << endl;
    return;
  }
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
