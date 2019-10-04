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

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
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
  int lclSuccess = 1;
  int gblSuccess = 0;

  out << "CrsMatrix: Test sumIntoLocalValues, StaticProfile, "
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
  crs_matrix_type A (rowAndColMap, rowAndColMap, maxNumEntPerRow,
                     Tpetra::StaticProfile);

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

    const double vals[2] = { 5.0, 6.0 };
    const LO inds[2] = { LO (1), LO (2) };
    const LO numToInsert = 2;

    for (LO lclRow = 0; lclRow < lclNumInds; ++lclRow) {
      A.insertLocalValues (lclRow, numToInsert, vals, inds);
      const size_t newNumEnt = A.getNumEntriesInLocalRow (lclRow);
      TEST_ASSERT( newNumEnt == 2 );

      Teuchos::ArrayView<const LO> inds_av;
      Teuchos::ArrayView<const double> vals_av;
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

    const double vals[2] = { 20.0 };
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
