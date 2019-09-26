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
// ************************************************************************
// @HEADER
*/

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include <algorithm>
#include <array>

namespace { // (anonymous)

using Teuchos::Comm;
using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using crs_matrix_type = Tpetra::CrsMatrix<>;
using map_type = Tpetra::Map<>;
using ST = crs_matrix_type::scalar_type;
using LO = map_type::local_ordinal_type;
using GO = map_type::global_ordinal_type;

RCP<const map_type>
makeRowMap (RCP<const Comm<int>> comm, const LO lclNumRows)
{
  const GO indexBase = 0;
  const GO gblNumRows = GO (comm->getSize ()) * GO (lclNumRows);
  return rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
}

RCP<const map_type>
makeColumnMap (RCP<const Comm<int>> comm)
{
  const GO indexBase = 0;
  const std::array<GO, 1> myGblInds {{indexBase}};
  const GO gblNumCols (comm->getSize ());
  return rcp (new map_type (gblNumCols, myGblInds.data (),
                            LO (myGblInds.size ()), indexBase, comm));
}

RCP<const map_type>
makeDomainMap (RCP<const Comm<int>> comm)
{
  const LO lclNumInds = (comm->getRank () == 0) ? 1 : 0;
  const GO gblNumInds = 1;
  const GO indexBase = 0;
  return rcp (new map_type (gblNumInds, lclNumInds, indexBase, comm));
}

RCP<const map_type>
makeRangeMap (RCP<const Comm<int>> comm, const LO lclNumRows)
{
  const GO indexBase = 0;
  const GO gblNumRows = GO (comm->getSize ()) * GO (lclNumRows);
  return rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
}

RCP<crs_matrix_type>
createTestMatrix (RCP<const Comm<int>> comm, const LO lclNumRows)
{
  RCP<crs_matrix_type> A;
#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  A = rcp (new crs_matrix_type (makeRowMap (comm, lclNumRows),
                                makeColumnMap (comm), size_t (1),
                                Tpetra::StaticProfile));
#else
  A = rcp (new crs_matrix_type (makeRowMap (comm, lclNumRows),
                                makeColumnMap (comm), size_t (1)));
#endif // TPETRA_ENABLE_DEPRECATED_CODE

  constexpr LO numEnt = 1;
  std::array<LO, 1> lclCols {{ 0 }};
  std::array<ST, 1> vals {{ 1.0 }};
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    A->insertLocalValues (lclRow, numEnt, vals.data (), lclCols.data ());
  }
  A->fillComplete (makeDomainMap (comm),
                   makeRangeMap (comm, lclNumRows));
  return A;
}

void
standardTransposeTests (bool& success,
                        Teuchos::FancyOStream& out,
                        const crs_matrix_type& A,
                        const crs_matrix_type& AT,
                        const Comm<int>& comm)
{
  int lclSuccess = 1;
  int gblSuccess = 0;

  auto G = A.getCrsGraph ();
  TEST_ASSERT( ! G.is_null () );
  auto GT = AT.getCrsGraph ();
  TEST_ASSERT( ! GT.is_null () );

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0;
  reduceAll (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_ASSERT( gblSuccess == 1 );
  if (! success) {
    return;
  }

  auto colMap = A.getColMap ();
  TEST_ASSERT( ! colMap.is_null () );
  auto domMap = A.getDomainMap ();
  TEST_ASSERT( ! domMap.is_null () );
  auto ranMap = A.getRangeMap ();
  TEST_ASSERT( ! ranMap.is_null () );
  auto rowMap = A.getRowMap ();
  TEST_ASSERT( ! rowMap.is_null () );

  auto colMap_at = AT.getColMap ();
  TEST_ASSERT( ! colMap_at.is_null () );
  auto domMap_at = AT.getDomainMap ();
  TEST_ASSERT( ! domMap_at.is_null () );
  auto ranMap_at = AT.getRangeMap ();
  TEST_ASSERT( ! ranMap_at.is_null () );
  auto rowMap_at = AT.getRowMap ();
  TEST_ASSERT( ! rowMap_at.is_null () );

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0;
  reduceAll (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_ASSERT( gblSuccess == 1 );
  if (! success) {
    return;
  }

  // createTranspose redistributes the result to have a nonshared row
  // Map (row Map == range Map), so the column Map of the transpose
  // may not be the same as the row Map of the original, and the row
  // Map of the transpose may not be the same as the column Map of the
  // original.

  TEST_ASSERT( domMap_at->isSameAs (*ranMap) );
  TEST_ASSERT( ranMap_at->isSameAs (*domMap) );

  auto G_imp = G->getImporter ();
  TEST_ASSERT( colMap->isSameAs (*domMap) || ! G_imp.is_null () );
  auto GT_imp = GT->getImporter ();
  TEST_ASSERT( colMap_at->isSameAs (*domMap_at) || ! GT_imp.is_null () );

  auto G_exp = G->getExporter ();
  TEST_ASSERT( rowMap->isSameAs (*ranMap) || ! G_exp.is_null () );
  auto GT_exp = GT->getExporter ();
  TEST_ASSERT( rowMap_at->isSameAs (*ranMap_at) || ! GT_exp.is_null () );

  lclSuccess = success ? 1 : 0;
  gblSuccess = 0;
  reduceAll (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_ASSERT( gblSuccess == 1 );
  if (! success) {
    return;
  }
}

void
testTranspose (bool& success,
               Teuchos::FancyOStream& out,
               RCP<const Comm<int>> comm,
               const LO lclNumRows)
{
  int lclSuccess = 1;
  int gblSuccess = 0;

  RCP<crs_matrix_type> A = createTestMatrix (comm, lclNumRows);
  {
    out << "Test with default \"sort\"" << std::endl;
    Teuchos::OSTab tab1 (out);

    RCP<crs_matrix_type> AT_unsorted = [&] () {
      Tpetra::RowMatrixTransposer<> transposer (A);
      return transposer.createTranspose ();
    } ();

    TEST_ASSERT( ! AT_unsorted.is_null () );
    standardTransposeTests (success, out, *A, *AT_unsorted, *comm);

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (! success) {
      return;
    }

    // By this point, we know that the row Maps are the same.
    auto rowMap_at = AT_unsorted->getRowMap ();
    const LO lclNumRows_at (rowMap_at->getNodeNumElements ());
    if (lclNumRows_at != 0) {
      std::vector<LO> lclColIndsBuf;
      std::vector<ST> valsBuf;
      for (LO lclRow = 0; lclRow < lclNumRows_at; ++lclRow) {
        Teuchos::ArrayView<const LO> lclColInds;
        Teuchos::ArrayView<const ST> vals;
        AT_unsorted->getLocalRowView (lclRow, lclColInds, vals);

        const GO gblNumRows = GO (lclNumRows) * GO (comm->getSize ());
        TEST_ASSERT( LO (lclColInds.size ()) == LO (gblNumRows) );
        TEST_ASSERT( LO (vals.size ()) == LO (gblNumRows) );

        if (success) {
          // Result's rows may not be sorted.
          lclColIndsBuf.resize (lclColInds.size ());
          valsBuf.resize (vals.size ());
          std::copy (lclColInds.begin (), lclColInds.end (),
                     lclColIndsBuf.begin ());
          std::copy (vals.begin (), vals.end (), valsBuf.begin ());
          Tpetra::sort2 (lclColIndsBuf.begin (), lclColIndsBuf.end (),
                         valsBuf.begin ());

          bool good = true;
          for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
            if (lclColIndsBuf[lclRow] != lclRow) {
              good = false;
              break;
            }
            else if (valsBuf[lclRow] != 1.0) {
              good = false;
              break;
            }
          }
          TEST_ASSERT( good );
        }
      }
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (! success) {
      return;
    }
  }

  {
    out << "Test with \"sort\" = true" << std::endl;
    Teuchos::OSTab tab1 (out);

    RCP<crs_matrix_type> AT_sorted = [&] () {
      auto params = Teuchos::parameterList ("Tpetra::RowMatrixTransposer");
      params->set ("sort", true);
      Tpetra::RowMatrixTransposer<> transposer (A);
      return transposer.createTranspose (params);
    } ();

    TEST_ASSERT( ! AT_sorted.is_null () );
    standardTransposeTests (success, out, *A, *AT_sorted, *comm);

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (! success) {
      return;
    }

    // By this point, we know that the row Maps are the same.
    auto rowMap_at = AT_sorted->getRowMap ();
    const LO lclNumRows_at (rowMap_at->getNodeNumElements ());
    if (lclNumRows_at != 0) {
      for (LO lclRow = 0; lclRow < lclNumRows_at; ++lclRow) {
        Teuchos::ArrayView<const LO> lclColInds;
        Teuchos::ArrayView<const ST> vals;
        AT_sorted->getLocalRowView (lclRow, lclColInds, vals);

        const GO gblNumRows = GO (lclNumRows) * GO (comm->getSize ());
        TEST_ASSERT( LO (lclColInds.size ()) == LO (gblNumRows) );
        TEST_ASSERT( LO (vals.size ()) == LO (gblNumRows) );

        if (success) {
          bool good = true;
          for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
            if (lclColInds[lclRow] != lclRow) {
              good = false;
              break;
            }
            else if (vals[lclRow] != 1.0) {
              good = false;
              break;
            }
          }
          TEST_ASSERT( good );
        }
      }
    }

    lclSuccess = success ? 1 : 0;
    gblSuccess = 0;
    reduceAll (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_ASSERT( gblSuccess == 1 );
    if (! success) {
      return;
    }
  }
}

TEUCHOS_UNIT_TEST( CrsMatrixTranspose, SortRows )
{
  const LO lclNumRows = 1000;
  RCP<const Comm<int>> comm = Tpetra::getDefaultComm ();
  testTranspose (success, out, comm, lclNumRows );
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  int errCode = 0;
  {
    Tpetra::ScopeGuard tpetraScope (&argc, &argv);
    errCode = Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  }
  return errCode;
}
