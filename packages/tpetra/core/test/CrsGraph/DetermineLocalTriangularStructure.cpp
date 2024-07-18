// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file DetermineLocalTriangularStructure.cpp
/// \brief Unit test for Tpetra::Details::determineLocalTriangularStructure

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_determineLocalTriangularStructure.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include <type_traits>

namespace { // (anonymous)

template<class LocalOrdinalType, class DeviceType>
using LocalCrsGraph = Kokkos::StaticCrsGraph<LocalOrdinalType, Kokkos::LayoutLeft, DeviceType>;

// template<class ScalarType, class LocalOrdinalType, class DeviceType>
// using LocalCrsMatrix =
//   KokkosSparse::CrsMatrix<ScalarType, LocalOrdinalType, DeviceType, void,
//                           typename KokkosCrsGraph<LocalOrdinalType, DeviceType>::size_type>;

template<class LocalOrdinalType, class GlobalOrdinalType, class DeviceType>
using LocalMap = typename Tpetra::Map<LocalOrdinalType, GlobalOrdinalType, DeviceType>::local_map_type;

template<class LO, class GO, class NT>
LocalCrsGraph<LO, typename NT::device_type>
makeDiagonalGraph (const LO lclNumRows,
                   const Tpetra::Map<LO, GO, NT>& rowMap,
                   const Tpetra::Map<LO, GO, NT>& colMap)
{
  using crs_graph_type = LocalCrsGraph<LO, typename NT::device_type>;
  using row_map_type = typename crs_graph_type::row_map_type::non_const_type;
  using entries_type = typename crs_graph_type::entries_type::non_const_type;

  row_map_type ptr ("ptr", lclNumRows+1);
  entries_type ind ("ind", lclNumRows);

  auto ptr_h = Kokkos::create_mirror_view (ptr);
  auto ind_h = Kokkos::create_mirror_view (ind);

  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    ptr_h[lclRow] = lclRow;
    const GO gblRow = rowMap.getGlobalElement (lclRow);
    const GO gblCol = gblRow;
    const LO lclCol = colMap.getLocalElement (gblCol);
    ind_h[lclRow] = lclCol;
  }
  ptr_h[lclNumRows] = lclNumRows;

  Kokkos::deep_copy (ptr, ptr_h);
  Kokkos::deep_copy (ind, ind_h);
  return crs_graph_type (ind, ptr);
}


template<class LO, class GO, class NT>
LocalCrsGraph<LO, typename NT::device_type>
makeTriangularGraph (const LO lclNumRows,
                     const Tpetra::Map<LO, GO, NT>& rowMap,
                     const Tpetra::Map<LO, GO, NT>& colMap,
                     const bool lowerTriangular, // false means upper
                     const bool explicitDiagonal)
{
  using crs_graph_type = LocalCrsGraph<LO, typename NT::device_type>;
  using row_map_type = typename crs_graph_type::row_map_type::non_const_type;
  using entries_type = typename crs_graph_type::entries_type::non_const_type;
  using size_type = typename crs_graph_type::size_type;

  // two entries per row, except one row (first row for lower
  // triangular; last row for upper triangular).
  constexpr LO ONE {1};
  const size_type lclNumEnt = [&] () {
    if (lclNumRows == 0) {
      return size_type (0);
    }
    if (lclNumRows == ONE) {
      return explicitDiagonal ? size_type (1) : size_type (0);
    }
    else {
      if (explicitDiagonal) {
        return size_type (lclNumRows + (lclNumRows - ONE));
      }
      else {
        return size_type (lclNumRows - ONE);
      }
    }
  } ();

  row_map_type ptr ("ptr", lclNumRows+1);
  entries_type ind ("ind", lclNumEnt);

  auto ptr_h = Kokkos::create_mirror_view (ptr);
  auto ind_h = Kokkos::create_mirror_view (ind);

  if (lclNumRows == 0) {
    ptr_h[0] = 0;
    Kokkos::deep_copy (ptr, ptr_h); // ind is empty
    return crs_graph_type (ind, ptr);
  }
  if (lclNumRows == 1) {
    ptr_h[0] = 0;
    ptr_h[1] = 1;
    ind_h[0] = colMap.getLocalElement (rowMap.getGlobalElement (0));
    Kokkos::deep_copy (ptr, ptr_h);
    Kokkos::deep_copy (ind, ind_h);
    return crs_graph_type (ind, ptr);
  }

  const LO lclStartRow = 1;
  const LO lclEndRow = lclNumRows - LO (1);
  size_type curPos = 0;

  {
    const LO lclRow = 0;
    ptr_h[lclRow] = curPos;
    const GO gblDiagCol = rowMap.getGlobalElement (lclRow);
    if (explicitDiagonal) {
      ind_h[curPos++] = colMap.getLocalElement (gblDiagCol);
    }
    if (! lowerTriangular) {
      ind_h[curPos++] = colMap.getLocalElement (gblDiagCol + 1);
    }
  }
  for (LO lclRow = lclStartRow; lclRow < lclEndRow; ++lclRow) {
    ptr_h[lclRow] = curPos;
    const GO gblDiagCol = rowMap.getGlobalElement (lclRow);
    if (lowerTriangular) {
      ind_h[curPos++] = colMap.getLocalElement (gblDiagCol - 1);
    }
    if (explicitDiagonal) {
      ind_h[curPos++] = colMap.getLocalElement (gblDiagCol);
    }
    if (! lowerTriangular) {
      ind_h[curPos++] = colMap.getLocalElement (gblDiagCol + 1);
    }
  }
  {
    const LO lclRow = lclNumRows - 1;
    ptr_h[lclRow] = curPos;
    const GO gblDiagCol = rowMap.getGlobalElement (lclRow);
    if (lowerTriangular) {
      ind_h[curPos++] = colMap.getLocalElement (gblDiagCol - 1);
    }
    if (explicitDiagonal) {
      ind_h[curPos++] = colMap.getLocalElement (gblDiagCol);
    }
  }
  ptr_h[lclNumRows] = curPos;

  TEUCHOS_TEST_FOR_EXCEPTION
    (curPos != lclNumEnt, std::logic_error, "Bug in test: "
     "curPos = " << curPos << " != lclNumEnt = " << lclNumEnt << ".");

  Kokkos::deep_copy (ptr, ptr_h);
  Kokkos::deep_copy (ind, ind_h);
  return crs_graph_type (ind, ptr);
}


template<class LO, class GO, class NT>
LocalCrsGraph<LO, typename NT::device_type>
makeTridiagonalGraph (const LO lclNumRows,
                      const Tpetra::Map<LO, GO, NT>& rowMap,
                      const Tpetra::Map<LO, GO, NT>& colMap)
{
  using crs_graph_type = LocalCrsGraph<LO, typename NT::device_type>;
  using row_map_type = typename crs_graph_type::row_map_type::non_const_type;
  using entries_type = typename crs_graph_type::entries_type::non_const_type;
  using size_type = typename crs_graph_type::size_type;

  // two entries per row, except one row (first row for lower
  // triangular; last row for upper triangular).
  const size_type lclNumRows_st = static_cast<size_type> (lclNumRows);
  constexpr LO ONE {1};
  constexpr LO TWO {2};
  const size_type lclNumEnt = (lclNumRows < TWO) ?
    size_type (1) :
    size_type (4) + size_type (3) * (lclNumRows_st - size_type (2));

  row_map_type ptr ("ptr", lclNumRows+1);
  entries_type ind ("ind", lclNumEnt);

  auto ptr_h = Kokkos::create_mirror_view (ptr);
  auto ind_h = Kokkos::create_mirror_view (ind);

  if (lclNumRows == 0) {
    ptr_h[0] = 0;
    Kokkos::deep_copy (ptr, ptr_h); // ind is empty
    return crs_graph_type (ind, ptr);
  }
  if (lclNumRows == 1) {
    ptr_h[0] = 0;
    ptr_h[1] = 1;
    ind_h[0] = colMap.getLocalElement (rowMap.getGlobalElement (0));
    Kokkos::deep_copy (ptr, ptr_h);
    Kokkos::deep_copy (ind, ind_h);
    return crs_graph_type (ind, ptr);
  }

  const LO lclStartRow {1};
  const LO lclEndRow {lclNumRows - ONE}; // exclusive
  size_type curPos {0};

  {
    const LO lclRow = 0;
    ptr_h[lclRow] = curPos;
    const GO gblDiagCol = rowMap.getGlobalElement (lclRow);
    ind_h[curPos++] = colMap.getLocalElement (gblDiagCol);
    ind_h[curPos++] = colMap.getLocalElement (gblDiagCol + 1);
  }
  for (LO lclRow = lclStartRow; lclRow < lclEndRow; ++lclRow) {
    ptr_h[lclRow] = curPos;
    const GO gblRow = rowMap.getGlobalElement (lclRow);
    const GO gblDiagCol = gblRow;
    ind_h[curPos++] = colMap.getLocalElement (gblDiagCol - 1);
    ind_h[curPos++] = colMap.getLocalElement (gblDiagCol);
    ind_h[curPos++] = colMap.getLocalElement (gblDiagCol + 1);
  }
  {
    const LO lclRow = lclNumRows - 1;
    ptr_h[lclRow] = curPos;
    const GO gblDiagCol = rowMap.getGlobalElement (lclRow);
    ind_h[curPos++] = colMap.getLocalElement (gblDiagCol - 1);
    ind_h[curPos++] = colMap.getLocalElement (gblDiagCol);
  }
  ptr_h[lclNumRows] = curPos;

  if (curPos != lclNumEnt) {
    throw std::logic_error ("BUG!");
  }

  Kokkos::deep_copy (ptr, ptr_h);
  Kokkos::deep_copy (ind, ind_h);
  return crs_graph_type (ind, ptr);
}

template<class LO, class GO, class DT>
void
testGraph (bool& success,
           Teuchos::FancyOStream& out,
           const LocalCrsGraph<LO, DT>& graph,
           const LocalMap<LO, GO, DT>& lclRowMap,
           const LocalMap<LO, GO, DT>& lclColMap,
           const bool ignoreMapsForTriangularStructure,
           const Tpetra::Details::LocalTriangularStructureResult<LO>& expectedResult)
{
  using Tpetra::Details::determineLocalTriangularStructure;
  using result_type = Tpetra::Details::LocalTriangularStructureResult<LO>;

  const result_type result =
    determineLocalTriangularStructure (graph, lclRowMap, lclColMap,
                                       ignoreMapsForTriangularStructure);
  TEST_EQUALITY( result.diagCount, expectedResult.diagCount );
  TEST_EQUALITY( result.maxNumRowEnt, expectedResult.maxNumRowEnt );
  TEST_EQUALITY( result.couldBeLowerTriangular, expectedResult.couldBeLowerTriangular );
  TEST_EQUALITY( result.couldBeUpperTriangular, expectedResult.couldBeUpperTriangular );
}

TEUCHOS_UNIT_TEST(DetermineLocalTriangularStructure, Test0)
{
  using Tpetra::Details::determineLocalTriangularStructure;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::endl;

  using LO = Tpetra::Map<>::local_ordinal_type;
  using GO = Tpetra::Map<>::global_ordinal_type;
  using NT = Tpetra::Map<>::node_type;
  using DT = Tpetra::Map<>::device_type;
  using map_type = Tpetra::Map<LO, GO, NT>;
  //using local_map_type = LocalMap<LO, GO, DT>;
  //using local_crs_graph_type = LocalCrsGraph<LO, DT>;
  using result_type = Tpetra::Details::LocalTriangularStructureResult<LO>;

  out << "Test Tpetra::Details::determineLocalTriangularStructure" << endl;
  Teuchos::OSTab tab1 (out);
  auto comm = Tpetra::getDefaultComm ();

  // Use an odd number, so that locally reversing the Map doesn't
  // change the diagonal entry.
  constexpr LO lclNumRows = 9;
  const GO gblNumRows = GO (lclNumRows) * GO (comm->getSize ());

  const GO indexBase = 0;
  RCP<const map_type> rowMap_contig =
    rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));
  RCP<const map_type> colMap_contig =
    rcp (new map_type (gblNumRows, lclNumRows, indexBase, comm));

  out << "Test with row Map and column Map both contiguous" << endl;
  {
    Teuchos::OSTab tab2 (out);

    RCP<const map_type> rowMap = rowMap_contig;
    RCP<const map_type> colMap = colMap_contig;
    auto lclRowMap = rowMap->getLocalMap ();
    auto lclColMap = colMap->getLocalMap ();
    {
      out << "Diagonal graph" << endl;
      Teuchos::OSTab tab3 (out);

      auto graph = makeDiagonalGraph<LO, GO, NT> (lclNumRows, *rowMap, *colMap);
      for (bool ignoreMapsForTriangularStructure : {false, true}) {
        const LO expectedLclNumDiag = lclNumRows;
        const LO expectedMaxNumRowEnt = (lclNumRows < LO (1)) ? lclNumRows : LO (1);
        const result_type expectedResult {expectedLclNumDiag,
            expectedMaxNumRowEnt, true, true};
        testGraph<LO, GO, DT> (success, out, graph, lclRowMap, lclColMap,
                               ignoreMapsForTriangularStructure,
                               expectedResult);
      }
    }

    // triangular graphs
    for (bool lowerTriangular : {false, true}) {
      const bool upperTriangular = ! lowerTriangular;
      for (bool explicitDiagonal: {false, true}) {
        out << (lowerTriangular ? "Lower" : "Upper") << " triangular graph, "
            << "with " << (explicitDiagonal ? "explicit" : "implicit unit")
            << " diagonal" << endl;
        Teuchos::OSTab tab3 (out);

        auto graph = makeTriangularGraph<LO, GO, NT> (lclNumRows, *rowMap, *colMap,
                                                      lowerTriangular, explicitDiagonal);
        for (bool ignoreMapsForTriangularStructure : {false, true}) {
          const LO expectedLclNumDiag = explicitDiagonal ? lclNumRows : LO (0);
          const LO expectedMaxNumRowEnt = [&] () {
              if (lclNumRows < LO (1)) {
                return LO (0);
              }
              else {
                return explicitDiagonal ? LO (2) : LO (1);
              }
            } ();
          const result_type expectedResult {expectedLclNumDiag,
              expectedMaxNumRowEnt, lowerTriangular, upperTriangular};
          testGraph<LO, GO, DT> (success, out, graph, lclRowMap, lclColMap,
                                 ignoreMapsForTriangularStructure,
                                 expectedResult);
        }
      }
    }

    {
      out << "Tridiagonal graph" << endl;
      Teuchos::OSTab tab3 (out);

      auto graph = makeTridiagonalGraph<LO, GO, NT> (lclNumRows, *rowMap, *colMap);
      for (bool ignoreMapsForTriangularStructure : {false, true}) {
        const LO expectedLclNumDiag = lclNumRows;
        const LO expectedMaxNumRowEnt = (lclNumRows < LO (3)) ? lclNumRows : LO (3);
        const result_type expectedResult {expectedLclNumDiag,
            expectedMaxNumRowEnt, false, false};
        testGraph<LO, GO, DT> (success, out, graph, lclRowMap, lclColMap,
                               ignoreMapsForTriangularStructure,
                               expectedResult);
      }
    }
  } // row Map and column Map both contiguous

  // In rowMap_noncontig, global row indices are in reverse order ON
  // EACH PROCESS.  On each process, any index in rowMap_noncontig is
  // also in rowMap on that process, and vice versa.  Ditto for
  // colMap_noncontig.
  RCP<const map_type> rowMap_noncontig = [&] () {
    std::vector<GO> rowMapInds (lclNumRows);
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const LO lclRow2 = (lclNumRows - LO (1)) - lclRow;
      rowMapInds[lclRow] = rowMap_contig->getGlobalElement (lclRow2);
    }
    return rcp (new map_type (gblNumRows, rowMapInds.data (),
                              lclNumRows, indexBase, comm));
  } ();
  RCP<const map_type> colMap_noncontig = rowMap_noncontig;

  // Test with rowMap contiguous, and colMap locally reversed (see above).
  {
    RCP<const map_type> rowMap = rowMap_contig;
    RCP<const map_type> colMap = colMap_noncontig;
    auto lclRowMap = rowMap->getLocalMap ();
    auto lclColMap = colMap->getLocalMap ();
    {
      out << "Diagonal graph" << endl;
      Teuchos::OSTab tab3 (out);

      auto graph = makeDiagonalGraph<LO, GO, NT> (lclNumRows, *rowMap, *colMap);
      for (bool ignoreMapsForTriangularStructure : {false, true}) {
        const LO expectedLclNumDiag = lclNumRows;
        const LO expectedMaxNumRowEnt = (lclNumRows < LO (1)) ? lclNumRows : LO (1);
        // Graph is correctly diagonal with respect to global indices.
        // If ignoring global indices, it's neither upper nor lower
        // triangular.
        const bool looksDiag = ! ignoreMapsForTriangularStructure;
        const result_type expectedResult {expectedLclNumDiag,
            expectedMaxNumRowEnt, looksDiag, looksDiag};
        testGraph<LO, GO, DT> (success, out, graph, lclRowMap, lclColMap,
                               ignoreMapsForTriangularStructure,
                               expectedResult);
      }
    }

    // triangular graphs
    for (bool lowerTriangular : {false, true}) {
      const bool upperTriangular = ! lowerTriangular;
      for (bool explicitDiagonal: {false, true}) {
        out << (lowerTriangular ? "Lower" : "Upper") << " triangular graph, "
            << "with " << (explicitDiagonal ? "explicit" : "implicit unit")
            << " diagonal" << endl;
        Teuchos::OSTab tab3 (out);

        auto graph = makeTriangularGraph<LO, GO, NT> (lclNumRows, *rowMap, *colMap,
                                                      lowerTriangular, explicitDiagonal);
        for (bool ignoreMapsForTriangularStructure : {false, true}) {
          const LO expectedLclNumDiag = explicitDiagonal ? lclNumRows : LO (0);
          const LO expectedMaxNumRowEnt = [&] () {
              if (lclNumRows < LO (1)) {
                return LO (0);
              }
              else {
                return explicitDiagonal ? LO (2) : LO (1);
              }
            } ();
          // If ignoring global indices, it's neither upper nor lower
          // triangular.
          const bool expectLowerTri =
            ignoreMapsForTriangularStructure ? false : lowerTriangular;
          const bool expectUpperTri =
            ignoreMapsForTriangularStructure ? false : upperTriangular;
          const result_type expectedResult {expectedLclNumDiag,
              expectedMaxNumRowEnt, expectLowerTri, expectUpperTri};
          testGraph<LO, GO, DT> (success, out, graph, lclRowMap, lclColMap,
                                 ignoreMapsForTriangularStructure,
                                 expectedResult);
        }
      }
    }

    {
      out << "Tridiagonal graph" << endl;
      Teuchos::OSTab tab3 (out);

      auto graph = makeTridiagonalGraph<LO, GO, NT> (lclNumRows, *rowMap, *colMap);
      for (bool ignoreMapsForTriangularStructure : {false, true}) {
        const LO expectedLclNumDiag = lclNumRows;
        const LO expectedMaxNumRowEnt = (lclNumRows < LO (3)) ? lclNumRows : LO (3);
        const result_type expectedResult {expectedLclNumDiag,
            expectedMaxNumRowEnt, false, false};
        testGraph<LO, GO, DT> (success, out, graph, lclRowMap, lclColMap,
                               ignoreMapsForTriangularStructure,
                               expectedResult);
      }
    }
  }
}

} // namespace (anonymous)

