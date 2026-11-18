// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <iostream>
#include <string>

namespace TpetraExamples {

template <typename V>
Kokkos::View<typename V::data_type, typename V::array_layout, typename V::device_type, Kokkos::MemoryTraits<Kokkos::Unmanaged>>
makeUnmanaged(const V& v) {
  return v;
}

// Return a pointer (RCP is like std::shared_ptr) to an output
// stream.  It prints on Process 0 of the given MPI communicator,
// but ignores all output on other MPI processes.
Teuchos::RCP<Teuchos::FancyOStream>
getOutputStream(const Teuchos::Comm<int>& comm) {
  using Teuchos::getFancyOStream;

  const int myRank = comm.getRank();
  if (0 == myRank) {
    // Process 0 of the given communicator prints to std::cout.
    return getFancyOStream(Teuchos::rcpFromRef(std::cout));
  } else {
    // A "black hole output stream" ignores all output directed to it.
    return getFancyOStream(Teuchos::rcp(new Teuchos::oblackholestream()));
  }
}

// Check that two maps have the same global indices, in the same local order,
// on the calling process.  Returns true on a local match.
template <class MapType>
bool sameLocalGIDs(const MapType& a, const MapType& b) {
  if (a.getLocalNumElements() != b.getLocalNumElements()) return false;
  auto aGIDs = a.getMyGlobalIndices();
  auto bGIDs = b.getMyGlobalIndices();
  for (size_t i = 0; i < aGIDs.extent(0); ++i)
    if (aGIDs(i) != bGIDs(i)) return false;
  return true;
}

// Check that two fill-complete CrsGraphs are identical on the calling process:
// same row, column, domain and range maps (column map compared GID by GID, in
// local order), and the same local row offsets and local column indices, in
// the same order within each row.  Prints the first difference found to
// std::cerr.  Returns true on a local match.
template <class GraphType>
bool identicalLocalGraphs(const GraphType& ref, const GraphType& test,
                          const std::string& label, const int myRank) {
  auto report = [&](const std::string& what) {
    std::cerr << "Proc " << myRank << ": " << label << ": " << what << std::endl;
    return false;
  };
  if (!sameLocalGIDs(*ref.getRowMap(), *test.getRowMap())) return report("row maps differ");
  if (!sameLocalGIDs(*ref.getColMap(), *test.getColMap())) return report("column maps differ");
  if (!sameLocalGIDs(*ref.getDomainMap(), *test.getDomainMap())) return report("domain maps differ");
  if (!sameLocalGIDs(*ref.getRangeMap(), *test.getRangeMap())) return report("range maps differ");

  auto refGraph  = ref.getLocalGraphHost();
  auto testGraph = test.getLocalGraphHost();
  if (refGraph.row_map.extent(0) != testGraph.row_map.extent(0)) return report("number of rows differs");
  for (size_t i = 0; i < refGraph.row_map.extent(0); ++i) {
    if (refGraph.row_map(i) != testGraph.row_map(i))
      return report("row offsets differ at local row " + std::to_string(i));
  }
  // The row offsets are equal, so only the first row_map(numRows) entries are
  // meaningful (entries may be over-allocated in unpacked storage).
  const size_t nnz = refGraph.row_map.extent(0) ? refGraph.row_map(refGraph.row_map.extent(0) - 1) : 0;
  for (size_t k = 0; k < nnz; ++k) {
    if (refGraph.entries(k) != testGraph.entries(k))
      return report("column indices differ at local entry " + std::to_string(k));
  }
  return true;
}

// Check that two fully-assembled FECrsGraphs (in owned mode) are identical,
// both their owned graphs and their owned+shared graphs.  Collective over comm.
// Returns true iff the graphs match on every process.
template <class FEGraphType>
bool identicalFEGraphs(FEGraphType& ref, FEGraphType& test, const Teuchos::Comm<int>& comm) {
  const int myRank = comm.getRank();
  bool ok          = identicalLocalGraphs(ref, test, "owned graph", myRank);
  // Look at the owned+shared graphs, then switch both back to owned mode.
  ref.switchActiveCrsGraph();
  test.switchActiveCrsGraph();
  ok = identicalLocalGraphs(ref, test, "owned+shared graph", myRank) && ok;
  ref.switchActiveCrsGraph();
  test.switchActiveCrsGraph();

  int lclOK = ok ? 1 : 0;
  int gblOK = 0;
  Teuchos::reduceAll<int, int>(comm, Teuchos::REDUCE_MIN, lclOK, Teuchos::outArg(gblOK));
  return gblOK == 1;
}

}  // namespace TpetraExamples

#endif  // TPETRAEXAMPLES_FEM_ASSEMBLY_UTILITY_HPP
