// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Build a Tpetra::CrsGraph for a two-dimensional five-point finite-difference
// stencil on a rectangular grid using Tpetra::Experimental::createCrsGraphFromStencil

#include "Tpetra_Core.hpp"
#include "Tpetra_CreateCrsGraphFromStencil.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>

int main(int argc, char* argv[]) {
  using map_type            = Tpetra::Map<>;
  using local_ordinal_type  = typename map_type::local_ordinal_type;
  using global_ordinal_type = typename map_type::global_ordinal_type;
  Tpetra::ScopeGuard scope(&argc, &argv);
  const auto comm = Tpetra::getDefaultComm();
  auto out        = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  global_ordinal_type nx = 20;
  global_ordinal_type ny = 20;
  bool verbose           = false;

  Teuchos::CommandLineProcessor commandLine(false, true);
  commandLine.setOption("nx", &nx, "Number of grid points in x");
  commandLine.setOption("ny", &ny, "Number of grid points in y");
  commandLine.setOption("verbose", "quiet", &verbose,
                        "Describe the completed graph");

  const auto parseResult = commandLine.parse(argc, argv);
  if (parseResult == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    return EXIT_SUCCESS;
  if (parseResult != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    return EXIT_FAILURE;

  int localSuccess = 1;
  try {
    TEUCHOS_TEST_FOR_EXCEPTION(nx <= 0 || ny <= 0, std::invalid_argument,
                               "nx and ny must both be positive");
    TEUCHOS_TEST_FOR_EXCEPTION(
        nx > std::numeric_limits<global_ordinal_type>::max() / ny,
        std::overflow_error, "nx * ny does not fit in the global ordinal type");

    const Tpetra::global_size_t numGlobalRows =
        static_cast<Tpetra::global_size_t>(nx) *
        static_cast<Tpetra::global_size_t>(ny);

    // Use the same two-dimensional Cartesian process grid inferred by
    // createCrsGraphFromStencil.  This example keeps the decomposition uniform
    // so that Tpetra's contiguous map gives every rank one rectangular block.
    int procsX = 1;
    while ((procsX + 1) * (procsX + 1) <= comm->getSize())
      ++procsX;
    int procsY = comm->getSize() / procsX;
    while (procsX * procsY != comm->getSize())
      procsY = comm->getSize() / (--procsX);
    TEUCHOS_TEST_FOR_EXCEPTION(
        nx % procsX != 0 || ny % procsY != 0, std::invalid_argument,
        "nx and ny must be divisible by the " << procsX << " x " << procsY
                                              << " process grid");

    const global_ordinal_type indexBase = 0;
    const Teuchos::RCP<const map_type> rowMap =
        Teuchos::rcp(new map_type(numGlobalRows, indexBase, comm));
    const Teuchos::Array<local_ordinal_type> localGridDimensions{
        static_cast<local_ordinal_type>(nx / procsX),
        static_cast<local_ordinal_type>(ny / procsY)};
    const Teuchos::Array<Tpetra::Experimental::StencilOffset> stencilOffsets{
        {0, 0, 0}, {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}};
    const local_ordinal_type dofsPerNode = 1;
    const bool debug                     = true;
    const auto graph                     = Tpetra::Experimental::createCrsGraphFromStencil(
                            rowMap, localGridDimensions(), stencilOffsets(), dofsPerNode,
                            Tpetra::Experimental::StencilNumbering::rankBlocked, debug);

    // A rectangular nx-by-ny grid has nx*ny diagonal entries, two directed
    // entries per horizontal edge, and two per vertical edge.
    const Tpetra::global_size_t expectedEntries =
        numGlobalRows +
        2 * static_cast<Tpetra::global_size_t>(nx - 1) *
            static_cast<Tpetra::global_size_t>(ny) +
        2 * static_cast<Tpetra::global_size_t>(ny - 1) *
            static_cast<Tpetra::global_size_t>(nx);
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph->getGlobalNumEntries() != expectedEntries, std::runtime_error,
        "Graph has " << graph->getGlobalNumEntries() << " entries, expected "
                     << expectedEntries);

    const global_ordinal_type localNx = localGridDimensions[0];
    const global_ordinal_type localNy = localGridDimensions[1];
    const int procX                   = comm->getRank() % procsX;
    const int procY                   = comm->getRank() / procsX;
    for (size_t localRow = 0; localRow < rowMap->getLocalNumElements();
         ++localRow) {
      const auto row = rowMap->getGlobalElement(
          static_cast<local_ordinal_type>(localRow));
      const global_ordinal_type localX = localRow % localNx;
      const global_ordinal_type localY = localRow / localNx;
      const global_ordinal_type x      = procX * localNx + localX;
      const global_ordinal_type y      = procY * localNy + localY;
      const size_t expectedRowLength =
          1 + static_cast<size_t>(x > 0) + static_cast<size_t>(x + 1 < nx) +
          static_cast<size_t>(y > 0) + static_cast<size_t>(y + 1 < ny);
      TEUCHOS_TEST_FOR_EXCEPTION(
          graph->getNumEntriesInGlobalRow(row) != expectedRowLength,
          std::runtime_error, "Incorrect stencil graph row " << row);
    }

    if (verbose)
      graph->describe(*out, Teuchos::VERB_EXTREME);
    *out << "Constructed a " << nx << " x " << ny
         << " five-point finite-difference graph with " << expectedEntries
         << " entries on " << comm->getSize() << " rank(s).\n";
  } catch (const std::exception& error) {
    std::cerr << "Rank " << comm->getRank() << ": " << error.what() << '\n';
    localSuccess = 0;
  }

  int globalSuccess = 0;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, localSuccess,
                     Teuchos::outArg(globalSuccess));
  if (globalSuccess) {
    *out << "End Result: TEST PASSED\n";
    return EXIT_SUCCESS;
  }
  *out << "End Result: TEST FAILED\n";
  return EXIT_FAILURE;
}
