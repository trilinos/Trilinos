// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Galeri_XpetraMaps.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using scalar_type = double;
using map_type    = Tpetra::Map<>;
using matrix_type = Tpetra::CrsMatrix<>;
using vector_type = Tpetra::MultiVector<>;
using LO          = typename map_type::local_ordinal_type;
using GO          = typename map_type::global_ordinal_type;

int run(int argc, char* argv[],
        const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
  int nx = 120;
  int ny = 120;
  int nz = 120;
  std::string watchrProblemName =
      "Galeri scalar 7-point stencil Galeri implementation";

  Teuchos::CommandLineProcessor clp(false);
  clp.setOption("nx", &nx, "global number of grid points in x");
  clp.setOption("ny", &ny, "global number of grid points in y");
  clp.setOption("nz", &nz, "global number of grid points in z");
  clp.setOption("watchr-problem-name", &watchrProblemName,
                "problem name for Watchr plot headers");
  clp.setDocString(
      "Measure Galeri construction of a scalar three-dimensional "
      "seven-point stencil matrix.");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
      break;
  }
  if (nx <= 0 || ny <= 0 || nz <= 0)
    throw std::invalid_argument("nx, ny, and nz must be positive");

  Teuchos::ParameterList parameters;
  parameters.set("nx", static_cast<GO>(nx));
  parameters.set("ny", static_cast<GO>(ny));
  parameters.set("nz", static_cast<GO>(nz));
  Teuchos::RCP<const map_type> rowMap =
      Galeri::Xpetra::CreateMap<LO, GO, map_type>(
          "Cartesian3D", comm, parameters);
  auto problem = Galeri::Xpetra::BuildProblem<
      scalar_type, LO, GO, map_type, matrix_type, vector_type>(
      "Laplace3D", rowMap, parameters);

  const std::string rootTimerName =
      "Galeri CreateCrsGraphFromStencil performance";
  auto stackedTimer =
      Teuchos::rcp(new Teuchos::StackedTimer(rootTimerName.c_str()));
  Teuchos::TimeMonitor::setStackedTimer(stackedTimer);

  Teuchos::RCP<matrix_type> matrix;
  comm->barrier();
  Kokkos::fence();
  {
    Teuchos::TimeMonitor timer(*Teuchos::TimeMonitor::getNewTimer(
        "Galeri: full matrix construction"));
    matrix = problem->BuildMatrix();
    Kokkos::fence();
  }

  stackedTimer->stop(rootTimerName.c_str());

  const Tpetra::global_size_t expectedRows =
      static_cast<Tpetra::global_size_t>(nx) *
      static_cast<Tpetra::global_size_t>(ny) *
      static_cast<Tpetra::global_size_t>(nz);
  if (matrix.is_null() || !matrix->isFillComplete() ||
      matrix->getGlobalNumRows() != expectedRows)
    throw std::runtime_error("Galeri did not construct the expected matrix");

  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction   = options.output_histogram =
      options.output_minmax = true;
  stackedTimer->report(std::cout, comm, options);
  auto xmlOut = stackedTimer->reportWatchrXML(watchrProblemName, comm);
  if (!xmlOut.empty() && comm->getRank() == 0)
    std::cout << "\nAlso created Watchr performance report " << xmlOut
              << '\n';

  if (comm->getRank() == 0)
    std::cout << "End Result: TEST PASSED\n";
  return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scope(&argc, &argv);
  const auto comm = Tpetra::getDefaultComm();
  try {
    return run(argc, argv, comm);
  } catch (const std::exception& error) {
    std::cerr << "Rank " << comm->getRank() << ": " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
