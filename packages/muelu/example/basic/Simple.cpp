// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <cstdio>
#include <iomanip>
#include <iostream>
#include <unistd.h>

#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Xpetra
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>

#include <MueLu_BaseClass.hpp>
#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif
#include <MueLu_Level.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
#include <MatrixLoad.hpp>
#include <DriverCore.hpp>

/*********************************************************************/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::coordinateType real_type;
  typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

  std::string xmlFileName = "";
  clp.setOption("xml", &xmlFileName, "read parameters from an xml file");
  std::string yamlFileName = "";
  clp.setOption("yaml", &yamlFileName, "read parameters from a yaml file");
  bool printTimings = true;
  clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
  std::string timingsFormat = "table-fixed";
  clp.setOption("time-format", &timingsFormat, "timings format (table-fixed | table-scientific | yaml)");
  std::string solveType = "belos";
  clp.setOption("solver", &solveType, "solve type: (none | belos)");
  std::string belosType = "cg";
  clp.setOption("belosType", &belosType, "belos solver type: (Pseudoblock CG | Block CG | Pseudoblock GMRES | Block GMRES | ...) see BelosSolverFactory.hpp for exhaustive list of solvers");
  bool computeCondEst = false;
  clp.setOption("condEst", "noCondEst", &computeCondEst, "compute condition number estimate (currently only available for Pseudoblock CG)");
  bool enforceBoundaryConditionsOnInitialGuess = true;
  clp.setOption("enforceBCs", "noEnforceBCs", &enforceBoundaryConditionsOnInitialGuess, "enforce Dirichlet boundary condition on initial guess");
  double tol = 1e-12;
  clp.setOption("tol", &tol, "solver convergence tolerance");
  bool binaryFormat = false;
  clp.setOption("binary", "ascii", &binaryFormat, "read matrices in binary format");
  std::string rowMapFile;
  clp.setOption("rowmap", &rowMapFile, "map data file");
  std::string colMapFile;
  clp.setOption("colmap", &colMapFile, "colmap data file");
  std::string domainMapFile;
  clp.setOption("domainmap", &domainMapFile, "domainmap data file");
  std::string rangeMapFile;
  clp.setOption("rangemap", &rangeMapFile, "rangemap data file");
  std::string matrixFile;
  clp.setOption("matrix", &matrixFile, "matrix data file");
  std::string rhsFile;
  clp.setOption("rhs", &rhsFile, "rhs data file");
  std::string coordFile;
  clp.setOption("coords", &coordFile, "coordinates data file");
  std::string coordMapFile;
  clp.setOption("coordsmap", &coordMapFile, "coordinates map data file");
  std::string nullFile;
  clp.setOption("nullspace", &nullFile, "nullspace data file");
  std::string materialFile;
  clp.setOption("material", &materialFile, "material data file");
  int maxIts = 200;
  clp.setOption("its", &maxIts, "maximum number of solver iterations");
  int numVectors = 1;
  clp.setOption("multivector", &numVectors, "number of rhs to solve simultaneously");
  bool scaleResidualHist = true;
  clp.setOption("scale", "noscale", &scaleResidualHist, "scaled Krylov residual history");
  bool solvePreconditioned = true;
  clp.setOption("solve-preconditioned", "no-solve-preconditioned", &solvePreconditioned, "use MueLu preconditioner in solve");
  int cacheSize = 0;
  clp.setOption("cachesize", &cacheSize, "cache size (in KB)");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName != "" && yamlFileName != "", std::runtime_error,
                             "Cannot provide both xml and yaml input files");

  // Instead of checking each time for rank, create a rank 0 stream
  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream &out       = *fancy;
  out.setOutputToRootOnly(0);

  ParameterList paramList;
  auto inst = xpetraParameters.GetInstantiation();

  if (yamlFileName != "") {
    Teuchos::updateParametersFromYamlFileAndBroadcast(yamlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  } else {
    xmlFileName = (xmlFileName != "" ? xmlFileName : "simple.xml");
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);
  }

  if (inst == Xpetra::COMPLEX_INT_INT && solveType == "belos") {
    belosType = "gmres";
    out << "WARNING: CG will not work with COMPLEX scalars, switching to GMRES" << std::endl;
  }

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  std::ostringstream galeriStream;

  comm->barrier();
  Teuchos::TimeMonitor::setStackedTimer(Teuchos::null);
  RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: S - Global Time")));
  RCP<TimeMonitor> tm                = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 1 - Matrix Build")));

  RCP<Matrix> A;
  RCP<const Map> map;
  RCP<RealValuedMultiVector> coordinates;
  RCP<MultiVector> nullspace;
  RCP<MultiVector> material;
  RCP<MultiVector> X, B;

  // Load the matrix off disk (or generate it via Galeri)
  MatrixLoad<SC, LO, GO, NO>(comm, lib, binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile, domainMapFile, rangeMapFile, coordFile, coordMapFile, nullFile, materialFile, map, A, coordinates, nullspace, material, X, B, numVectors, galeriParameters, xpetraParameters, galeriStream);
  comm->barrier();
  tm = Teuchos::null;
  out << galeriStream.str();

  // =========================================================================
  // Preconditioner construction
  // =========================================================================
  bool useML = paramList.isParameter("use external multigrid package") && (paramList.get<std::string>("use external multigrid package") == "ml");
  out << "*********** MueLu ParameterList ***********" << std::endl;
  out << paramList;
  out << "*******************************************" << std::endl;

  RCP<Hierarchy> H;
  RCP<Operator> Prec;
  {
    comm->barrier();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Driver: 2 - MueLu Setup")));
    PreconditionerSetup(A, coordinates, nullspace, material, paramList, false, false, useML, false, 0, H, Prec);
    comm->barrier();
    tm = Teuchos::null;
  }

  // =========================================================================
  // System solution (Ax = b)
  // =========================================================================
  {
    comm->barrier();
    SystemSolve(A, X, B, H, Prec, out, solveType, belosType, false, false, useML, cacheSize, 0, scaleResidualHist, solvePreconditioned, maxIts, tol, computeCondEst, enforceBoundaryConditionsOnInitialGuess);
    comm->barrier();
  }

  tm                = Teuchos::null;
  globalTimeMonitor = Teuchos::null;

  if (printTimings) {
    RCP<ParameterList> reportParams = rcp(new ParameterList);
    if (timingsFormat == "yaml") {
      reportParams->set("Report format", "YAML");  // "Table" or "YAML"
      reportParams->set("YAML style", "compact");  // "spacious" or "compact"
    }
    reportParams->set("How to merge timer sets", "Union");
    reportParams->set("alwaysWriteLocal", false);
    reportParams->set("writeGlobalStats", true);
    reportParams->set("writeZeroTimers", false);

    const std::string filter = "";

    std::ios_base::fmtflags ff(out.flags());
    if (timingsFormat == "table-fixed")
      out << std::fixed;
    else
      out << std::scientific;
    TimeMonitor::report(comm.ptr(), out, filter, reportParams);
    out << std::setiosflags(ff);
  }

  TimeMonitor::clearCounters();
  out << std::endl;

  fflush(NULL);
  comm->barrier();

  return EXIT_SUCCESS;
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
