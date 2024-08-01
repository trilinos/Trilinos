// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

/*
   Call MueLu via the Stratimikos interface.

Usage:
./MueLu_Stratimikos.exe : use xml configuration file stratimikos_ParameterList.xml

Note:
The source code is not MueLu specific and can be used with any Stratimikos strategy.
*/

// Teuchos includes
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>

// Stratimikos includes
#include <Stratimikos_LinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>

// Xpetra include
#include <Xpetra_Parameters.hpp>

// MueLu includes
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include <MatrixLoad.hpp>

// Galeri includes
#include <Galeri_XpetraParameters.hpp>

template <typename Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename STS::coordinateType real_type;
  typedef Xpetra::MultiVector<real_type, LocalOrdinal, GlobalOrdinal, Node> RealValuedMultiVector;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
  try {
    //
    // MPI initialization
    //
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    //
    // Parameters
    //
    // manage parameters of the test case
    Galeri::Xpetra::Parameters<GlobalOrdinal> matrixParameters(clp, 100, 100, 100, "Laplace2D");
    // manage parameters of Xpetra
    Xpetra::Parameters xpetraParameters(clp);

    // command line parameters
    std::string xmlFileName = "stratimikos_ParameterList.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from an xml file");
    std::string yamlFileName = "";
    clp.setOption("yaml", &yamlFileName, "read parameters from a yaml file");
    bool printTimings = false;
    clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
    bool use_stacked_timer = false;
    clp.setOption("stacked-timer", "no-stacked-timer", &use_stacked_timer, "Run with or without stacked timer output");
    std::string timingsFormat = "table-fixed";
    clp.setOption("time-format", &timingsFormat, "timings format (table-fixed | table-scientific | yaml)");
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
    int numVectors = 1;
    clp.setOption("multivector", &numVectors, "number of rhs to solve simultaneously");
    int numSolves = 1;
    clp.setOption("numSolves", &numSolves, "number of times the system should be solved");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream &out       = *fancy;
    out.setOutputToRootOnly(0);

    // Set up timers
    Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
    if (use_stacked_timer)
      stacked_timer = rcp(new Teuchos::StackedTimer("Main"));
    TimeMonitor::setStackedTimer(stacked_timer);

    // Read in parameter list
    TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName == "" && yamlFileName == "", std::runtime_error,
                               "Need to provide xml or yaml input file");
    RCP<ParameterList> paramList = rcp(new ParameterList("params"));
    if (yamlFileName != "")
      Teuchos::updateParametersFromYamlFileAndBroadcast(yamlFileName, paramList.ptr(), *comm);
    else
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, paramList.ptr(), *comm);

    //
    // Construct the problem
    //

    RCP<Matrix> A;
    RCP<const Map> map;
    RCP<RealValuedMultiVector> coordinates;
    RCP<MultiVector> nullspace;
    RCP<MultiVector> material;
    RCP<MultiVector> X, B;

    std::ostringstream galeriStream;
    MatrixLoad<SC, LocalOrdinal, GlobalOrdinal, Node>(comm, lib, binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile, domainMapFile, rangeMapFile, coordFile, coordMapFile, nullFile, materialFile, map, A, coordinates, nullspace, material, X, B, numVectors, matrixParameters, xpetraParameters, galeriStream);
    out << galeriStream.str();
    X->putScalar(0);

    //
    // Build Thyra linear algebra objects
    //

    RCP<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> > xpCrsA = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(A);

    RCP<const Thyra::LinearOpBase<Scalar> > thyraA    = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(xpCrsA->getCrsMatrix());
    RCP<Thyra::MultiVectorBase<Scalar> > thyraX       = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraMultiVector(X));
    RCP<const Thyra::MultiVectorBase<Scalar> > thyraB = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraMultiVector(B);

    //
    // Build Stratimikos solver
    //

    // This is the Stratimikos main class (= factory of solver factory).
    Stratimikos::LinearSolverBuilder<Scalar> linearSolverBuilder;
    // Register MueLu as a Stratimikos preconditioner strategy.
    Stratimikos::enableMueLu<Scalar, LocalOrdinal, GlobalOrdinal, Node>(linearSolverBuilder);

    // add coordinates and nullspace to parameter list
    if (paramList->isSublist("Preconditioner Types") &&
        paramList->sublist("Preconditioner Types").isSublist("MueLu")) {
      ParameterList &userParamList = paramList->sublist("Preconditioner Types").sublist("MueLu").sublist("user data");
      userParamList.set<RCP<RealValuedMultiVector> >("Coordinates", coordinates);
      userParamList.set<RCP<MultiVector> >("Nullspace", nullspace);
    }

    // Setup solver parameters using a Stratimikos parameter list.
    linearSolverBuilder.setParameterList(paramList);

    // Build a new "solver factory" according to the previously specified parameter list.
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    auto precFactory                                                = solverFactory->getPreconditionerFactory();
    RCP<Thyra::PreconditionerBase<Scalar> > prec;
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<Scalar> > thyraInverseA;
    if (!precFactory.is_null()) {
      prec = precFactory->createPrec();

      // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver.
      Thyra::initializePrec<Scalar>(*precFactory, thyraA, prec.ptr());
      thyraInverseA = solverFactory->createOp();
      Thyra::initializePreconditionedOp<Scalar>(*solverFactory, thyraA, prec, thyraInverseA.ptr());
    } else {
      thyraInverseA = Thyra::linearOpWithSolve(*solverFactory, thyraA);
    }

    // Solve Ax = b.
    Thyra::SolveStatus<Scalar> status = Thyra::solve<Scalar>(*thyraInverseA, Thyra::NOTRANS, *thyraB, thyraX.ptr());

    success = (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);

    for (int solveno = 1; solveno < numSolves; solveno++) {
      if (!precFactory.is_null())
        Thyra::initializePrec<Scalar>(*precFactory, thyraA, prec.ptr());
      thyraX->assign(0.);

      status = Thyra::solve<Scalar>(*thyraInverseA, Thyra::NOTRANS, *thyraB, thyraX.ptr());

      success = success && (status.solveStatus == Thyra::SOLVE_STATUS_CONVERGED);
    }

    // print timings
    if (printTimings) {
      if (use_stacked_timer) {
        stacked_timer->stop("Main");
        Teuchos::StackedTimer::OutputOptions options;
        options.output_fraction = options.output_histogram = options.output_minmax = true;
        stacked_timer->report(out, comm, options);
      } else {
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
    }

    TimeMonitor::clearCounters();
    out << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
