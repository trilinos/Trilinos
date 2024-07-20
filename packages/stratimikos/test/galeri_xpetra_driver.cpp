// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

/*
   Call Ifpack2 via the Stratimikos interface.

Usage:
./Ifpack2_Stratimikos.exe : use xml configuration file stratimikos_ParameterList.xml

Note:
The source code is not MueLu specific and can be used with any Stratimikos strategy.
*/

// Teuchos includes
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_StackedTimer.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_YamlParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include "Teuchos_AbstractFactoryStd.hpp"

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>

// Stratimikos includes
#include <Stratimikos_LinearSolverBuilder.hpp>
#include <Stratimikos_InternalConfig.h>

// Xpetra include
#include <Xpetra_Parameters.hpp>
#include <Xpetra_ThyraUtils.hpp>

// Galeri includes
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>


template <class Scalar>
int
main_(int argc, char *argv[], Teuchos::CommandLineProcessor& clp) {
  typedef Tpetra::Map<> map_type;
  typedef map_type::local_ordinal_type LocalOrdinal;
  typedef map_type::global_ordinal_type GlobalOrdinal;
  typedef map_type::node_type Node;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::TimeMonitor;
  #include <Xpetra_UseShortNames.hpp>

  bool success = false;
  bool verbose = true;
  try {

    //
    // MPI initialization
    //
    // Teuchos::CommandLineProcessor clp(false);
    const auto comm = Teuchos::DefaultComm<int>::getComm ();

    //
    // Parameters
    //
    // manage parameters of the test case
    Galeri::Xpetra::Parameters<GlobalOrdinal> galeriParameters(clp, 100, 100, 100, "Laplace2D");
    // manage parameters of Xpetra
    Xpetra::Parameters                        xpetraParameters(clp);

    // command line parameters
    std::string xmlFileName       = "stratimikos_ParameterList.xml"; clp.setOption("xml",      &xmlFileName,       "read parameters from an xml file");
    std::string yamlFileName      = "";                 clp.setOption("yaml",                  &yamlFileName,      "read parameters from a yaml file");
    bool        printTimings      = false;              clp.setOption("timings", "notimings",  &printTimings,      "print timings to screen");
    bool        use_stacked_timer = false;              clp.setOption("stacked-timer", "no-stacked-timer", &use_stacked_timer, "Run with or without stacked timer output");
    std::string timingsFormat     = "table-fixed";      clp.setOption("time-format",           &timingsFormat,     "timings format (table-fixed | table-scientific | yaml)");
    int         numVectors        = 1;                  clp.setOption("multivector",           &numVectors,        "number of rhs to solve simultaneously");
    int         numSolves         = 1;                  clp.setOption("numSolves",             &numSolves,         "number of times the system should be solved");

    switch (clp.parse(argc,argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *fancy;
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

    RCP<Matrix>                A;
    RCP<const Map>             map;
    RCP<MultiVector>           X, B;

    std::ostringstream galeriStream;
    Teuchos::ParameterList galeriList = galeriParameters.GetParameterList();
    galeriStream << "========================================================\n" << xpetraParameters;
    galeriStream << galeriParameters;

    // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
    //                                 d1  d2  d3
    //                                 d4  d5  d6
    //                                 d7  d8  d9
    //                                 d10 d11 d12
    // A perfect distribution is only possible when the #processors is a perfect square.
    // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
    // size. For example, np=14 will give a 7-by-2 distribution.
    // If you don't want Galeri to do this, specify mx or my on the galeriList.
    std::string matrixType = galeriParameters.GetMatrixType();

    // Create map
    if (matrixType == "Laplace1D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);

    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" ||
               matrixType == "BigStar2D" || matrixType == "AnisotropicDiffusion" || matrixType == "Elasticity2D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);

    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 2);
    if (matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, 3);

    galeriStream << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
                 << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
                 << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
                 << "========================================================" << std::endl;

    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
      // Our default test case for elasticity: all boundaries of a square/cube have Neumann b.c. except left which has Dirichlet
      galeriList.set("right boundary" , "Neumann");
      galeriList.set("bottom boundary", "Neumann");
      galeriList.set("top boundary"   , "Neumann");
      galeriList.set("front boundary" , "Neumann");
      galeriList.set("back boundary"  , "Neumann");
    }

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
    A = Pr->BuildMatrix();

    if (matrixType == "Elasticity2D" ||
        matrixType == "Elasticity3D") {
      A->SetFixedBlockSize((galeriParameters.GetMatrixType() == "Elasticity2D") ? 2 : 3);
    }

    out << galeriStream.str();
    X = MultiVectorFactory::Build(map, numVectors);
    B = MultiVectorFactory::Build(map, numVectors);
    B->putScalar(1);
    X->putScalar(0);

    //
    // Build Thyra linear algebra objects
    //

    RCP<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xpCrsA = Teuchos::rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(A);

    RCP<const Thyra::LinearOpBase<Scalar> >    thyraA = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(xpCrsA->getCrsMatrix());
    RCP<      Thyra::MultiVectorBase<Scalar> > thyraX = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<Scalar> >(Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraMultiVector(X));
    RCP<const Thyra::MultiVectorBase<Scalar> > thyraB = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyraMultiVector(B);

    //
    // Build Stratimikos solver
    //

    // This is the Stratimikos main class (= factory of solver factory).
    Stratimikos::LinearSolverBuilder<Scalar> linearSolverBuilder;

    // Setup solver parameters using a Stratimikos parameter list.
    linearSolverBuilder.setParameterList(paramList);

    // Build a new "solver factory" according to the previously specified parameter list.
    RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
    auto precFactory = solverFactory->getPreconditionerFactory();
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
          reportParams->set("Report format",             "YAML");            // "Table" or "YAML"
          reportParams->set("YAML style",                "compact");         // "spacious" or "compact"
        }
        reportParams->set("How to merge timer sets",   "Union");
        reportParams->set("alwaysWriteLocal",          false);
        reportParams->set("writeGlobalStats",          true);
        reportParams->set("writeZeroTimers",           false);

        const std::string filter = "";

        std::ios_base::fmtflags ff(out.flags());
        if (timingsFormat == "table-fixed") out << std::fixed;
        else                                out << std::scientific;
        TimeMonitor::report(comm.ptr(), out, filter, reportParams);
        out << std::setiosflags(ff);
      }
    }

    TimeMonitor::clearCounters();
    out << std::endl;

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}


enum scalarType {
  DOUBLE,
  FLOAT,
  COMPLEX_DOUBLE,
  COMPLEX_FLOAT
};

int
main(int argc, char *argv[]) {
  Teuchos::GlobalMPISession session (&argc, &argv, NULL);

  Teuchos::CommandLineProcessor clp(false);
  scalarType scalar = DOUBLE;
  std::vector<const char*> availableScalarTypeStrings;
  std::vector<scalarType> availableScalarTypes;
#ifdef HAVE_TPETRA_INST_DOUBLE
  availableScalarTypeStrings.push_back("double");
  availableScalarTypes.push_back(DOUBLE);
#endif
#ifdef HAVE_TPETRA_INST_FLOAT
  availableScalarTypeStrings.push_back("float");
  availableScalarTypes.push_back(FLOAT);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  availableScalarTypeStrings.push_back("complex<double>");
  availableScalarTypes.push_back(COMPLEX_DOUBLE);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  availableScalarTypeStrings.push_back("complex<float>");
  availableScalarTypes.push_back(COMPLEX_FLOAT);
#endif
  clp.setOption("scalarType", &scalar, availableScalarTypes.size(), availableScalarTypes.data(), availableScalarTypeStrings.data(), "scalar type");
  clp.recogniseAllOptions(false);
  switch (clp.parse(argc, argv, NULL)) {
    case Teuchos::CommandLineProcessor::PARSE_ERROR:                return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:         break;
  }

#ifdef HAVE_TPETRA_INST_DOUBLE
  if (scalar == DOUBLE)
    return main_<double>(argc, argv, clp);
#endif
#ifdef HAVE_TPETRA_INST_FLOAT
  if (scalar == FLOAT)
    return main_<float>(argc, argv, clp);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  if (scalar == COMPLEX_DOUBLE)
    return main_<std::complex<double> >(argc, argv, clp);
#endif
#ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
  if (scalar == COMPLEX_FLOAT)
    return main_<std::complex<float> >(argc, argv, clp);
#endif
}
