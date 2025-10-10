// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <MueLu_ConfigDefs.hpp>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_Utilities.hpp>

#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

#if defined(HAVE_MUELU_ML)
#include <ml_MultiLevelPreconditioner.h>
#endif

#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

// Default problem is Laplace1D with nx = 8748. Use --help to list available options.

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool success = false;
  try {
    // MPI initialization using Teuchos
    RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
    // int MyPID = comm->getRank();
    // int NumProc = comm->getSize();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream &fancyout  = *fancy;
    fancyout.setOutputToRootOnly(0);

    // Convenient definitions
    using STS            = Teuchos::ScalarTraits<SC>;
    using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
    // const SC zero = Teuchos::ScalarTraits<SC>::zero();
    // const SC one = Teuchos::ScalarTraits<SC>::one();

    // Initialize and read parameters from command line
    std::string xmlFileName;
    clp.setOption("xml", &xmlFileName, "read parameters from a file. Otherwise, this example uses by default an hard-coded parameter list.");
    int num_iters = 9;
    clp.setOption("numIters", &num_iters, "Max number of iterations");

    Xpetra::Parameters xpetraParameters(clp);
    Galeri::Xpetra::Parameters<GO> matrixParameters(clp);

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    // Construct the problem
    RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector>> Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    RCP<Matrix> A = Pr->BuildMatrix();

    // Preconditioner configuration
    RCP<Teuchos::ParameterList> params;
    if (xmlFileName != "") {
      fancyout << "Reading " << xmlFileName << " ..." << std::endl;
      //! [ReadParametersFromXMLFile begin]
      params = Teuchos::getParametersFromXmlFile(xmlFileName);
      //! [ReadParametersFromXMLFile end]
    } else {
      fancyout << "Using hard-coded parameter list:" << std::endl;
      //! [ParameterList begin]
      params = rcp(new Teuchos::ParameterList());

      params->set("ML output", 10);
      params->set("max levels", 2);
      params->set("smoother: type", "symmetric Gauss-Seidel");
      params->set("coarse: type", "Amesos-KLU");
      //! [ParameterList end]
    }

    fancyout << "Initial parameter list" << std::endl;
    fancyout << *params << std::endl;

    // Construct a multigrid preconditioner

    // Build default null space
    LocalOrdinal numPDEs = 1;
    if (A->IsView("stridedMaps") == true) {
      Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps");
      numPDEs                     = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize();
      oldView                     = A->SwitchToView(oldView);
    }

    //! [BuildDefaultNullSpace begin]
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getDomainMap(), numPDEs);
    for (int i = 0; i < numPDEs; ++i) {
      Teuchos::ArrayRCP<Scalar> nsValues = nullspace->getDataNonConst(i);
      const int numBlocks                = nsValues.size() / numPDEs;

      for (int j = 0; j < numBlocks; ++j)
        nsValues[j * numPDEs + i] = STS::one();
    }
    //! [BuildDefaultNullSpace end]

    //! [MultigridHierarchy begin]
    MLParameterListInterpreter mueLuFactory(*params);
    RCP<Hierarchy> hierarchy = mueLuFactory.CreateHierarchy();
    //! [MultigridHierarchy end]

    //! [FeedInInformation begin]
    hierarchy->GetLevel(0)->Set("Nullspace", nullspace);
    hierarchy->GetLevel(0)->Set("A", A);
    //! [FeedInInformation end]

    // build hierarchy
    //! [CallSetupRoutine begin]
    mueLuFactory.SetupHierarchy(*hierarchy);
    //! [CallSetupRoutine end]

    // Setup vectors X and B to complete the linear system Ax = b
    RCP<Vector> x_vec = VectorFactory::Build(map);
    RCP<Vector> b_vec = VectorFactory::Build(map);

    x_vec->putScalar(STS::zero());
    b_vec->setSeed(846930886);
    b_vec->randomize();

    // Solve Ax = b with AMG as a standalone solver
    hierarchy->IsPreconditioner(false);
    hierarchy->Iterate(*b_vec, *x_vec, num_iters);

    // Print relative residual norm
    magnitude_type residualNorms = Utilities::ResidualNorm(*A, *x_vec, *b_vec)[0];
    fancyout << "||Residual|| = " << residualNorms << std::endl;

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main_

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
