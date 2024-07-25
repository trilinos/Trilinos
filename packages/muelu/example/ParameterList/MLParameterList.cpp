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

#include <Tpetra_KokkosCompat_DefaultNode.hpp>  // For Epetra only runs this points to FakeKokkos in Xpetra

#include <Teuchos_XMLParameterListHelpers.hpp>  // getParametersFromXmlFile()

#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
#include <Epetra_CrsMatrix.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif

#ifdef HAVE_MUELU_AZTECOO
#include <AztecOO.h>
#endif

#ifdef HAVE_MUELU_EPETRA
#include <MueLu_EpetraOperator.hpp>
#endif

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

#include <MueLu_ML2MueLuParameterTranslator.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Default problem is Laplace1D with nx = 8748. Use --help to list available options.

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;

  //
  // MPI initialization using Teuchos
  //

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    //
    // Parameters
    //

    // TODO: FIXME: option by default does not work for MueLu/Tpetra

    int nIts = 9;

    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 256);  // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                   // manage parameters of xpetra

    std::string xmlFileName;
    clp.setOption("xml", &xmlFileName, "read parameters from a file. Otherwise, this example uses by default an hard-coded parameter list.");
    int translatedmuelu = true;
    clp.setOption("muelu2", &translatedmuelu, "use muelu through XML parameter translation and ParameterListInterpreter");
    int ml = true;
#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA)
    clp.setOption("ml", &ml, "use ml");
#endif

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    if (comm->getRank() == 0) {
      std::cout << xpetraParameters << matrixParameters;
    }

    // choose ML and Tpetra
    if (ml && xpetraParameters.GetLib() == Xpetra::UseTpetra) {
      ml = false;
      std::cout << "ML preconditionner can only be built if --linAlgebra=Epetra. Option --ml ignored" << std::endl;
    }

    //
    // Construct the problem
    //

    RCP<const Map> map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);
    RCP<Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    RCP<Matrix> A = Pr->BuildMatrix();

    //
    // Preconditioner configuration
    //

    // ML parameter list
    RCP<Teuchos::ParameterList> params;
    if (xmlFileName != "") {
      std::cout << "Reading " << xmlFileName << " ..." << std::endl;
      params = Teuchos::getParametersFromXmlFile(xmlFileName);

    } else {
      std::cout << "Using hard-coded parameter list:" << std::endl;
      params = rcp(new Teuchos::ParameterList());

      params->set("ML output", 10);
      params->set("max levels", 2);
      params->set("smoother: type", "symmetric Gauss-Seidel");

      if (xpetraParameters.GetLib() == Xpetra::UseTpetra)  // TODO: remove 'if' when Amesos2-KLU becomes available
        params->set("coarse: type", "Amesos-Superlu");
      else
        params->set("coarse: type", "Amesos-KLU");
    }

    std::cout << "Initial parameter list" << std::endl;
    std::cout << *params << std::endl;

    if (translatedmuelu) {
      //
      // Construct a multigrid preconditioner
      //

      RCP<Teuchos::ParameterList> mueluParams = Teuchos::getParametersFromXmlString(MueLu::ML2MueLuParameterTranslator::translate(*params, "SA"));

      std::cout << MueLu::ML2MueLuParameterTranslator::translate(*params, "SA") << std::endl;

      if (xpetraParameters.GetLib() == Xpetra::UseEpetra)
        mueluParams->set("use kokkos refactor", false);

      // Multigrid Hierarchy
      ParameterListInterpreter mueLuFactory(*mueluParams);
      RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();

      // build default null space
      LocalOrdinal numPDEs = 1;
      if (A->IsView("stridedMaps") == true) {
        Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps");  // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
        numPDEs                     = Teuchos::rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize();
        oldView                     = A->SwitchToView(oldView);
      }

      RCP<MultiVector> nullspace = MultiVectorFactory::Build(A->getDomainMap(), numPDEs);

      for (int i = 0; i < numPDEs; ++i) {
        Teuchos::ArrayRCP<Scalar> nsValues = nullspace->getDataNonConst(i);
        int numBlocks                      = nsValues.size() / numPDEs;
        for (int j = 0; j < numBlocks; ++j) {
          nsValues[j * numPDEs + i] = 1.0;
        }
      }

      H->GetLevel(0)->Set("Nullspace", nullspace);
      H->GetLevel(0)->Set("A", A);

      //
      // build hierarchy
      //
      mueLuFactory.SetupHierarchy(*H);

      //
      // Solve Ax = b
      //

      RCP<Vector> X = VectorFactory::Build(map);
      RCP<Vector> B = VectorFactory::Build(map);

      X->putScalar((Scalar)0.0);
      B->setSeed(846930886);
      B->randomize();

      // AMG as a standalone solver
      H->IsPreconditioner(false);
      H->Iterate(*B, *X, nIts);

      // Print relative residual norm
      typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = Utilities::ResidualNorm(*A, *X, *B)[0];
      if (comm->getRank() == 0)
        std::cout << "||Residual|| = " << residualNorms << std::endl;

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AZTECOO) && !defined(HAVE_MUELU_CUDA) && !defined(HAVE_MUELU_HIP) && !defined(HAVE_MUELU_SYCL)
      // TODO TAW: 4/8/2016
      // temporarely deactivate this due to runtime error on perseus:
      // Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed
      // if SERIAL=OFF, OPENMP=OFF, PTHREAD=ON, CUDA=OFF
      // probably a fix necessary in EpetraOperator (which only supports
      // SERIAL or OPENMP, but not PTHREAD of course).

      if (xpetraParameters.GetLib() == Xpetra::UseEpetra) {  // TODO: should be doable with Tpetra too

        // AMG as a preconditioner

        // TODO: name mueluPrec and mlPrec not

        H->IsPreconditioner(true);
        MueLu::EpetraOperator mueluPrec(H);  // Wrap MueLu preconditioner into an Epetra Operator

        //
        // Solve Ax = b
        //
        RCP<Epetra_CrsMatrix> eA;  // duplicate code
        {                          // TODO: simplify this
          RCP<CrsMatrixWrap> xCrsOp    = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true);
          RCP<CrsMatrix> xCrsMtx       = xCrsOp->getCrsMatrix();
          RCP<EpetraCrsMatrix> eCrsMtx = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xCrsMtx, true);
          eA                           = eCrsMtx->getEpetra_CrsMatrixNonConst();
        }

        RCP<Epetra_Vector> eX = rcp(new Epetra_Vector(eA->RowMap()));
        RCP<Epetra_Vector> eB = rcp(new Epetra_Vector(eA->RowMap()));

        eX->PutScalar((Scalar)0.0);
        eB->SetSeed(846930886);
        eB->Random();

        Epetra_LinearProblem eProblem(eA.get(), eX.get(), eB.get());

        // AMG as a standalone solver
        AztecOO solver(eProblem);
        solver.SetPrecOperator(&mueluPrec);
        solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
        solver.SetAztecOption(AZ_output, 1);

        solver.Iterate(nIts, 1e-10);

        {  // TODO: simplify this
          RCP<Vector> mueluX = rcp(new Xpetra::EpetraVectorT<int, Node>(eX));
          RCP<Vector> mueluB = rcp(new Xpetra::EpetraVectorT<int, Node>(eB));
          // Print relative residual norm
          typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms2 = Utilities::ResidualNorm(*A, *mueluX, *mueluB)[0];
          if (comm->getRank() == 0)
            std::cout << "||Residual|| = " << residualNorms2 << std::endl;
        }

        // TODO: AMG as a preconditioner (AZ_cg)
      }
#endif  // HAVE_MUELU_AZTECOO
    }   // if (translatedmuelu)

#if defined(HAVE_MUELU_ML) && defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AZTECOO) && !defined(HAVE_MUELU_CUDA) && !defined(HAVE_MUELU_HIP) && !defined(HAVE_MUELU_SYCL)
    if (ml) {
      std::cout << std::endl
                << std::endl
                << std::endl
                << std::endl
                << "**** ML ml ML ml ML" << std::endl
                << std::endl
                << std::endl
                << std::endl;

      //
      // Construct a multigrid preconditioner
      //

      // Multigrid Hierarchy
      RCP<CrsMatrixWrap> crsOp                  = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true);
      RCP<CrsMatrix> crsMtx                     = crsOp->getCrsMatrix();
      RCP<EpetraCrsMatrix> epetraCrsMtx         = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(crsMtx, true);
      RCP<const Epetra_CrsMatrix> epetra_CrsMtx = epetraCrsMtx->getEpetra_CrsMatrix();

      RCP<Epetra_CrsMatrix> eA;
      {  // TODO: simplify this
        RCP<CrsMatrixWrap> xCrsOp    = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(A, true);
        RCP<CrsMatrix> xCrsMtx       = xCrsOp->getCrsMatrix();
        RCP<EpetraCrsMatrix> eCrsMtx = Teuchos::rcp_dynamic_cast<EpetraCrsMatrix>(xCrsMtx, true);
        eA                           = eCrsMtx->getEpetra_CrsMatrixNonConst();
      }

      RCP<ML_Epetra::MultiLevelPreconditioner> mlPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*eA, *params));

      //
      // Solve Ax = b
      //

      RCP<Epetra_Vector> eX = rcp(new Epetra_Vector(eA->RowMap()));
      RCP<Epetra_Vector> eB = rcp(new Epetra_Vector(eA->RowMap()));

      eX->PutScalar((Scalar)0.0);
      eB->SetSeed(846930886);
      eB->Random();

      Epetra_LinearProblem eProblem(eA.get(), eX.get(), eB.get());

      // AMG as a standalone solver
      AztecOO solver(eProblem);
      solver.SetPrecOperator(mlPrec.get());
      solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
      solver.SetAztecOption(AZ_output, 1);

      solver.Iterate(nIts, 1e-10);

      {  // TODO: simplify this
        RCP<Vector> mueluX = rcp(new Xpetra::EpetraVectorT<int, Node>(eX));
        RCP<Vector> mueluB = rcp(new Xpetra::EpetraVectorT<int, Node>(eB));
        // Print relative residual norm
        typename Teuchos::ScalarTraits<SC>::magnitudeType residualNorms = Utilities::ResidualNorm(*A, *mueluX, *mueluB)[0];
        if (comm->getRank() == 0)
          std::cout << "||Residual|| = " << residualNorms << std::endl;
      }

      std::cout << "Parameter list after ML run" << std::endl;
      const Teuchos::ParameterList &paramsAfterML = mlPrec->GetList();
      std::cout << paramsAfterML << std::endl;

    }  // if (ml)

#endif  // HAVE_MUELU_ML && HAVE_MUELU_EPETRA

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
  bool success = false;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize(argc, argv);
  try {
    const bool throwExceptions     = false;
    const bool recogniseAllOptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions, recogniseAllOptions);
    Xpetra::Parameters xpetraParameters(clp);

    std::string node = "";
    clp.setOption("node", &node, "node type (serial | openmp | cuda | hip)");

    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double, int, int, Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }
    if (lib == Xpetra::UseTpetra) {
      std::cout << "Skipped tests for Tpetra. We officially only support the MLParameterListInterpreter for Epetra. It is supposed to be a transition from ML with Epetra to MueLu. Furthermore, there is only support for Epetra and not for Epetra64. That is, only GO=int allowed." << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  Kokkos::finalize();

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}
