// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Amesos2.hpp>
#include <Amesos2_config.h>

#include "Xpetra_ConfigDefs.hpp"
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_IO.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosStatusTestCombo.hpp>
#include <BelosXpetraStatusTestGenResSubNorm.hpp>

#include <Galeri_XpetraMaps.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ConfigDefs.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <Xpetra_Map.hpp>

#include <MatrixMarket_Tpetra.hpp>

//#{

#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Teuchos_FancyOStream.hpp>

void printMultiVector(const Teuchos::RCP<Tpetra::MultiVector<>> mv) {
  auto comm = mv->getMap()->getComm();
  Teuchos::RCP<Teuchos::FancyOStream> out =
      Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  // Prints in rank order and includes map info + entries at EXTREME
  mv->describe(*out, Teuchos::VERB_EXTREME);
}

/*Usage:
comm->barrier();
std::cout<<"\n//#exactSolution:\n\n";
comm->barrier();
printMultiVector(exactSolution);
comm->barrier();
*/

//#}

int main(int argc, char* argv[]) {
#include "MueLu_UseShortNames.hpp"
  bool success = false;
  try {
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    {
      //! [TpetraTemplateParameters begin]
      using SC = Tpetra::MultiVector<>::scalar_type;
      using LO = Tpetra::MultiVector<>::local_ordinal_type;
      using GO = Tpetra::MultiVector<>::global_ordinal_type;
      using NO = Tpetra::MultiVector<>::node_type;
      //! [TpetraTemplateParameters end]

      // ================================
      // Convenient definitions
      // ================================
      //! [UsingStatements begin]
      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;

      using CrsMatrix   = Tpetra::CrsMatrix<>;
      using Map         = Tpetra::Map<>;
      using MultiVector = Tpetra::MultiVector<>;
      using Operator    = Tpetra::Operator<>;

      using STS            = Teuchos::ScalarTraits<SC>;
      using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

      using real_type             = typename STS::coordinateType;
      using RealValuedMultiVector = Tpetra::MultiVector<real_type, LO, GO, NO>;

      const SC one = Teuchos::ScalarTraits<SC>::one();
      //! [UsingStatements end]

      //! [CommunicatorObject begin]
      RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
      int MyPID                          = comm->getRank();
      int NumProc                        = comm->getSize();
      (void)MyPID;    // unused void pointer cast to avoid unused variable warnings
      (void)NumProc;  // unused void pointer cast to avoid unused variable warnings
      //! [CommunicatorObject end]

      // Instead of checking each time for rank, create a rank 0 stream
      RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      Teuchos::FancyOStream& fancyout  = *fancy;
      fancyout.setOutputToRootOnly(0);

      // ================================
      // Parameters initialization
      // ================================
      Teuchos::CommandLineProcessor clp(false);

      std::string xmlFileName = "xml/muelu_ParameterList.xml";
      clp.setOption("xml", &xmlFileName, "read parameters from a file [default = 'xml/muelu_ParameterList.xml']");

      int globalNumDofs = 0;  // 7020;
      clp.setOption("globalNumDofs", &globalNumDofs, "global number of degrees of freedom [has to be set by user, default = 0 -> error]");
      int nDofsPerNode = 1;
      clp.setOption("nDofsPerNode", &nDofsPerNode, "number of degrees of freedom per node [has to be set by user, default = 1]");
      int nProcs             = comm->getSize();
      std::string dsolveType = "cg";
      clp.setOption("solver", &dsolveType, "solve type: (none | cg | gmres | standalone) [default = cg]");  //# none vs standalone, whats the difference?
      double dtol = 1e-12;
      clp.setOption("tol", &dtol, "solver convergence tolerance [default = 1e-12]");
      std::string problemFile = "stru2d";
      clp.setOption("problem", &problemFile, "prefix for problem files (e.g. 'stru2d' expects 'stru2d_A.txt', 'stru2d_b.txt', and 'stru2d_ns.txt' files to be provided)");
      std::string coordsFile = "";
      clp.setOption("coordinates", &coordsFile, "file name containing coordinates in matrix market format");

      switch (clp.parse(argc, argv)) {
        case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS; break;
        case Teuchos::CommandLineProcessor::PARSE_ERROR:
        case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
        case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
      }

      // ================================
      // Validation of input parameters
      // ================================
      TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName == "", std::runtime_error,
                                 "You need to specify the xml-file via the command line argument '--xml=<path/to/xml_file>'.");
      if (globalNumDofs == 0) {
        std::cout << "Please specify '--globalNumDofs'! Simulation cannot run without that parameter correctly set" << std::endl;
        return EXIT_FAILURE;
      }

      // ================================
      // Problem construction
      // ================================
      int nLocalDofs     = (int)globalNumDofs / nProcs;
      nLocalDofs         = nLocalDofs - (nLocalDofs % nDofsPerNode);
      int nCumulatedDofs = 0;
      MueLu_sumAll(comm, nLocalDofs, nCumulatedDofs);

      if (comm->getRank() == nProcs - 1) {
        nLocalDofs += globalNumDofs - nCumulatedDofs;
      }

      Map tmap(globalNumDofs, nLocalDofs, 0, comm);  //# possible error

      // read in problem
      fancyout << "========================================================\nReading matrix market files.\n========================================================" << std::endl;

      Teuchos::RCP<const CrsMatrix> tA    = Tpetra::MatrixMarket::Reader<CrsMatrix>::readSparseFile(problemFile + "_A.txt", comm, tmap, tmap, tmap);  //# possibly dont pass comm, and colmap is output I guess?
      Teuchos::RCP<const MultiVector> tB  = Tpetra::MatrixMarket::Reader<MultiVector>::readVectorFile(problemFile + "_B.txt", comm, tmap);
      Teuchos::RCP<const MultiVector> tNS = Tpetra::MatrixMarket::Reader<MultiVector>::readDenseFile(problemFile + "_NS.txt", comm, tmap);

      // read in coordinates
      RCP<RealValuedMultiVector> tCoords;
      if (coordsFile != "") {
        Map coordsMap(globalNumDofs / nDofsPerNode, nLocalDofs / nDofsPerNode, 0, comm);
        tCoords = Tpetra::MatrixMarket::Reader<RealValuedMultiVector>::readDenseFile(fileNS, comm, coordsMap);  //#
        //#?xCoords                          = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<int, Node>(epCoords));
      }

      ////////////////////

      //#?
      RCP<const Map> nodeMap = Xpetra::toTpetra(xpetra_map);
      RCP<const Map> dofMap  = nodeMap;
      //#?
      // Some safety checks to see, if Galeri delivered valid output
      TEUCHOS_ASSERT(!nodeMap.is_null());
      TEUCHOS_ASSERT(!dofMap.is_null());
      TEUCHOS_ASSERT(!coordinates.is_null());
      TEUCHOS_ASSERT(!matrix.is_null());

      fancyout << "========================================================\nGaleri complete.\n========================================================" << std::endl;  //#?

      // ================================
      // Preconditioner construction
      // ================================

      //! [ReadMueLuParamsFromXmlFile begin]
      // Read MueLu parameter list from xml file
      RCP<ParameterList> mueluParams = Teuchos::getParametersFromXmlFile(xmlFileName);
      //! [ReadMueLuParamsFromXmlFile end]

      //! [InsertNullspaceInUserData begin]
      // Register nullspace as user data in the MueLu parameter list
      ParameterList& userDataList = mueluParams->sublist("user data");
      userDataList.set<RCP<MultiVector>>("Nullspace", nullspace);
      //! [InsertNullspaceInUserData end]

      //! [CreateTpetraPreconditioner begin]
      // Create the MueLu preconditioner based on the Tpetra stack

      // Tpetra to Xpetra
      RCP<Xpetra::CrsMatrix</*SC,LO,GO,NO//#can we omit these as we do in the usings for Tpetra?*/>> xA = Xpetra::toXpetra(tA);
      RCP<Xpetra::Matrix<SC, LO, GO, NO>> xMat                                                          = rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA));  //# automatically rcp_dynamic_cast?
      xMat->->SetFixedBlockSize(nDofsPerNode);

      RCP<Xpetra::MultiVector<SC, LO, GO, NO>> xNS = rcp(new Xpetra::TpetraMultiVector<SC, LO, GO, NO>(tNS));

      const RCP<const Map> map = Xpetra::toXpetra<GO, Node>(tMap);

      ParameterListInterpreter mueLuFactory(xmlFileName, *comm);
      RCP<Hierarchy> H         = mueLuFactory.CreateHierarchy();
      RCP<MueLu::Level> Finest = H->GetLevel(0);
      Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
      Finest->Set("A", xMat);
      Finest->Set("Nullspace", xNS);
      if (xCoords != Teuchos::null) Finest->Set("Coordinates", xCoords);

      mueLuFactory.SetupHierarchy(*H);

      RCP<MueLu::TpetraOperator<SC, LO, GO, NO>> mueLuPreconditioner(H);
      //! [CreateTpetraPreconditioner end]
      // Solve Ax = b using AMG as a preconditioner in Belos
      //! [MueLuAsPrecCreateSolutionVector begin]
      //! [MueLuAsPrecCreateSolutionVector end]
      RCP<MultiVector> tX = rcp(new MultiVector(dofMap, 1, true));
      /*// Initialize solution vector with random values //# we dont actually do this in Challenge.cpp. Maybe add an optional cli arg randomInitialGuess=false by default?//#
      RCP<MultiVector> tX = rcp(new MultiVector(dofMap, 1, true));
      STS::seedrandom(100);
      tX->randomize();*/

      //! [MueLuAsPrecSetupLinearSystem begin]
      // Construct a Belos LinearProblem object and hand-in the MueLu preconditioner
      RCP<Belos::LinearProblem<SC, MultiVector, Operator>> belosProblem =
          rcp(new Belos::LinearProblem<SC, MultiVector, Operator>(tA, tX, tB));
      belosProblem->setLeftPrec(mueLuPreconditioner);
      bool set = belosProblem->setProblem();
      //! [MueLuAsPrecSetupLinearSystem end]

      if (set == false) {
        fancyout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return EXIT_FAILURE;
      }

      //! [MueLuAsPrecConfigureAndCreateBelosSolver begin]
      // Belos parameter list
      RCP<ParameterList> belosList = Teuchos::parameterList();
      belosList->set("Maximum Iterations", 500);      // Maximum number of iterations allowed
      belosList->set("Convergence Tolerance", dtol);  // Relative convergence tolerance requested
      belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList->set("Output Frequency", 1);
      belosList->set("Output Style", Belos::Brief);

      // Create an iterative solver manager
      Belos::SolverFactory<SC, MultiVector, Operator> solverFactory;
      RCP<Belos::SolverManager<SC, MultiVector, Operator>> solver;  // = solverFactory.create("Block GMRES", belosList);
      if (dsolveType == "cg") {
        solver = solverFactory.create("CG", belosList);
      } else if (dsolveType == "gmres") {
        solver = solverFactory.create("BLOCK GMRES", belosList);
      } else {
        solver = solverFactory.create("FIXED POINT", belosList);
      }

      // Pass the linear problem to the solver
      solver->setProblem(belosProblem);
      //! [MueLuAsPrecConfigureAndCreateBelosSolver end]

      //! [MueLuAsPrecSolve begin]
      // Perform solve
      Belos::ReturnType retStatus = Belos::Unconverged;
      retStatus                   = solver->solve();
      //! [MueLuAsPrecSolve end]

      // Get the number of iterations for this solve
      fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;  //# wasnt in Challenge.cpp, but we can keep it

      // Print relative residual norm
      RCP<Xpetra::MultiVector<SC, LO, GO, NO>> xX = rcp(new Xpetra::TpetraMultiVector<SC, LO, GO, NO>(tX));
      auto norms                                  = MueLu::Utilities<SC, LO, GO, NO>::ResidualNorm(*xA, *xX, *xB);  // per column
      if (comm->getRank() == 0) fancyout << "||Residual|| = " << norms[0] << "\n";

      // Check convergence status
      if (retStatus != Belos::Converged)
        fancyout << std::endl
                 << "ERROR:  Belos did not converge! " << std::endl;
      else
        fancyout << std::endl
                 << "SUCCESS:  Belos converged!" << std::endl;
    }  // end of Tpetra::ScopeGuard
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main
