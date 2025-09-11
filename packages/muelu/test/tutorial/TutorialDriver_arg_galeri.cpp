//# follow TutorialDriver_xml.cpp closely, but with --args rather than an xml

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

#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <Xpetra_Map.hpp>

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
      using Teuchos::TimeMonitor;

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

      std::string availableMatrixTypes[] = {"Laplace2D", "Laplace3D", "Recirc2D"};
      std::string matrixType             = availableMatrixTypes[0];
      std::string matrixTypeDesc         = "matrix type which defines problem. Currently available options: [";
      for (size_t i = 0; i < std::size(availableMatrixTypes); ++i) {
        matrixTypeDesc += availableMatrixTypes[i];
        if (i < std::size(availableMatrixTypes) - 1) matrixTypeDesc += ", ";
      }
      matrixTypeDesc += "]";
      clp.setOption("matrixType", &matrixType, matrixTypeDesc.c_str());

      GO nx = 100;
      clp.setOption("nx", &nx, "mesh size in x direction");
      GO ny = 100;
      clp.setOption("ny", &ny, "mesh size in y direction");
      GO nz = 100;
      clp.setOption("nz", &ny, "mesh size in z direction");

      std::string xmlFileName = "";
      clp.setOption("xml", &xmlFileName, "read parameters from a file");
      int mgridSweeps = 1;
      clp.setOption("mgridSweeps", &mgridSweeps, "number of multigrid sweeps within multigrid solver");
      std::string printTimings = "no";
      clp.setOption("timings", &printTimings, "print timings to screen [yes/no]");
      double tol = 1e-12;
      clp.setOption("tol", &tol, "solver convergence tolerance");

      double diffusion = 1e-5;
      clp.setOption("diffusion", &diffusion, "diffusion coefficient, also called epsilon");

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
      TEUCHOS_TEST_FOR_EXCEPTION([&] {
        for (std::string e : availableMatrixTypes)
          if (matrixType == e) return false;
        return true;
      }(),
                                 std::runtime_error, "Invalid matrixType.");

      RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time"))), tm;

      // ================================
      // Problem construction
      // ================================
      comm->barrier();
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

      Teuchos::ParameterList galeriList;
      galeriList.set("nx", nx);
      galeriList.set("ny", ny);
      galeriList.set("mx", comm->getSize());
      galeriList.set("my", 1);

      if (matrixType == "Laplace3D") {
        galeriList.set("nz", nz);
        galeriList.set("mz", 1);
      }

      if (matrixType == "Recirc2D") {
        galeriList.set("lx", 1.0);  // length of x-axis
        galeriList.set("ly", 1.0);  // length of y-axis
        galeriList.set("diff", diffusion);
        galeriList.set("conv", 1.0);
      }

      // Create node map (equals dof map, since one dof per node)
      RCP<const Xpetra::Map<LO, GO, NO>> xpetra_map;
      if (matrixType == "Laplace2D" || matrixType == "Recirc2D")
        xpetra_map = Galeri::Xpetra::CreateMap<LO, GO, NO>(Xpetra::UseTpetra, "Cartesian2D", comm, galeriList);
      else if (matrixType == "Laplace3D")
        xpetra_map = Galeri::Xpetra::CreateMap<LO, GO, NO>(Xpetra::UseTpetra, "Cartesian3D", comm, galeriList);
      RCP<const Map> nodeMap = Xpetra::toTpetra(xpetra_map);
      RCP<const Map> dofMap  = nodeMap;

      // Create coordinates
      RCP<const RealValuedMultiVector> coordinates;
      if (matrixType == "Laplace2D" || matrixType == "Recirc2D")
        coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("2D", nodeMap, galeriList);
      else if (matrixType == "Laplace3D")
        coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC, LO, GO, Map, RealValuedMultiVector>("3D", nodeMap, galeriList);

      // Create the matrix
      RCP<Galeri::Xpetra::Problem<Map, CrsMatrix, MultiVector>> galeriProblem =
          Galeri::Xpetra::BuildProblem<SC, LO, GO, Map, CrsMatrix, MultiVector>(matrixType, dofMap, galeriList);
      RCP<CrsMatrix> matrix = galeriProblem->BuildMatrix();

      // Some safety checks to see, if Galeri delivered valid output
      TEUCHOS_ASSERT(!nodeMap.is_null());
      TEUCHOS_ASSERT(!dofMap.is_null());
      TEUCHOS_ASSERT(!coordinates.is_null());
      TEUCHOS_ASSERT(!matrix.is_null());

      //! [RhsAndSolutionVector begin]
      // Create right-hand side (with all ones)
      RCP<MultiVector> B = rcp(new MultiVector(dofMap, 1, true));
      B->putScalar(one);

      // Initialize solution vector with random values
      RCP<MultiVector> X0 = rcp(new MultiVector(dofMap, 1, true));
      STS::seedrandom(100);
      X0->randomize();
      //! [RhsAndSolutionVector end]

      //! [BuildNullSpaceVector begin]
      // Build null space vector for a scalar problem, i.e. vector with all ones
      RCP<MultiVector> nullspace = rcp(new MultiVector(dofMap, 1));
      nullspace->putScalar(one);
      //! [BuildNullSpaceVector end]

      comm->barrier();
      tm = Teuchos::null;

      fancyout << "========================================================\nGaleri complete.\n========================================================" << std::endl;

      // ================================
      // Preconditioner construction
      // ================================

      comm->barrier();
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1.5 - MueLu read XML")));

      //! [ReadMueLuParamsFromXmlFile begin]
      // Read MueLu parameter list from xml file
      RCP<ParameterList> mueluParams = Teuchos::getParametersFromXmlFile(xmlFileName);
      //! [ReadMueLuParamsFromXmlFile end]

      comm->barrier();
      tm = Teuchos::null;

      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

      //! [InsertNullspaceInUserData begin]
      // Register nullspace as user data in the MueLu parameter list
      ParameterList& userDataList = mueluParams->sublist("user data");
      userDataList.set<RCP<MultiVector>>("Nullspace", nullspace);
      //! [InsertNullspaceInUserData end]

      //! [CreateTpetraPreconditioner begin]
      // Create the MueLu preconditioner based on the Tpetra stack
      RCP<MueLu::TpetraOperator<SC, LO, GO, NO>> mueLuPreconditioner =
          MueLu::CreateTpetraPreconditioner(Teuchos::rcp_dynamic_cast<Operator>(matrix), *mueluParams);
      //! [CreateTpetraPreconditioner end]

      comm->barrier();
      tm = Teuchos::null;

      // Generate exact solution using a direct solver
      //! [ExactSolutionVector begin]
      RCP<MultiVector> exactSolution = rcp(new MultiVector(dofMap, 1, true));
      //! [ExactSolutionVector end]
      {
        fancyout << "========================================================\nCalculate exact solution." << std::endl;
        tm                                                        = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - direct solve")));
        RCP<Amesos2::Solver<CrsMatrix, MultiVector>> directSolver = Amesos2::create<CrsMatrix, MultiVector>("KLU2", matrix, exactSolution, B);
        directSolver->solve();

        comm->barrier();
        std::cout << "\n//#exactSolution:\n\n";
        comm->barrier();
        printMultiVector(exactSolution);
        comm->barrier();

        comm->barrier();
        tm = Teuchos::null;
      }

      // Solve Ax = b using AMG as a preconditioner in Belos
      //! [MueLuAsPrecCreateSolutionVector begin]
      RCP<MultiVector> precSolVec = rcp(new MultiVector(dofMap, 1, true));
      //! [MueLuAsPrecCreateSolutionVector end]
      {
        fancyout << "========================================================\nUse multigrid hierarchy as preconditioner within CG." << std::endl;
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - AMG as preconditioner")));

        precSolVec->update(1.0, *X0, 0.0);

        //! [MueLuAsPrecSetupLinearSystem begin]
        // Construct a Belos LinearProblem object and hand-in the MueLu preconditioner
        RCP<Belos::LinearProblem<SC, MultiVector, Operator>> belosProblem =
            rcp(new Belos::LinearProblem<SC, MultiVector, Operator>(matrix, precSolVec, B));
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
        belosList->set("Maximum Iterations", 50);      // Maximum number of iterations allowed
        belosList->set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
        belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
        belosList->set("Output Frequency", 1);
        belosList->set("Output Style", Belos::Brief);

        // Create an iterative solver manager
        Belos::SolverFactory<SC, MultiVector, Operator> solverFactory;
        RCP<Belos::SolverManager<SC, MultiVector, Operator>> solver = solverFactory.create("Block GMRES", belosList);

        // Pass the linear problem to the solver
        solver->setProblem(belosProblem);
        //! [MueLuAsPrecConfigureAndCreateBelosSolver end]

        //! [MueLuAsPrecSolve begin]
        // Perform solve
        Belos::ReturnType retStatus = Belos::Unconverged;
        retStatus                   = solver->solve();
        //! [MueLuAsPrecSolve end]

        comm->barrier();
        tm = Teuchos::null;

        // Get the number of iterations for this solve
        fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

        comm->barrier();
        std::cout << "\n//#precSolVec:\n\n";
        comm->barrier();
        printMultiVector(precSolVec);
        comm->barrier();

        // Check convergence status
        if (retStatus != Belos::Converged)
          fancyout << std::endl
                   << "ERROR:  Belos did not converge! " << std::endl;
        else
          fancyout << std::endl
                   << "SUCCESS:  Belos converged!" << std::endl;
      }

      // Solve Ax = b using AMG as a solver
      //! [MueLuAsSolverCreateSolutionVector begin]
      RCP<MultiVector> multigridSolVec = rcp(new MultiVector(dofMap, 1, true));
      //! [MueLuAsSolverCreateSolutionVector end]
      {
        fancyout << "========================================================\nUse multigrid hierarchy as solver." << std::endl;
        tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 5 - Multigrid Solve")));

        multigridSolVec->update(1.0, *X0, 0.0);

        //! [ExtractHierarchyFromTpetraPrec begin]
        // Extract the underlying MueLu hierarchy
        RCP<MueLu::Hierarchy<SC, LO, GO, NO>> hierarchy = mueLuPreconditioner->GetHierarchy();
        //! [ExtractHierarchyFromTpetraPrec end]

        //! [MueLuAsSolverSetSolverMode begin]
        // Configure MueLu to be used as solver
        hierarchy->IsPreconditioner(false);
        //! [MueLuAsSolverSetSolverMode end]

        //! [MueLuAsSolverIterate begin]
        // Solve
        hierarchy->Iterate(*Xpetra::toXpetra(B), *Xpetra::toXpetra(multigridSolVec), mgridSweeps);
        //! [MueLuAsSolverIterate end]

        comm->barrier();
        std::cout << "\n//#multigridSolVec:\n\n";
        comm->barrier();
        printMultiVector(multigridSolVec);
        comm->barrier();

        comm->barrier();
        tm = Teuchos::null;
      }

      // Write results into file
      /*{
        fancyout << "========================================================\nExport results.\n========================================================" << std::endl;
        std::ofstream myfile;
        std::stringstream ss;
        ss << "example" << MyPID << ".txt";
        myfile.open(ss.str().c_str());

        // loop over all procs
        for (int iproc = 0; iproc < NumProc; iproc++) {
          if (MyPID == iproc) {
            int NumVectors1               = 2;
            int NumMyElements1            = coordinates->getMap()->getLocalNumElements();
            int MaxElementSize1           = 1;
            int* FirstPointInElementList1 = NULL;
            if (MaxElementSize1 != 1) FirstPointInElementList1 = coordinates->getMap().FirstPointInElementList();
            double** A_Pointers = coordinates->Pointers();

            if (MyPID == 0) {
              myfile.width(8);
              myfile << "#     MyPID";
              myfile << "    ";
              myfile.width(12);
              if (MaxElementSize1 == 1)
                myfile << "GID  ";
              else
                myfile << "     GID/Point";
              for (int j = 0; j < NumVectors1; j++) {
                myfile.width(20);
                myfile << "Value  ";
              }
              myfile << std::endl;
            }
            for (int i = 0; i < NumMyElements1; i++) {
              for (int ii = 0; ii < coordinates->Map().ElementSize(i); ii++) {
                int iii;
                myfile.width(10);
                myfile << MyPID;
                myfile << "    ";
                myfile.width(10);
                if (MaxElementSize1 == 1) {
                  if (coordinates->Map().GlobalIndicesInt()) {
                    int* MyGlobalElements1 = coordinates->Map().MyGlobalElements();
                    myfile << MyGlobalElements1[i] << "    ";
                  }

                  iii = i;
                } else {
                  if (coordinates->Map().GlobalIndicesInt()) {
                    int* MyGlobalElements1 = coordinates->Map().MyGlobalElements();
                    myfile << MyGlobalElements1[i] << "/" << ii << "    ";
                  }

                  iii = FirstPointInElementList1[i] + ii;
                }
                for (int j = 0; j < NumVectors1; j++) {
                  myfile.width(20);
                  myfile << A_Pointers[j][iii];
                }

                myfile.precision(18);  // set high precision for output

                // add exact solution vector entry
                myfile.width(25);
                myfile << exactSolution->getData(iii).get();

                // add preconditioned solution vector entry
                myfile.width(25);
                myfile << precSolVec->getData(iii).get();

                // add multigrid solution vector entry
                myfile.width(25);
                myfile << multigridSolVec->getData(iii).get();

                myfile.precision(6);  // set default precision
                myfile << std::endl;
              }
            }  // end loop over all lines on current proc
            myfile << std::flush;

            // syncronize procs // Why three times?
            comm->barrier();
            comm->barrier();
            comm->barrier();

          }  // end myProc
        } // end loop over all procs

        myfile.close();

        comm->barrier();
        tm                = Teuchos::null;
        globalTimeMonitor = Teuchos::null;

        if (printTimings == "yes") {
          TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
        }
      }*/

    }  // end of Tpetra::ScopeGuard
    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
}  // main
