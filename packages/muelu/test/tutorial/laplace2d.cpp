// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
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

int main (int argc, char *argv[])
{
  bool success = false;
  try
  {
    Tpetra::ScopeGuard tpetraScope (&argc, &argv);
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

      using CrsMatrix = Tpetra::CrsMatrix<>;
      using Map = Tpetra::Map<>;
      using MultiVector = Tpetra::MultiVector<>;
      using Operator = Tpetra::Operator<>;

      using STS = Teuchos::ScalarTraits<SC>;
      using magnitude_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;

      using real_type = typename STS::coordinateType;
      using RealValuedMultiVector = Tpetra::MultiVector<real_type,LO,GO,NO>;

      const SC one = Teuchos::ScalarTraits<SC>::one();
      //! [UsingStatements end]

      //! [CommunicatorObject begin]
      RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
      int MyPID = comm->getRank();
      int NumProc = comm->getSize();
      (void) MyPID;   // unused void pointer cast to avoid unused variable warnings
      (void) NumProc; // unused void pointer cast to avoid unused variable warnings
      //! [CommunicatorObject end]

      // Instead of checking each time for rank, create a rank 0 stream
      RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      Teuchos::FancyOStream& fancyout = *fancy;
      fancyout.setOutputToRootOnly(0);

      // ================================
      // Parameters initialization
      // ================================
      Teuchos::CommandLineProcessor clp(false);
      GO nx                    = 100;   clp.setOption("nx",                       &nx, "mesh size in x direction");
      GO ny                    = 100;   clp.setOption("ny",                       &ny, "mesh size in y direction");
      std::string xmlFileName  = "";    clp.setOption("xml",             &xmlFileName, "read parameters from a file");
      int mgridSweeps          = 1;     clp.setOption("mgridSweeps",     &mgridSweeps, "number of multigrid sweeps within Multigrid solver.");
      std::string printTimings = "no";  clp.setOption("timings",        &printTimings, "print timings to screen [yes/no]");
      double tol               = 1e-12; clp.setOption("tol",                     &tol, "solver convergence tolerance");

      switch (clp.parse(argc,argv)) {
        case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
        case Teuchos::CommandLineProcessor::PARSE_ERROR:
        case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
        case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
      }

      // ================================
      // Validation of input parameters
      // ================================
      TEUCHOS_TEST_FOR_EXCEPTION(xmlFileName=="", std::runtime_error,
          "You need to specify the xml-file via the command line argument '--xml=<path/to/xml_file>'.");

      // ================================
      // Problem construction
      // ================================

      Teuchos::ParameterList galeriList;
      galeriList.set("nx", nx);
      galeriList.set("ny", ny);
      galeriList.set("mx", comm->getSize());
      galeriList.set("my", 1);

      // Create node map (equals dof map, since one dof per node)
      RCP<const Xpetra::Map<LO, GO, NO>> xpetra_map = Galeri::Xpetra::CreateMap<LO,GO,NO>(Xpetra::UseTpetra, "Cartesian2D", comm, galeriList);
      RCP<const Map> nodeMap = Xpetra::toTpetra(xpetra_map);
      RCP<const Map> dofMap = nodeMap;

      // Create coordinates
      RCP<const RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,RealValuedMultiVector>("2D", nodeMap, galeriList);

      // Create the matrix
      RCP<Galeri::Xpetra::Problem<Map,CrsMatrix,MultiVector>> galeriProblem =
          Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrix,MultiVector>("Laplace2D", dofMap, galeriList);
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

      // Initilize solution vector with random values
      RCP<MultiVector> X = rcp(new MultiVector(dofMap, 1, true));
      STS::seedrandom(100);
      X->randomize();
      //! [RhsAndSolutionVector end]

      //! [BuildNullSpaceVector begin]
      // Build null space vector for a scalar problem, i.e. vector with all ones
      RCP<MultiVector> nullspace = rcp(new MultiVector(dofMap, 1));
      nullspace->putScalar(one);
      //! [BuildNullSpaceVector end]

      // ================================
      // Preconditioner construction
      // ================================

      //! [ReadMueLuParamsFromXmlFile begin]
      // Read MueLu parameter list from xml file
      RCP<ParameterList> mueluParams =  Teuchos::getParametersFromXmlFile(xmlFileName);
      //! [ReadMueLuParamsFromXmlFile end]

      //! [InsertNullspaceInUserData begin]
      // Register nullspace as user data in the MueLu parameter list
      ParameterList& userDataList = mueluParams->sublist("user data");
      userDataList.set<RCP<MultiVector>>("Nullspace", nullspace);
      //! [InsertNullspaceInUserData end]

      //! [CreateTpetraPreconditioner begin]
      // Create the MueLu preconditioner based on the Tpetra stack
      RCP<MueLu::TpetraOperator<SC,LO,GO,NO>> mueLuPreconditioner =
          MueLu::CreateTpetraPreconditioner(Teuchos::rcp_dynamic_cast<Operator>(matrix), *mueluParams);
      //! [CreateTpetraPreconditioner end]

      // Generate exact solution using a direct solver
      //! [ExactSolutionVector begin]
      RCP<MultiVector> exactSolution = rcp(new MultiVector(dofMap, 1, true));
      //! [ExactSolutionVector end]
      {
        exactSolution->update(0.0, *X, 1.0);

        RCP<Amesos2::Solver<CrsMatrix,MultiVector>> directSolver = Amesos2::create<CrsMatrix,MultiVector>("KLU2", matrix, X, B);
        directSolver->solve();
      }

      // Solve Ax = b using AMG as a preconditioner in Belos
      {
        //! [MueLuAsPrecCreateSolutionVector begin]
        // Create solution vector and set initial guess (already stored in X)
        RCP<MultiVector> precSolVec = rcp(new MultiVector(dofMap, 1, true));
        precSolVec->update(0.0, *X, 1.0);
        //! [MueLuAsPrecCreateSolutionVector end]

        //! [MueLuAsPrecSetupLinearSystem begin]
        // Construct a Belos LinearProblem object and hand-in the MueLu preconditioner
        RCP<Belos::LinearProblem<SC,MultiVector,Operator>> belosProblem =
            rcp(new Belos::LinearProblem<SC,MultiVector,Operator>(matrix, precSolVec, B));
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
        belosList->set("Maximum Iterations", 50); // Maximum number of iterations allowed
        belosList->set("Convergence Tolerance", tol); // Relative convergence tolerance requested
        belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
        belosList->set("Output Frequency", 1);
        belosList->set("Output Style", Belos::Brief);

        // Create an iterative solver manager
        Belos::SolverFactory<SC,MultiVector,Operator> solverFactory;
        RCP<Belos::SolverManager<SC,MultiVector,Operator>> solver = solverFactory.create("Block GMRES", belosList);

        // Pass the linear problem to the solver
        solver->setProblem(belosProblem);
        //! [MueLuAsPrecConfigureAndCreateBelosSolver end]

        //! [MueLuAsPrecSolve begin]
        // Perform solve
        Belos::ReturnType retStatus = Belos::Unconverged;
        retStatus = solver->solve();
        //! [MueLuAsPrecSolve end]

        // Get the number of iterations for this solve
        fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

        // Check convergence status
        if (retStatus != Belos::Converged)
          fancyout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
        else
          fancyout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
      }

      // Solve Ax = b using AMG as a solver
      {
        //! [MueLuAsSolverCreateSolutionVector begin]
        // Create solution vector and set initial guess (already stored in X)
        RCP<MultiVector> multigridSolVec = rcp(new MultiVector(dofMap, 1, true));
        multigridSolVec->update(0.0, *X, 1.0);
        //! [MueLuAsSolverCreateSolutionVector end]

        //! [ExtractHierarchyFromTpetraPrec begin]
        // Extract the underlying MueLu hierarchy
        RCP<MueLu::Hierarchy<SC,LO,GO,NO>> hierarchy = mueLuPreconditioner->GetHierarchy();
        //! [ExtractHierarchyFromTpetraPrec end]

        //! [MueLuAsSolverSetSolverMode begin]
        // Configure MueLu to be used as solver
        hierarchy->IsPreconditioner(false);
        //! [MueLuAsSolverSetSolverMode end]

        //! [MueLuAsSolverIterate begin]
        // Solve
        hierarchy->Iterate(*Xpetra::toXpetra(B), *Xpetra::toXpetra(multigridSolVec), mgridSweeps);
        //! [MueLuAsSolverIterate end]
      }

    } // end of Tpetra::ScopeGuard
  success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // main
