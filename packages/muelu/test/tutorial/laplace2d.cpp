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
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

#include <Amesos2.hpp>
#include <Amesos2_config.h>

#include <BelosBlockCGSolMgr.hpp>
#include <BelosMueLuAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosStatusTestCombo.hpp>
#include <BelosXpetraAdapter.hpp>
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
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

int main (int argc, char *argv[])
{
  bool success = false;
  try
  {
    Tpetra::ScopeGuard tpetraScope (&argc, &argv);
    {
      // ================================
      // Convenient definitions
      // ================================
      using scalar_type = Tpetra::MultiVector<>::scalar_type;
      using local_ordinal_type = Tpetra::MultiVector<>::local_ordinal_type;
      using global_ordinal_type = Tpetra::MultiVector<>::global_ordinal_type;
      using node_type = Tpetra::MultiVector<>::node_type;

      using Scalar = scalar_type;
      using LocalOrdinal = local_ordinal_type;
      using GlobalOrdinal = global_ordinal_type;
      using Node = node_type;

      using Teuchos::ParameterList;
      using Teuchos::RCP;
      using Teuchos::rcp;
      using Teuchos::TimeMonitor;

      using crs_matrix_type = Tpetra::CrsMatrix<>;
      using map_type = Tpetra::Map<>;
      using multi_vector_type = Tpetra::MultiVector<>;
      using operator_type = Tpetra::Operator<>;

      using STS = Teuchos::ScalarTraits<scalar_type>;
      using magnitude_type = typename Teuchos::ScalarTraits<scalar_type>::magnitudeType;

      using real_type = typename STS::coordinateType;
      using RealValuedMultiVector = Tpetra::MultiVector<real_type,local_ordinal_type,global_ordinal_type,node_type>;

      const scalar_type zero = Teuchos::ScalarTraits<scalar_type>::zero();
      const scalar_type one = Teuchos::ScalarTraits<scalar_type>::one();

      //! [CommunicatorObject begin]
      RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
      int MyPID = comm->getRank();
      int NumProc = comm->getSize();
      //! [CommunicatorObject end]

      // Instead of checking each time for rank, create a rank 0 stream
      RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      Teuchos::FancyOStream& fancyout = *fancy;
      fancyout.setOutputToRootOnly(0);

      // ================================
      // Parameters initialization
      // ================================
      Teuchos::CommandLineProcessor clp(false);
      global_ordinal_type nx   = 100;   clp.setOption("nx",                       &nx, "mesh size in x direction");
      global_ordinal_type ny   = 100;   clp.setOption("ny",                       &ny, "mesh size in y direction");
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

      //! [2DLaplacianOperator begin]
      // Create node map (equals dof map, since one dof per node)
      RCP<const Xpetra::Map<local_ordinal_type, global_ordinal_type, node_type>> xpetra_map = Galeri::Xpetra::CreateMap<local_ordinal_type, global_ordinal_type, node_type>(Xpetra::UseTpetra, "Cartesian2D", comm, galeriList);
      RCP<const map_type> nodeMap = Xpetra::toTpetra(xpetra_map);
      RCP<const map_type> dofMap = nodeMap;

      // Create coordinates
      RCP<const RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<scalar_type,local_ordinal_type,global_ordinal_type,map_type,RealValuedMultiVector>("2D", nodeMap, galeriList);

      // Create the matrix
      RCP<Galeri::Xpetra::Problem<map_type,crs_matrix_type,multi_vector_type>> galeriProblem =
          Galeri::Xpetra::BuildProblem<scalar_type,local_ordinal_type,global_ordinal_type,map_type,crs_matrix_type,multi_vector_type>("Laplace2D", dofMap, galeriList);
      RCP<crs_matrix_type> matrix = galeriProblem->BuildMatrix();
      // //! [2DLaplacianOperator end]

      // Some safety checks to see, if Galeri delivered valid output
      TEUCHOS_ASSERT(!nodeMap.is_null());
      TEUCHOS_ASSERT(!dofMap.is_null());
      TEUCHOS_ASSERT(!coordinates.is_null());
      TEUCHOS_ASSERT(!matrix.is_null());

      //! [RhsAndSolutionVector begin]
      // Create right-hand side (with all ones)
      RCP<multi_vector_type> B = rcp(new multi_vector_type(dofMap, 1, true));
      B->putScalar(one);

      // Initilize solution vector with random values
      RCP<multi_vector_type> X = rcp(new multi_vector_type(dofMap, 1, true));
      STS::seedrandom(100);
      X->randomize();
      //! [RhsAndSolutionVector end]

      //! [BuildNullSpaceVector begin]
      // Build null space vector for a scalar problem, i.e. vector with all ones
      RCP<multi_vector_type> nullspace = rcp(new multi_vector_type(dofMap, 1));
      nullspace->putScalar(one);
      //! [BuildNullSpaceVector end]

      // ================================
      // Preconditioner construction
      // ================================

      // Read MueLu parameter list from xml file
      RCP<ParameterList> mueluParams =  Teuchos::getParametersFromXmlFile(xmlFileName);

      // Register nullspace as user data in the MueLu parameter list
      ParameterList& userDataList = mueluParams->sublist("user data");
      userDataList.set<RCP<multi_vector_type>>("Nullspace", nullspace);

      // Create the MueLu preconditioner based on the Tpetra stack
      RCP<MueLu::TpetraOperator<scalar_type,local_ordinal_type,global_ordinal_type,node_type>> mueLuPreconditioner =
          MueLu::CreateTpetraPreconditioner(Teuchos::rcp_dynamic_cast<operator_type>(matrix), *mueluParams);

      // Generate exact solution using a direct solver
      //! [ExactSolutionVector begin]
      RCP<multi_vector_type> exactSolution = rcp(new multi_vector_type(dofMap, 1, true));
      //! [ExactSolutionVector end]
      {
        exactSolution->update(0.0, *X, 1.0);

        RCP<Amesos2::Solver<crs_matrix_type,multi_vector_type>> directSolver = Amesos2::create<crs_matrix_type,multi_vector_type>("KLU2", matrix, X, B);
        directSolver->solve();
      }

      //! [MueLuHierarchyAsPreconditionerWithinBelos begin]
      // Solve Ax = b using AMG as a preconditioner in Belos
      RCP<multi_vector_type> precSolVec = rcp(new multi_vector_type(dofMap, 1, true));
      {
        // Set initial guess
        precSolVec->update(0.0, *X, 1.0);

        // Construct a Belos LinearProblem object and hand-in the MueLu preconditioner
        RCP<Belos::LinearProblem<scalar_type,multi_vector_type,operator_type>> belosProblem =
            rcp(new Belos::LinearProblem<scalar_type,multi_vector_type,operator_type>(matrix, precSolVec, B));
        belosProblem->setLeftPrec(mueLuPreconditioner);

        bool set = belosProblem->setProblem();
        if (set == false) {
          fancyout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
          return EXIT_FAILURE;
        }

        // Belos parameter list
        RCP<ParameterList> belosList = Teuchos::parameterList();
        belosList->set("Maximum Iterations", 50); // Maximum number of iterations allowed
        belosList->set("Convergence Tolerance", tol); // Relative convergence tolerance requested
        belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
        belosList->set("Output Frequency", 1);
        belosList->set("Output Style", Belos::Brief);

        // Create an iterative solver manager
        Belos::SolverFactory<scalar_type,multi_vector_type,operator_type> solverFactory;
        RCP<Belos::SolverManager<scalar_type,multi_vector_type,operator_type>> solver = solverFactory.create("Block GMRES", belosList);
        solver->setProblem(belosProblem);

        // Perform solve
        Belos::ReturnType retStatus = Belos::Unconverged;
        retStatus = solver->solve();

        // Get the number of iterations for this solve
        fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

        // Check convergence status
        if (retStatus != Belos::Converged)
          fancyout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
        else
          fancyout << std::endl << "SUCCESS:  Belos converged!" << std::endl;

      }
      //! [MueLuHierarchyAsPreconditionerWithinBelos end]

      //! [UseMultigridHierarchyAsSolver begin]
      // Solve Ax = b using AMG as a solver
      RCP<multi_vector_type> multigridSolVec = rcp(new multi_vector_type(dofMap, 1, true));
      {
        // Set initial guess
        multigridSolVec->update(0.0, *X, 1.0);

        // Extract the underlying MueLu hierarchy
        RCP<MueLu::Hierarchy<scalar_type,local_ordinal_type,global_ordinal_type,node_type>> hierarchy = mueLuPreconditioner->GetHierarchy();

        // Configure MueLu to be used as solver
        hierarchy->IsPreconditioner(false);

        // Solve
        hierarchy->Iterate(*Xpetra::toXpetra(B), *Xpetra::toXpetra(multigridSolVec), mgridSweeps);
      }
      //! [UseMultigridHierarchyAsSolver end]

    } // end of Tpetra::ScopeGuard
  success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);


  //   // ================================
  //   // Preconditioner construction
  //   // ================================

  //   ParameterListInterpreter mueLuFactory(xmlFileName, *comm);

  //   RCP<Hierarchy> hierarchy = mueLuFactory.CreateHierarchy();
  //   hierarchy->IsPreconditioner(true);
  //   hierarchy->GetLevel(0)->Set("A", matrix);
  //   hierarchy->GetLevel(0)->Set("Nullspace", nullspace);
  //   // hierarchy->GetLevel(0)->Set("Coordinates", coordinates);

  //   mueLuFactory.SetupHierarchy(*hierarchy);

  //   // Generate exact solution using a direct solver
  //   //! [ExactSolutionVector begin]
  //   RCP<Vector> exactSolution = VectorFactory::Build(dofMap, true);
  //   //! [EpetraSolutionVector end]
  //   {
  //     exactSolution->update(1.0, *X, 1.0);

  //     // Amesos2 works with Tpetra objects directly, so we need to convert
  //     RCP<Tpetra_CrsMatrix> tA = Utilities::Op2NonConstTpetraCrs(matrix);
  //     RCP<Tpetra_MultiVector> tX = Utilities::MV2NonConstTpetraMV2(*X);
  //     RCP<Tpetra_MultiVector> tB = Utilities::MV2NonConstTpetraMV2(*B);

  //     RCP<Amesos2::Solver<Tpetra_CrsMatrix,Tpetra_MultiVector>> directSolver = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>("KLU2", tA, tX, tB);
  //     directSolver->solve();
  //   }

  //   //! [MueLuHierarchyAsPreconditionerWithinBelos begin]
  //   // Solve Ax = b using AMG as a preconditioner in Belos
  //   RCP<Vector> precSolVec = VectorFactory::Build(dofMap, true);
  //   {
  //     // Set initial guess
  //     precSolVec->update(0.0, *X, 1.0);

  //     // Configure MueLu to be used as a preconditioner
  //     hierarchy->IsPreconditioner(true);

  //     // Turn the Xpetra::Matrix object into a Belos operator
  //     Teuchos::RCP<OP> belosOp = Teuchos::rcp(new Belos::XpetraOp<SC,LO,GO,NO>(matrix));

  //     // Turns a MueLu::Hierarchy object into a Belos operator
  //     Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC,LO,GO,NO>(hierarchy));

  //     // Construct a Belos LinearProblem object
  //     RCP<Belos::LinearProblem<SC,MV,OP> > belosProblem =
  //         rcp(new Belos::LinearProblem<SC,MV,OP>(belosOp, precSolVec, B));

  //     bool set = belosProblem->setProblem();
  //     if (set == false) {
  //       fancyout << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
  //       return EXIT_FAILURE;
  //     }

  //     // Belos parameter list
  //     RCP<ParameterList> belosList = Teuchos::parameterList();
  //     belosList->set("Maximum Iterations", 50); // Maximum number of iterations allowed
  //     belosList->set("Convergence Tolerance", tol); // Relative convergence tolerance requested
  //     belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  //     belosList->set("Output Frequency", 1);
  //     belosList->set("Output Style", Belos::Brief);

  //     // Create an iterative solver manager
  //     Belos::SolverFactory<SC,MV,OP> solverFactory;
  //     RCP<Belos::SolverManager<SC,MV,OP>> solver = solverFactory.create("Block GMRES", belosList);
  //     solver->setProblem(belosProblem);

  //     // Perform solve
  //     Belos::ReturnType retStatus = Belos::Unconverged;
  //     retStatus = solver->solve();

  //     // Get the number of iterations for this solve
  //     fancyout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

  //     // Check convergence status
  //     if (retStatus != Belos::Converged)
  //       fancyout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
  //     else
  //       fancyout << std::endl << "SUCCESS:  Belos converged!" << std::endl;

  //   }
  //   //! [MueLuHierarchyAsPreconditionerWithinBelos end]

  //   //! [UseMultigridHierarchyAsSolver begin]
  //   // Solve Ax = b using AMG as a solver
  //   RCP<Vector> multigridSolVec = VectorFactory::Build(dofMap, true);
  //   {
  //     // Set initial guess
  //     multigridSolVec->update(0.0, *X, 1.0);

  //     // Configure MueLu to be used as solver
  //     hierarchy->IsPreconditioner(false);

  //     // Solve
  //     hierarchy->Iterate(*B, *multigridSolVec, mgridSweeps);
  //   }
  //   //! [UseMultigridHierarchyAsSolver end]

  //   success = true;
  // }
  // TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // main
