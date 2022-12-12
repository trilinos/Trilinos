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

#include <MueLu_ConfigDefs.hpp>
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp
#include <MueLu_Utilities.hpp>

#include <Teuchos_StandardCatchMacros.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib& lib, int argc, char *argv[])
{
#include <MueLu_UseShortNames.hpp>

  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  using Tpetra_CrsMatrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Tpetra_MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  using MV = MultiVector;
  using OP = Belos::OperatorT<MV>;


  //! [CommunicatorObject begin]
  // Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  try {
    RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
    int MyPID = comm->getRank();
    int NumProc = comm->getSize();

    //! [CommunicatorObject end]
    // ================================
    // Convenient definitions
    // ================================
    using STS = Teuchos::ScalarTraits<SC>;
    using magnitude_type = typename Teuchos::ScalarTraits<Scalar>::magnitudeType;
    using real_type = typename STS::coordinateType;
    using RealValuedMultiVector = Xpetra::MultiVector<real_type,LO,GO,NO>;
    const SC zero = Teuchos::ScalarTraits<SC>::zero();
    const SC one = Teuchos::ScalarTraits<SC>::one();

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
    // galeriList.set("lx", 1.0); // length of x-axis
    // galeriList.set("ly", 1.0); // length of y-axis

    //! [2DLaplacianOperator begin]
    // Create node map (equals dof map, since one dof per node)
    RCP<Map> nodeMap = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian2D", comm, galeriList);
    RCP<Map> dofMap = nodeMap;

    // Create coordinates
    RCP<RealValuedMultiVector> coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<double,LO,GO,Map,RealValuedMultiVector>("2D", nodeMap, galeriList);

    // Create the matrix
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > galeriProblem =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>("Laplace2D", dofMap, galeriList);
    RCP<Matrix> matrix = galeriProblem->BuildMatrix();
    matrix->SetFixedBlockSize(1);
    //! [2DLaplacianOperator end]

    // Some safety checks to see, if Galeri delived valid output
    TEUCHOS_ASSERT(!nodeMap.is_null());
    TEUCHOS_ASSERT(!dofMap.is_null());
    TEUCHOS_ASSERT(!coordinates.is_null());
    TEUCHOS_ASSERT(!matrix.is_null());

    //! [RhsAndSolutionVector begin]
    // Create right-hand side (with all ones)
    RCP<MultiVector> B = MultiVectorFactory::Build(dofMap, 1, true);
    B->putScalar(one);

    // Initilize solution vector with random values
    RCP<MultiVector> X = MultiVectorFactory::Build(dofMap, 1);
    X->setSeed(100);
    X->randomize();
    //! [RhsAndSolutionVector end]

    //! [BuildNullSpaceVector begin]
    // Build null space vector for a scalar problem, i.e. vector with all ones
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(dofMap, 1);
    nullspace->putScalar(one);
    //! [BuildNullSpaceVector end]

    // ================================
    // Preconditioner construction
    // ================================

    ParameterListInterpreter mueLuFactory(xmlFileName, *comm);

    RCP<Hierarchy> hierarchy = mueLuFactory.CreateHierarchy();
    hierarchy->IsPreconditioner(true);
    hierarchy->GetLevel(0)->Set("A", matrix);
    hierarchy->GetLevel(0)->Set("Nullspace", nullspace);
    // hierarchy->GetLevel(0)->Set("Coordinates", coordinates);

    mueLuFactory.SetupHierarchy(*hierarchy);

    // Generate exact solution using a direct solver
    //! [ExactSolutionVector begin]
    RCP<Vector> exactSolution = VectorFactory::Build(dofMap, true);
    //! [EpetraSolutionVector end]
    {
      exactSolution->update(1.0, *X, 1.0);

      // Amesos2 works with Tpetra objects directly, so we need to convert
      RCP<Tpetra_CrsMatrix> tA = Utilities::Op2NonConstTpetraCrs(matrix);
      RCP<Tpetra_MultiVector> tX = Utilities::MV2NonConstTpetraMV2(*X);
      RCP<Tpetra_MultiVector> tB = Utilities::MV2NonConstTpetraMV2(*B);

      RCP<Amesos2::Solver<Tpetra_CrsMatrix,Tpetra_MultiVector>> directSolver = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>("KLU2", tA, tX, tB);
      directSolver->solve();
    }

    //! [MueLuHierarchyAsPreconditionerWithinBelos begin]
    // Solve Ax = b using AMG as a preconditioner in Belos
    RCP<Vector> precSolVec = VectorFactory::Build(dofMap, true);
    {
      // Set initial guess
      precSolVec->update(0.0, *X, 1.0);

      // Configure MueLu to be used as a preconditioner
      hierarchy->IsPreconditioner(true);

      // Turn the Xpetra::Matrix object into a Belos operator
      Teuchos::RCP<OP> belosOp = Teuchos::rcp(new Belos::XpetraOp<SC,LO,GO,NO>(matrix));

      // Turns a MueLu::Hierarchy object into a Belos operator
      Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC,LO,GO,NO>(hierarchy));

      // Construct a Belos LinearProblem object
      RCP<Belos::LinearProblem<SC,MV,OP> > belosProblem =
          rcp(new Belos::LinearProblem<SC,MV,OP>(belosOp, precSolVec, B));

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
      Belos::SolverFactory<SC,MV,OP> solverFactory;
      RCP<Belos::SolverManager<SC,MV,OP>> solver = solverFactory.create("Block GMRES", belosList);
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
    RCP<Vector> multigridSolVec = VectorFactory::Build(dofMap, true);
    {
      // Set initial guess
      multigridSolVec->update(0.0, *X, 1.0);

      // Configure MueLu to be used as solver
      hierarchy->IsPreconditioner(false);

      // Solve
      hierarchy->Iterate(*B, *multigridSolVec, mgridSweeps);
    }
    //! [UseMultigridHierarchyAsSolver end]

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
} // main_

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}
