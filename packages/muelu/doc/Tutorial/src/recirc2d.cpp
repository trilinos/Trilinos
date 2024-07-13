// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>

#include <Teuchos_StandardCatchMacros.hpp>

#include <Epetra_RowMatrix.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>

// Galeri
#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>
#include <Galeri_Utils.h>
//

//#include <Teuchos_ArrayRCP.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_ParameterListInterpreter.hpp>  // TODO: move into MueLu.hpp

#include <MueLu_Utilities.hpp>

#include <MueLu_MutuallyExclusiveTime.hpp>

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#if defined(HAVE_MUELU_EPETRA)
#include <MueLu_EpetraOperator.hpp>

// prescribe types
// run plain Epetra
typedef double Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef Xpetra::EpetraNode Node;  // Epetra needs SerialNode
#endif

int main(int argc, char* argv[]) {
#if defined(HAVE_MUELU_EPETRA)
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;  // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    int MyPID                           = comm->getRank();
    int NumProc                         = comm->getSize();

    const Teuchos::RCP<Epetra_Comm> epComm = Teuchos::rcp_const_cast<Epetra_Comm>(Xpetra::toEpetra(comm));

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    // SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& fancyout  = *fancy;
    fancyout.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    Teuchos::CommandLineProcessor clp(false);
    GO nx = 100, ny = 100;
    clp.setOption("nx", &nx, "mesh size in x direction");
    clp.setOption("ny", &ny, "mesh size in y direction");
    std::string xmlFileName = "xml/s3a.xml";
    clp.setOption("xml", &xmlFileName, "read parameters from a file. Otherwise, this example uses by default 'tutorial1a.xml'");
    int mgridSweeps = 1;
    clp.setOption("mgridSweeps", &mgridSweeps, "number of multigrid sweeps within Multigrid solver.");
    std::string printTimings = "no";
    clp.setOption("timings", &printTimings, "print timings to screen [yes/no]");
    double tol = 1e-12;
    clp.setOption("tol", &tol, "solver convergence tolerance");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    // =========================================================================
    // Problem construction
    // =========================================================================
    RCP<TimeMonitor> globalTimeMonitor = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time"))), tm;

    comm->barrier();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));

    Teuchos::ParameterList GaleriList;
    GaleriList.set("nx", nx);
    GaleriList.set("ny", ny);
    GaleriList.set("mx", epComm->NumProc());
    GaleriList.set("my", 1);
    GaleriList.set("lx", 1.0);  // length of x-axis
    GaleriList.set("ly", 1.0);  // length of y-axis
    GaleriList.set("diff", 1e-5);
    GaleriList.set("conv", 1.0);

    // create map
    Teuchos::RCP<Epetra_Map> epMap = Teuchos::rcp(Galeri::CreateMap("Cartesian2D", *epComm, GaleriList));

    // create coordinates
    Teuchos::RCP<Epetra_MultiVector> epCoord = Teuchos::rcp(Galeri::CreateCartesianCoordinates("2D", epMap.get(), GaleriList));

    // create matrix
    Teuchos::RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(Galeri::CreateCrsMatrix("Recirc2D", epMap.get(), GaleriList));

    // Epetra -> Xpetra
    Teuchos::RCP<CrsMatrix> exA         = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<int, Node>(epA));
    Teuchos::RCP<CrsMatrixWrap> exAWrap = Teuchos::rcp(new CrsMatrixWrap(exA));

    RCP<Matrix> A = Teuchos::rcp_dynamic_cast<Matrix>(exAWrap);
    int numPDEs   = 1;
    A->SetFixedBlockSize(numPDEs);

    // set rhs and solution vector
    RCP<Epetra_Vector> B = Teuchos::rcp(new Epetra_Vector(*epMap));
    RCP<Epetra_Vector> X = Teuchos::rcp(new Epetra_Vector(*epMap));
    B->PutScalar(1.0);
    X->PutScalar(0.0);

    // Epetra -> Xpetra
    RCP<Vector> xB          = Teuchos::rcp(new Xpetra::EpetraVectorT<int, Node>(B));
    RCP<Vector> xX          = Teuchos::rcp(new Xpetra::EpetraVectorT<int, Node>(X));
    RCP<MultiVector> coords = Teuchos::rcp(new Xpetra::EpetraMultiVectorT<int, Node>(epCoord));

    xX->setSeed(100);
    xX->randomize();

    // build null space vector
    RCP<const Map> map         = A->getRowMap();
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, numPDEs);

    for (int i = 0; i < numPDEs; ++i) {
      Teuchos::ArrayRCP<Scalar> nsValues = nullspace->getDataNonConst(i);
      int numBlocks                      = nsValues.size() / numPDEs;
      for (int j = 0; j < numBlocks; ++j) {
        nsValues[j * numPDEs + i] = 1.0;
      }
    }

    comm->barrier();
    tm = Teuchos::null;

    fancyout << "========================================================\nGaleri complete.\n========================================================" << std::endl;

    // =========================================================================
    // Preconditioner construction
    // =========================================================================
    comm->barrier();
    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1.5 - MueLu read XML")));
    ParameterListInterpreter mueLuFactory(xmlFileName, *comm);
    comm->barrier();
    tm = Teuchos::null;

    tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 2 - MueLu Setup")));

    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();

    H->GetLevel(0)->Set("A", A);
    H->GetLevel(0)->Set("Nullspace", nullspace);
    H->GetLevel(0)->Set("Coordinates", coords);

    mueLuFactory.SetupHierarchy(*H);

    comm->barrier();
    tm = Teuchos::null;

    // =========================================================================
    // System solution (Ax = b)
    // =========================================================================

    //
    // generate exact solution using a direct solver
    //
    RCP<Epetra_Vector> exactLsgVec = rcp(new Epetra_Vector(X->Map()));
    {
      fancyout << "========================================================\nCalculate exact solution." << std::endl;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 3 - direct solve")));
      exactLsgVec->PutScalar(0.0);
      exactLsgVec->Update(1.0, *X, 1.0);
      Epetra_LinearProblem epetraProblem(epA.get(), exactLsgVec.get(), B.get());

      Amesos amesosFactory;
      RCP<Amesos_BaseSolver> rcp_directSolver = Teuchos::rcp(amesosFactory.Create("Amesos_Klu", epetraProblem));
      rcp_directSolver->SymbolicFactorization();
      rcp_directSolver->NumericFactorization();
      rcp_directSolver->Solve();

      comm->barrier();
      tm = Teuchos::null;
    }

    //
    // Solve Ax = b using AMG as a preconditioner in AztecOO
    //

    RCP<Epetra_Vector> precLsgVec = rcp(new Epetra_Vector(X->Map()));
    {
      fancyout << "========================================================\nUse multigrid hierarchy as preconditioner within CG." << std::endl;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 4 - AMG as preconditioner")));

      precLsgVec->PutScalar(0.0);
      precLsgVec->Update(1.0, *X, 1.0);
      Epetra_LinearProblem epetraProblem(epA.get(), precLsgVec.get(), B.get());

      AztecOO aztecSolver(epetraProblem);
      aztecSolver.SetAztecOption(AZ_solver, AZ_gmres);

      MueLu::EpetraOperator aztecPrec(H);
      aztecSolver.SetPrecOperator(&aztecPrec);

      int maxIts = 50;

      aztecSolver.Iterate(maxIts, tol);

      comm->barrier();
      tm = Teuchos::null;
    }

    //////////////////

    // use multigrid hierarchy as solver
    RCP<Vector> mgridLsgVec = VectorFactory::Build(map);
    mgridLsgVec->putScalar(0.0);
    {
      fancyout << "========================================================\nUse multigrid hierarchy as solver." << std::endl;
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 5 - Multigrid Solve")));
      mgridLsgVec->update(1.0, *xX, 1.0);
      H->IsPreconditioner(false);
      H->Iterate(*xB, *mgridLsgVec, mgridSweeps);
      comm->barrier();
      tm = Teuchos::null;
    }

    //////////////////

    fancyout << "========================================================\nExport results.\n========================================================" << std::endl;
    std::ofstream myfile;
    std::stringstream ss;
    ss << "example" << MyPID << ".txt";
    myfile.open(ss.str().c_str());

    //////////////////

    // loop over all procs
    for (int iproc = 0; iproc < NumProc; iproc++) {
      if (MyPID == iproc) {
        int NumVectors1               = 2;
        int NumMyElements1            = epCoord->Map().NumMyElements();
        int MaxElementSize1           = epCoord->Map().MaxElementSize();
        int* FirstPointInElementList1 = NULL;
        if (MaxElementSize1 != 1) FirstPointInElementList1 = epCoord->Map().FirstPointInElementList();
        double** A_Pointers = epCoord->Pointers();

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
          for (int ii = 0; ii < epCoord->Map().ElementSize(i); ii++) {
            int iii;
            myfile.width(10);
            myfile << MyPID;
            myfile << "    ";
            myfile.width(10);
            if (MaxElementSize1 == 1) {
              if (epCoord->Map().GlobalIndicesInt()) {
                int* MyGlobalElements1 = epCoord->Map().MyGlobalElements();
                myfile << MyGlobalElements1[i] << "    ";
              }

              iii = i;
            } else {
              if (epCoord->Map().GlobalIndicesInt()) {
                int* MyGlobalElements1 = epCoord->Map().MyGlobalElements();
                myfile << MyGlobalElements1[i] << "/" << ii << "    ";
              }

              iii = FirstPointInElementList1[i] + ii;
            }
            for (int j = 0; j < NumVectors1; j++) {
              myfile.width(20);
              myfile << A_Pointers[j][iii];
            }

            myfile.precision(18);  // set high precision for output

            // add solution vector entry
            myfile.width(25);
            myfile << (*exactLsgVec)[iii];

            // add preconditioned solution vector entry
            myfile.width(25);
            myfile << (*precLsgVec)[iii];

            Teuchos::ArrayRCP<SC> mgridLsgVecData = mgridLsgVec->getDataNonConst(0);
            myfile.width(25);
            myfile << mgridLsgVecData[iii];

            myfile.precision(6);  // set default precision
            myfile << std::endl;
          }
        }  // end loop over all lines on current proc
        myfile << std::flush;

        // syncronize procs
        comm->barrier();
        comm->barrier();
        comm->barrier();

      }  // end myProc
    }

    ////////////
    myfile.close();

    comm->barrier();
    tm                = Teuchos::null;
    globalTimeMonitor = Teuchos::null;

    if (printTimings == "yes") {
      TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, false, true, false, Teuchos::Union, "", true);
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
#else
  return EXIT_SUCCESS;
#endif  // #if defined(HAVE_MUELU_EPETRA) and defined(HAVE_MUELU_SERIAL)
}  // main
