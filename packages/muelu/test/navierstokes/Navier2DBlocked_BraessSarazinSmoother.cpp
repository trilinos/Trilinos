// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Navier2D_blocked_bssmoother.cpp
 *
 *  Created on: Jun 18, 2012
 *      Author: wiesner
 */

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_StridedMap.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_BraessSarazinSmoother.hpp"
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_TopSmootherFactory.hpp"
#include "MueLu_HierarchyUtils.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include "Navier2D_Helpers.h"

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 *  using block matrices
 */

int main(int argc, char* argv[]) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLuTests;
  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
    RCP<FancyOStream> out      = fancyOStream(rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

    // Timing
    Time myTime("global");
    TimeMonitor MM(myTime);

    // read in input parameters

    // default parameters
    LO BS_nSweeps           = 100;
    Scalar BS_omega         = 1.7;
    LO SC_nSweeps           = 1;
    Scalar SC_omega         = 1.0;
    int SC_bUseDirectSolver = 0;

    // Note: use --help to list available options.
    CommandLineProcessor clp(false);
    clp.setOption("BraessSarazin_sweeps", &BS_nSweeps, "number of sweeps with BraessSarazin smoother");
    clp.setOption("BraessSarazin_omega", &BS_omega, "scaling factor for BraessSarazin smoother");
    clp.setOption("SchurComp_sweeps", &SC_nSweeps, "number of sweeps for BraessSarazin internal SchurComp solver/smoother (GaussSeidel)");
    clp.setOption("SchurComp_omega", &SC_omega, "damping parameter for BraessSarazin internal SchurComp solver/smoother (GaussSeidel)");
    clp.setOption("SchurComp_solver", &SC_bUseDirectSolver, "if 1: use direct solver for SchurComp equation, otherwise use GaussSeidel smoother (=default)");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    int globalNumDofs = 1500;  // used for the maps
    // int nDofsPerNode = 3;      // used for generating the fine level null-space

    // build strided maps
    // striding information: 2 velocity dofs and 1 pressure dof = 3 dofs per node
    std::vector<size_t> stridingInfo;
    stridingInfo.push_back(2);
    stridingInfo.push_back(1);

    /////////////////////////////////////// build strided maps
    // build strided maps:
    // xstridedfullmap: full map (velocity and pressure dof gids), continous
    // xstridedvelmap: only velocity dof gid maps (i.e. 0,1,3,4,6,7...)
    // xstridedpremap: only pressure dof gid maps (i.e. 2,5,8,...)
    Xpetra::UnderlyingLib lib       = Xpetra::UseEpetra;
    RCP<StridedMap> xstridedfullmap = StridedMapFactory::Build(lib, globalNumDofs, 0, stridingInfo, comm, -1);
    RCP<StridedMap> xstridedvelmap  = StridedMapFactory::Build(xstridedfullmap, 0);
    RCP<StridedMap> xstridedpremap  = StridedMapFactory::Build(xstridedfullmap, 1);

    /////////////////////////////////////// transform Xpetra::Map objects to Epetra
    // this is needed for our splitting routine
    const RCP<const Epetra_Map> fullmap = rcpFromRef(Xpetra::toEpetra(*xstridedfullmap));
    RCP<const Epetra_Map> velmap        = rcpFromRef(Xpetra::toEpetra(*xstridedvelmap));
    RCP<const Epetra_Map> premap        = rcpFromRef(Xpetra::toEpetra(*xstridedpremap));

    /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

    // read in problem
    Epetra_CrsMatrix* ptrA    = 0;
    Epetra_Vector* ptrf       = 0;
    Epetra_MultiVector* ptrNS = 0;

    *out << "Reading matrix market file" << std::endl;

    EpetraExt::MatrixMarketFileToCrsMatrix("A_re1000_5932.txt", *fullmap, *fullmap, *fullmap, ptrA);
    EpetraExt::MatrixMarketFileToVector("b_re1000_5932.txt", *fullmap, ptrf);

    RCP<Epetra_CrsMatrix> epA    = rcp(ptrA);
    RCP<Epetra_Vector> epv       = rcp(ptrf);
    RCP<Epetra_MultiVector> epNS = rcp(ptrNS);

    /////////////////////////////////////// split system into 2x2 block system

    *out << "Split matrix into 2x2 block matrix" << std::endl;

    // split fullA into A11,..., A22
    RCP<Epetra_CrsMatrix> A11;
    RCP<Epetra_CrsMatrix> A12;
    RCP<Epetra_CrsMatrix> A21;
    RCP<Epetra_CrsMatrix> A22;

    if (SplitMatrix2x2(epA, *velmap, *premap, A11, A12, A21, A22) == false)
      *out << "Problem with splitting matrix" << std::endl;

    /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

    // build Xpetra objects from Epetra_CrsMatrix objects
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA11 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A11));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA12 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A12));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA21 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A21));
    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA22 = rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A22));

    /////////////////////////////////////// generate MapExtractor object

    std::vector<RCP<const Xpetra::Map<LO, GO, Node> > > xmaps;

    xmaps.push_back(xstridedvelmap);
    xmaps.push_back(xstridedpremap);

    RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(xstridedfullmap, xmaps);

    /////////////////////////////////////// build blocked transfer operator
    // using the map extractor
    RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(map_extractor, map_extractor, 10));
    bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA11)));
    bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA12)));
    bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA21)));
    bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA22)));

    bOp->fillComplete();
    //////////////////////////////////////////////////////// finest Level
    RCP<MueLu::Level> Finest = rcp(new Level());
    Finest->setDefaultVerbLevel(VERB_NONE);
    Finest->Set("A", rcp_dynamic_cast<Matrix>(bOp));

    ///////////////////////////////////
    // Test Braess Sarazin Smoother as a solver

    *out << "Test: Creating Braess Sarazin Smoother" << std::endl;
    *out << "Test: Omega for BraessSarazin = " << BS_omega << std::endl;
    *out << "Test: Number of sweeps for BraessSarazin = " << BS_nSweeps << std::endl;
    *out << "Test: Omega for Schur Complement solver= " << SC_omega << std::endl;
    *out << "Test: Number of Schur Complement solver= " << SC_nSweeps << std::endl;
    *out << "Test: Setting up Braess Sarazin Smoother" << std::endl;

    // define BraessSarazin Smoother with BS_nSweeps and BS_omega as scaling factor
    // AFact_ = null (= default) for the 2x2 blocked operator
    RCP<BraessSarazinSmoother> BraessSarazinSm = rcp(new BraessSarazinSmoother());
    BraessSarazinSm->SetParameter("Sweeps", Teuchos::ParameterEntry(BS_nSweeps));
    BraessSarazinSm->SetParameter("Damping factor", Teuchos::ParameterEntry(BS_omega));

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(BraessSarazinSm));

    /*note that omega must be the same in the SchurComplementFactory and in the BraessSarazinSmoother*/
    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // and the scaling/damping factor omega that is used for BraessSarazin
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    RCP<SchurComplementFactory> SFact = rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", ParameterEntry(BS_omega));
    SFact->SetFactory("A", MueLu::NoFactory::getRCP());

    // define smoother/solver for BraessSarazin
    RCP<SmootherPrototype> smoProtoSC = null;
    if (SC_bUseDirectSolver != 1) {
      // Smoother Factory, using SFact as a factory for A
      std::string ifpackSCType;
      ParameterList ifpackSCList;
      ifpackSCList.set("relaxation: sweeps", SC_nSweeps);
      ifpackSCList.set("relaxation: damping factor", SC_omega);
      ifpackSCType = "RELAXATION";
      ifpackSCList.set("relaxation: type", "Gauss-Seidel");
      smoProtoSC = rcp(new TrilinosSmoother(ifpackSCType, ifpackSCList, 0));
      smoProtoSC->SetFactory("A", SFact);
    } else {
      ParameterList ifpackDSList;
      std::string ifpackDSType;
      smoProtoSC = rcp(new DirectSolver(ifpackDSType, ifpackDSList));
      smoProtoSC->SetFactory("A", SFact);
    }

    RCP<SmootherFactory> SmooSCFact = rcp(new SmootherFactory(smoProtoSC));

    // define temporary FactoryManager that is used as input for BraessSarazin smoother
    RCP<FactoryManager> MB = rcp(new FactoryManager());
    MB->SetFactory("A", SFact);              // SchurComplement operator for correction step (defined as "A")
    MB->SetFactory("Smoother", SmooSCFact);  // solver/smoother for correction step
    MB->SetFactory("PreSmoother", SmooSCFact);
    MB->SetFactory("PostSmoother", SmooSCFact);
    MB->SetIgnoreUserData(true);                // always use data from factories defined in factory manager
    BraessSarazinSm->AddFactoryManager(MB, 0);  // set temporary factory manager in BraessSarazin smoother

    // setup main factory manager
    RCP<FactoryManager> M = rcp(new FactoryManager());
    M->SetFactory("A", MueLu::NoFactory::getRCP());  // this is the 2x2 blocked operator
    M->SetFactory("Smoother", smootherFact);         // BraessSarazin block smoother
    M->SetFactory("PreSmoother", smootherFact);
    M->SetFactory("PostSmoother", smootherFact);

    MueLu::SetFactoryManager SFMCoarse(Finest, M);
    Finest->Request(MueLu::TopSmootherFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>(M, "Smoother"));

    // call setup (= extract blocks and extract diagonal of F)
    BraessSarazinSm->Setup(*Finest);

    RCP<MultiVector> xtest = MultiVectorFactory::Build(xstridedfullmap, 1);
    xtest->putScalar((SC)0.0);

    RCP<Vector> xR = rcp(new Xpetra::EpetraVectorT<GO, Node>(epv));
    // calculate initial (absolute) residual
    Array<ScalarTraits<SC>::magnitudeType> norms(1);

    xR->norm2(norms);
    *out << "Test: ||x_0|| = " << norms[0] << std::endl;
    *out << "Test: Applying Braess-Sarazin Smoother" << std::endl;
    *out << "Test: START DATA" << std::endl;
    *out << "iterations\tVelocity_residual\tPressure_residual" << std::endl;
    BraessSarazinSm->Apply(*xtest, *xR);
    xtest->norm2(norms);
    *out << "Test: ||x_1|| = " << norms[0] << std::endl;

    Array<ScalarTraits<Scalar>::magnitudeType> test = MueLu::Utilities<Scalar, LO, GO, Node>::ResidualNorm(*bOp, *xtest, *xR);
    *out << "residual norm: " << test[0] << std::endl;

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
#else
  std::cout << "Epetra (and/or EpetraExt) are not available. Skip test." << std::endl;
  return EXIT_SUCCESS;
#endif  // #if defined(HAVE_MUELU_SERIAL) && defined(HAVE_MUELU_EPETRA)
}
