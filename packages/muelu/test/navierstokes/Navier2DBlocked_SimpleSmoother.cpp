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
#include "MueLu_SimpleSmoother.hpp"
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_TopSmootherFactory.hpp"
#include "MueLu_HierarchyUtils.hpp"
#include "MueLu_InverseApproximationFactory.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include "Navier2D_Helpers.h"

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 *  using block matrices
 */

int main(int argc, char* argv[]) {
#if defined(HAVE_MUELU_EPETRA)
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

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    RCP<Teuchos::FancyOStream> out      = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    out->setOutputToRootOnly(0);
    *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

    // Timing
    Teuchos::Time myTime("global");
    Teuchos::TimeMonitor MM(myTime);

    // read in input parameters

    // default parameters
    LO SIMPLE_nSweeps   = 600;
    Scalar SIMPLE_omega = 0.5;
    LO SC_nSweeps       = 10;
    Scalar SC_omega     = 1.0;
    LO PRED_nSweeps     = 3;
    Scalar PRED_omega   = 1.0;
    LO useSIMPLEC       = 0;

    int SC_bUseDirectSolver = 0;

    // Note: use --help to list available options.
    Teuchos::CommandLineProcessor clp(false);
    clp.setOption("SIMPLE_sweeps", &SIMPLE_nSweeps, "number of sweeps with SIMPLE smoother");
    clp.setOption("SIMPLE_omega", &SIMPLE_omega, "scaling factor for SIMPLE smoother");
    clp.setOption("Predict_sweeps", &PRED_nSweeps, "number of sweeps for SIMPLE internal velocity prediction smoother (GaussSeidel)");
    clp.setOption("Predict_omega", &PRED_omega, "damping parameter for SIMPLE internal velocity prediction smoother (GaussSeidel)");
    clp.setOption("SchurComp_sweeps", &SC_nSweeps, "number of sweeps for SIMPLE internal SchurComp solver/smoother (GaussSeidel)");
    clp.setOption("SchurComp_omega", &SC_omega, "damping parameter for SIMPLE internal SchurComp solver/smoother (GaussSeidel)");
    clp.setOption("SchurComp_solver", &SC_bUseDirectSolver, "if 1: use direct solver for SchurComp equation, otherwise use GaussSeidel smoother");
    clp.setOption("useSIMPLEC", &useSIMPLEC, "if 1: use SIMPLEC instead of SIMPLE (default = 0 (SIMPLE))");

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    int globalNumDofs = 1500;  // used for the maps

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
    Xpetra::UnderlyingLib lib             = Xpetra::UseEpetra;
    RCP<const StridedMap> xstridedfullmap = StridedMapFactory::Build(lib, globalNumDofs, 0, stridingInfo, comm, -1);
    RCP<const StridedMap> xstridedvelmap  = StridedMapFactory::Build(xstridedfullmap, 0);
    RCP<const StridedMap> xstridedpremap  = StridedMapFactory::Build(xstridedfullmap, 1);

    /////////////////////////////////////// transform Xpetra::Map objects to Epetra
    // this is needed for our splitting routine
    const RCP<const Epetra_Map> fullmap = Teuchos::rcpFromRef(Xpetra::toEpetra(*xstridedfullmap));
    RCP<const Epetra_Map> velmap        = Teuchos::rcpFromRef(Xpetra::toEpetra(*xstridedvelmap));
    RCP<const Epetra_Map> premap        = Teuchos::rcpFromRef(Xpetra::toEpetra(*xstridedpremap));

    /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

    // read in problem
    Epetra_CrsMatrix* ptrA    = 0;
    Epetra_Vector* ptrf       = 0;
    Epetra_MultiVector* ptrNS = 0;

    *out << "Reading matrix market file" << std::endl;

    EpetraExt::MatrixMarketFileToCrsMatrix("A_re1000_5932.txt", *fullmap, *fullmap, *fullmap, ptrA);
    EpetraExt::MatrixMarketFileToVector("b_re1000_5932.txt", *fullmap, ptrf);

    RCP<Epetra_CrsMatrix> epA    = Teuchos::rcp(ptrA);
    RCP<Epetra_Vector> epv       = Teuchos::rcp(ptrf);
    RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);

    /////////////////////////////////////// split system into 2x2 block system

    *out << "Split matrix into 2x2 block matrix" << std::endl;

    // split fullA into A11,..., A22
    Teuchos::RCP<Epetra_CrsMatrix> A11;
    Teuchos::RCP<Epetra_CrsMatrix> A12;
    Teuchos::RCP<Epetra_CrsMatrix> A21;
    Teuchos::RCP<Epetra_CrsMatrix> A22;

    if (SplitMatrix2x2(epA, *velmap, *premap, A11, A12, A21, A22) == false)
      *out << "Problem with splitting matrix" << std::endl;

    /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

    // build Xpetra objects from Epetra_CrsMatrix objects
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA11 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A11));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA12 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A12));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA21 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A21));
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > xA22 = Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<GO, Node>(A22));

    /////////////////////////////////////// generate MapExtractor object

    std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, Node> > > xmaps;

    xmaps.push_back(xstridedvelmap);
    xmaps.push_back(xstridedpremap);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar, LO, GO, Node>::Build(xstridedfullmap, xmaps);

    /////////////////////////////////////// build blocked transfer operator
    // using the map extractor
    Teuchos::RCP<Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsMatrix<Scalar, LO, GO, Node>(map_extractor, map_extractor, 10));
    bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA11)));
    bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA12)));
    bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA21)));
    bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(xA22)));

    bOp->fillComplete();
    //////////////////////////////////////////////////////// finest Level
    RCP<MueLu::Level> Finest = rcp(new Level());
    Finest->setDefaultVerbLevel(Teuchos::VERB_NONE);
    Finest->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));

    ///////////////////////////////////
    // Test Braess Sarazin Smoother as a solver

    *out << "Test: Creating SIMPLE Smoother" << std::endl;
    *out << "Test: Omega for SIMPLE = " << SIMPLE_omega << std::endl;
    *out << "Test: Number of sweeps for SIMPLE = " << SIMPLE_nSweeps << std::endl;
    *out << "Test: Omega for Schur Complement solver= " << SC_omega << std::endl;
    *out << "Test: Number of Schur Complement solver= " << SC_nSweeps << std::endl;
    *out << "Test: Setting up Braess Sarazin Smoother" << std::endl;

    // define SIMPLE Smoother with SIMPLE_nSweeps and SIMPLE_omega as scaling factor
    // AFact_ = Teuchos::null (= default) for the 2x2 blocked operator
    RCP<SimpleSmoother> SimpleSm = rcp(new SimpleSmoother());
    SimpleSm->SetParameter("Sweeps", Teuchos::ParameterEntry(SIMPLE_nSweeps));
    SimpleSm->SetParameter("Damping factor", Teuchos::ParameterEntry(SIMPLE_omega));
    if (useSIMPLEC == 1) SimpleSm->SetParameter("UseSIMPLEC", Teuchos::ParameterEntry(true));

    RCP<SmootherFactory> smootherFact = rcp(new SmootherFactory(SimpleSm));

    // define smoother for velocity prediction
    // RCP<SubBlockAFactory> A00Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
    RCP<SubBlockAFactory> A00Fact = rcp(new SubBlockAFactory());
    A00Fact->SetFactory("A", MueLu::NoFactory::getRCP());
    A00Fact->SetParameter("block row", Teuchos::ParameterEntry(0));
    A00Fact->SetParameter("block col", Teuchos::ParameterEntry(0));
    RCP<SmootherPrototype> smoProtoPredict = Teuchos::null;
    std::string ifpackPredictType;
    Teuchos::ParameterList ifpackPredictList;
    ifpackPredictList.set("relaxation: sweeps", PRED_nSweeps);
    ifpackPredictList.set("relaxation: damping factor", PRED_omega);
    ifpackPredictType = "RELAXATION";
    ifpackPredictList.set("relaxation: type", "Gauss-Seidel");
    smoProtoPredict = rcp(new TrilinosSmoother(ifpackPredictType, ifpackPredictList, 0));
    smoProtoPredict->SetFactory("A", A00Fact);
    RCP<SmootherFactory> SmooPredictFact = rcp(new SmootherFactory(smoProtoPredict));
    // define temporary FactoryManager that is used as input for BraessSarazin smoother
    RCP<FactoryManager> MPredict = rcp(new FactoryManager());
    MPredict->SetFactory("A", A00Fact);
    MPredict->SetFactory("Smoother", SmooPredictFact);  // solver/smoother for correction step
    MPredict->SetFactory("PreSmoother", SmooPredictFact);
    MPredict->SetFactory("PostSmoother", SmooPredictFact);
    MPredict->SetIgnoreUserData(true);                        // always use data from factories defined in factory manager
    SimpleSm->SetVelocityPredictionFactoryManager(MPredict);  // set temporary factory manager in BraessSarazin smoother

    // define SchurComplement Factory
    // SchurComp gets a RCP to AFact_ which has to be the 2x2 blocked operator
    // It stores the resulting SchurComplement operator as "A" generated by the SchurComplementFactory
    // Instead of F^{-1} it uses the approximation \hat{F}^{-1} with \hat{F} = diag(F)
    // InverseApproximation
    Teuchos::RCP<InverseApproximationFactory> AinvFact = Teuchos::rcp(new InverseApproximationFactory());
    AinvFact->SetFactory("A", A00Fact);
    if (useSIMPLEC == 1) AinvFact->SetParameter("inverse: approximation type", Teuchos::ParameterEntry(std::string("lumping")));

    RCP<SchurComplementFactory> SFact = Teuchos::rcp(new SchurComplementFactory());
    SFact->SetParameter("omega", Teuchos::ParameterEntry(1.0));  // for Simple, omega is always 1.0 in the SchurComplement
    SFact->SetFactory("A", MueLu::NoFactory::getRCP());

    // define smoother/solver for BraessSarazin
    RCP<SmootherPrototype> smoProtoSC = Teuchos::null;
    if (SC_bUseDirectSolver != 1) {
      // Smoother Factory, using SFact as a factory for A
      std::string ifpackSCType;
      Teuchos::ParameterList ifpackSCList;
      ifpackSCList.set("relaxation: sweeps", SC_nSweeps);
      ifpackSCList.set("relaxation: damping factor", SC_omega);
      ifpackSCType = "RELAXATION";
      ifpackSCList.set("relaxation: type", "Gauss-Seidel");
      smoProtoSC = rcp(new TrilinosSmoother(ifpackSCType, ifpackSCList, 0));
      smoProtoSC->SetFactory("A", SFact);
    } else {
      Teuchos::ParameterList ifpackDSList;
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
    MB->SetIgnoreUserData(true);               // always use data from factories defined in factory manager
    SimpleSm->SetSchurCompFactoryManager(MB);  // set temporary factory manager in BraessSarazin smoother

    // setup main factory manager
    RCP<FactoryManager> M = rcp(new FactoryManager());
    M->SetFactory("A", MueLu::NoFactory::getRCP());  // this is the 2x2 blocked operator
    M->SetFactory("Smoother", smootherFact);         // BraessSarazin block smoother
    M->SetFactory("PreSmoother", smootherFact);
    M->SetFactory("PostSmoother", smootherFact);

    MueLu::SetFactoryManager SFMCoarse(Finest, M);
    Finest->Request(MueLu::TopSmootherFactory<Scalar, LO, GO, Node>(M, "Smoother"));

    // call setup (= extract blocks and extract diagonal of F)
    SimpleSm->Setup(*Finest);

    RCP<MultiVector> xtest = MultiVectorFactory::Build(xstridedfullmap, 1);
    xtest->putScalar((Scalar)0.0);

    RCP<Vector> xR = Teuchos::rcp(new Xpetra::EpetraVectorT<int, Node>(epv));
    // calculate initial (absolute) residual
    Teuchos::Array<Teuchos::ScalarTraits<Scalar>::magnitudeType> norms(1);

    xR->norm2(norms);
    *out << "Test: ||x_0|| = " << norms[0] << std::endl;
    *out << "Test: Applying Simple Smoother" << std::endl;
    *out << "Test: START DATA" << std::endl;
    *out << "iterations\tVelocity_residual\tPressure_residual" << std::endl;
    SimpleSm->Apply(*xtest, *xR);
    xtest->norm2(norms);
    *out << "Test: ||x_1|| = " << norms[0] << std::endl;

    Teuchos::Array<Teuchos::ScalarTraits<Scalar>::magnitudeType> test = MueLu::Utilities<Scalar, LO, GO, Node>::ResidualNorm(*bOp, *xtest, *xR);
    *out << "residual norm: " << test[0] << std::endl;

    success = (test[0] < 1.0e-7);
    if (!success)
      *out << "no convergence" << std::endl;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return (success ? EXIT_SUCCESS : EXIT_FAILURE);
#else
  std::cout << "Epetra needs Serial node. Please recompile MueLu with the Serial node enabled." << std::endl;
  return EXIT_SUCCESS;
#endif  // #if defined(HAVE_MUELU_SERIAL) && defined(HAVE_MUELU_EPETRA)
}
