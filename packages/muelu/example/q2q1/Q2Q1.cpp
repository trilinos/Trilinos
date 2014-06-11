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
#include <cstdio>
#include <unistd.h>
#include <iostream>

#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_EasyParameterListInterpreter.hpp>

#include <MueLu_BlockedDirectSolver.hpp>
#include <MueLu_BlockedPFactory.hpp>
#include <MueLu_BlockedRAPFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_ConstraintFactory.hpp>
#include <MueLu_EminPFactory.hpp>
#include <MueLu_FactoryManager.hpp>
#include <MueLu_FilteredAFactory.hpp>
#include <MueLu_GenericRFactory.hpp>
#include <MueLu_PatternFactory.hpp>
#include <MueLu_Q2Q1PFactory.hpp>
#include <MueLu_Q2Q1uPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_SubBlockAFactory.hpp>

#include <MueLu_Utilities.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include "MueLu_SmootherFactory.hpp"
#include <MueLu_Ifpack2Smoother.hpp>
#include "MueLu_TrilinosSmoother.hpp"

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif


namespace MueLuTests {

#include <MueLu_UseShortNames.hpp>
  Teuchos::RCP<Matrix> FilterMatrix(Matrix& A, SC dropTol);

  void SetDependencyTree(FactoryManager& M, const std::string& mode1, const std::string& mode2);

  Teuchos::RCP<MultiVector> BuildCoords(Teuchos::RCP<const Map> map, int NDim);

}

#define PRESSURE_FIRST

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using Teuchos::Array;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::null;

  using namespace MueLuTests;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  TEUCHOS_TEST_FOR_EXCEPTION(comm->getSize() > 1, MueLu::Exceptions::RuntimeError,
                             "For now, Q2Q1 only works in the serial mode. "
                             "We are working on the parallel implementation.");

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  Teuchos::CommandLineProcessor clp(false);

  Xpetra::Parameters xpetraParameters(clp);

  std::string xmlFileName  = "driver.xml";   clp.setOption("xml",      &xmlFileName,   "read parameters from a file [default = 'driver.xml']");
  double      tol          = 1e-12;          clp.setOption("tol",      &tol,           "solver convergence tolerance");
  int         n            = 65;              clp.setOption("n",        &n,             "problem size (1D)");
  int         maxLevels    = 3;              clp.setOption("nlevels",  &maxLevels,     "max num levels");
  int         compare      = 0;              clp.setOption("compare",  &compare,       "compare block and point hierarchies");


  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);

  // =========================================================================
  // Problem construction
  // =========================================================================
  // Step 1: construct maps
  Xpetra::global_size_t numVelElements  = 2*(2*n-1)*(2*n-1);
  Xpetra::global_size_t numPresElements = n*n;
  Xpetra::global_size_t numElements     = numVelElements + numPresElements;

  const GO indexBase = 0;
  std::vector<size_t> stridingInfo(1,1);
  int stridedBlockId = -1;

  Array<GO> elementList(numElements);
  for (LO i = 0; i < numElements; i++)
    elementList[i] = i;
  RCP<Map> fullMap = StridedMapFactory::Build(lib, numElements, elementList(), indexBase, stridingInfo, comm);

  std::vector<RCP<const Map> > partMaps(2);
  partMaps[0] = StridedMapFactory::Build(lib, numVelElements,  elementList(0,              numVelElements),  indexBase,
                                         stridingInfo, comm);
  partMaps[1] = StridedMapFactory::Build(lib, numPresElements, elementList(numVelElements, numPresElements), indexBase,
                                         stridingInfo, comm, stridedBlockId, numVelElements);
#ifdef PRESSURE_FIRST
  std::vector<RCP<const Map> > partMaps1 = partMaps;
  std::swap(partMaps1[0], partMaps1[1]);
  RCP<const MapExtractor> mapExtractor = MapExtractorFactory::Build(fullMap, partMaps1);
#else
  RCP<const MapExtractor> mapExtractor = MapExtractorFactory::Build(fullMap, partMaps);
#endif

  // Step 2: read in matrices
  std::string matrixPrefix = "Q2Q1_" + MueLu::toString(n) + "x" + MueLu::toString(n);

  const bool binaryFormat = true;
  RCP<Matrix>    A_11     = Utils ::Read(matrixPrefix + "_A.dat", lib, comm, binaryFormat);
  RCP<Matrix>    A_21     = Utils ::Read(matrixPrefix + "_B.dat", partMaps[1], partMaps[0], partMaps[0], partMaps[1], true, binaryFormat);
  RCP<Matrix>    A_12     = Utils2::Transpose(*A_21, true);
  RCP<Matrix>    A_22     = Teuchos::null;

  RCP<CrsMatrix> A_11_crs = rcp_dynamic_cast<CrsMatrixWrap>(A_11)->getCrsMatrix();
  RCP<CrsMatrix> A_12_crs = rcp_dynamic_cast<CrsMatrixWrap>(A_12)->getCrsMatrix();
  RCP<CrsMatrix> A_21_crs = rcp_dynamic_cast<CrsMatrixWrap>(A_21)->getCrsMatrix();
  RCP<CrsMatrix> A_22_crs = Teuchos::null;

  RCP<BlockedCrsMatrix> A = rcp(new BlockedCrsMatrix(mapExtractor, mapExtractor, 10));
#ifdef PRESSURE_FIRST
  A->setMatrix(0, 0, A_22_crs);
  A->setMatrix(0, 1, A_21_crs);
  A->setMatrix(1, 0, A_12_crs);
  A->setMatrix(1, 1, A_11_crs);
#else
  A->setMatrix(0, 0, A_11_crs);
  A->setMatrix(0, 1, A_12_crs);
  A->setMatrix(1, 0, A_21_crs);
  A->setMatrix(1, 1, A_22_crs);
#endif
  A->fillComplete();
  A->Merge();           // rst: do we need this merge? 

  // Step 3: construct coordinates
  const int NDim = 2;
  RCP<MultiVector> coords1 = BuildCoords(partMaps[0], NDim);
  RCP<MultiVector> coords2 = BuildCoords(partMaps[1], NDim);

  // =========================================================================
  // Preconditioner construction - I (block)
  // =========================================================================
  RCP<Matrix> BBt = Utils::Multiply(*A_21, false, *A_12, false, out);

  // Filter matrices
  SC dropTol = 0.06;
  RCP<Matrix> filteredA = FilterMatrix(*A_11, dropTol);
  RCP<Matrix> filteredB = FilterMatrix(*BBt,  dropTol);

  RCP<CrsMatrix> fA_11_crs = rcp_dynamic_cast<CrsMatrixWrap>(filteredA)->getCrsMatrix();
  RCP<CrsMatrix> fA_12_crs = Teuchos::null;
  RCP<CrsMatrix> fA_21_crs = Teuchos::null;
  RCP<CrsMatrix> fA_22_crs = rcp_dynamic_cast<CrsMatrixWrap>(filteredB)->getCrsMatrix();

  RCP<BlockedCrsMatrix> fA = rcp(new BlockedCrsMatrix(mapExtractor, mapExtractor, 10));
#ifdef PRESSURE_FIRST
  fA->setMatrix(0, 0, fA_22_crs);
  fA->setMatrix(0, 1, fA_12_crs);
  fA->setMatrix(1, 0, fA_21_crs);
  fA->setMatrix(1, 1, fA_11_crs);
#else
  fA->setMatrix(0, 0, fA_11_crs);
  fA->setMatrix(0, 1, fA_12_crs);
  fA->setMatrix(1, 0, fA_21_crs);
  fA->setMatrix(1, 1, fA_22_crs);
#endif
  fA->fillComplete();

  // Construct the hierarchy
  FactoryManager M;
#ifdef PRESSURE_FIRST
  SetDependencyTree(M, "pressure", "velocity");
#else
  SetDependencyTree(M, "velocity", "pressure");
#endif

  std::vector<RCP<Hierarchy> > H(compare+1);
  H[0] = rcp(new Hierarchy);
  H[0]->GetLevel(0)->Set("A", rcp_dynamic_cast<Matrix>(fA));
  H[0]->GetLevel(0)->Set("CoordinatesVelocity", coords1);
  H[0]->GetLevel(0)->Set("CoordinatesPressure", coords2);
  H[0]->SetMaxCoarseSize(1);

  // The 1st invocation of Setup() builds the hierarchy using
  // the filtered matrix. This build includes the grid transfers
  // but not the creation of the smoothers. Note: we need to indicate
  // what should be kept from the 1st invocation for the 2nd
  // invocation, which then focuses on building the smoothers
  // for the unfiltered matrix.

  H[0]->Keep("P", M.GetFactory("P").get());
  H[0]->Keep("R", M.GetFactory("R").get());
  H[0]->Keep("Ptent", M.GetFactory("Ptent").get());
  H[0]->Setup(M, 0, maxLevels);

  // =========================================================================
  // Preconditioner construction - II (point)
  // =========================================================================
  if (compare) {
    ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);

    RCP<HierarchyManager> mueLuFactory = rcp(new EasyParameterListInterpreter(paramList));
    H[1] = mueLuFactory->CreateHierarchy();
    H[1]->GetLevel(0)->Set("A", rcp_dynamic_cast<Matrix>(rcp(new CrsMatrixWrap(A->Merge()))));
    mueLuFactory->SetupHierarchy(*H[1]);
  }

  // Set up Vanka smoothing via a combination of Schwarz and
  // block relaxation. 

  Teuchos::ParameterList OuterList;
  OuterList.set("schwarz: overlap level", (int) 0);
  OuterList.set("schwarz: zero starting solution", false);
  OuterList.set("subdomain solver name", "Block_Relaxation");
  Teuchos::ParameterList& InnerList = OuterList.sublist("subdomain solver parameters");
  InnerList.set("partitioner: type","user");
  InnerList.set("partitioner: overlap"    ,(int) 1);
  InnerList.set("relaxation: type"        ,"Gauss-Seidel");
  InnerList.set("relaxation: sweeps"      ,(int) 1);
  InnerList.set("relaxation: damping factor",0.5);
//InnerList.set("relaxation: backward mode",true);  NOT SUPPORTED YET
  InnerList.set("relaxation: zero starting solution", false);
  std::string ifpackType = "SCHWARZ";
  RCP<SmootherPrototype> smootherPrototype = rcp(new TrilinosSmoother(ifpackType, OuterList));
  M.SetFactory("Smoother", rcp(new SmootherFactory(smootherPrototype)));
  M.SetFactory("CoarseSolver", rcp(new SmootherFactory(smootherPrototype)));

  // For the smoother setup use the unfiltered matrix

  H[0]->GetLevel(0)->Set("A", rcp_dynamic_cast<Matrix>(A));
  H[0]->Setup(M, 0, H[0]->GetNumLevels());

  // =========================================================================
  // System solution (Ax = b) - I (block)
  // =========================================================================
  RCP<Vector> X = VectorFactory::Build(fullMap);
  RCP<Vector> B = VectorFactory::Build(fullMap);
  {
    // we set seed for reproducibility
    Utils::SetRandomSeed(*comm);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

    Teuchos::Array<STS::magnitudeType> norms(1);
    B->norm2(norms);
    B->scale(one/norms[0]);
    X->putScalar(zero);
  }

  // Operator and Multivector type that will be used with Belos
  typedef MultiVector          MV;
  typedef Belos::OperatorT<MV> OP;

  // Define Belos Operator
  Teuchos::RCP<OP> belosOp = rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(A)); // Turns a Xpetra::Matrix object into a Belos operator

  // Belos parameter list
  int maxIts = 2000;
  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations",    20);
  belosList.set("Convergence Tolerance", 1e-12);
  belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  belosList.set("Output Frequency",      1);
  belosList.set("Output Style",          Belos::Brief);

  for (int i = 0; i <= compare; i++) {
    H[i]->IsPreconditioner(true);

    // Define Belos Preconditioner
    Teuchos::RCP<OP> belosPrec = rcp(new Belos::MueLuOp <SC, LO, GO, NO, LMO>(H[i])); // Turns a MueLu::Hierarchy object into a Belos operator

    // Construct a Belos LinearProblem object
    RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
    belosProblem->setRightPrec(belosPrec);

    bool set = belosProblem->setProblem();
    if (set == false) {
      out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // Create an iterative solver manager
    // We use GMRES because it is a saddle point problem
    RCP< Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

    // Perform solve
    Belos::ReturnType ret = Belos::Unconverged;
    try {
      ret = solver->solve();

      // Get the number of iterations for this solve.
      out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;

    } catch (...) {
      out << "\nERROR:  Belos threw an error! " << std::endl;
    }

    // Check convergence
    if (ret != Belos::Converged)
      out << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    else
      out << std::endl << "SUCCESS:  Belos converged!" << std::endl;
  }

  return 0;
}

namespace MueLuTests {

#include <MueLu_UseShortNames.hpp>
  Teuchos::RCP<Matrix> FilterMatrix(Matrix& A, SC dropTol) {
    using Teuchos::RCP;

    Level level;
    level.SetLevelID(1);
    level.Set<RCP<Matrix> >("A", rcpFromRef(A));

    FactoryManager M;
    level.SetFactoryManager(rcpFromRef(M));

    RCP<CoalesceDropFactory> dropFactory = rcp(new CoalesceDropFactory());
    ParameterList dropParams = *(dropFactory->GetValidParameterList());
    dropParams.set("lightweight wrap",      true);
    dropParams.set("algorithm",             "classical");
    dropParams.set("aggregation threshold", dropTol);
    // dropParams.set("Dirichlet detection threshold", <>);
    dropFactory->SetParameterList(dropParams);
    M.SetFactory("Graph",     dropFactory);
    M.SetFactory("Filtering", dropFactory);

    RCP<FilteredAFactory> filterFactory = rcp(new FilteredAFactory());
    ParameterList filterParams = *(filterFactory->GetValidParameterList());
    filterParams.set("filtered matrix: reuse graph", false);
    filterFactory->SetParameterList(filterParams);
    filterFactory->SetFactory("Graph", dropFactory);

    // Build
    level.Request("A", filterFactory.get());
    filterFactory->Build(level);

    RCP<Matrix> filteredA;
    level.Get("A", filteredA, filterFactory.get());

    return filteredA;
  }

  void SetBlockDependencyTree(FactoryManager& M, int row, int col, const std::string&);

  void SetDependencyTree(FactoryManager& M, const std::string& mode1, const std::string& mode2) {
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<FactoryManager> M11 = rcp(new FactoryManager);
    SetBlockDependencyTree(*M11, 0, 0, mode1);

    RCP<FactoryManager> M22 = rcp(new FactoryManager);
    SetBlockDependencyTree(*M22, 1, 1, mode2);

    RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());
    PFact->AddFactoryManager(M11);
    PFact->AddFactoryManager(M22);
    M.SetFactory("P", PFact);

    RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
    RFact->SetFactory("P", PFact);
    M.SetFactory("R", RFact);

    RCP<Factory> AcFact = rcp(new BlockedRAPFactory());
    AcFact->SetFactory("P", PFact);
    AcFact->SetFactory("R", RFact);
    M.SetFactory("A", AcFact);

    M.SetFactory("Smoother",     Teuchos::null);
    M.SetFactory("CoarseSolver", Teuchos::null);

    RCP<Factory> coarseFact = rcp(new SmootherFactory(rcp(new BlockedDirectSolver()), Teuchos::null));
    // M.SetFactory("CoarseSolver", coarseFact);
    M.SetFactory("CoarseSolver", Teuchos::null);
  }

  void SetBlockDependencyTree(FactoryManager& M, int row, int col, const std::string& mode) {
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef MueLu::Q2Q1PFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Q2Q1PFactory;
    typedef MueLu::Q2Q1uPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> Q2Q1uPFactory;

    RCP<SubBlockAFactory> AFact = rcp(new SubBlockAFactory());
    AFact->SetFactory  ("A",         MueLu::NoFactory::getRCP());
    AFact->SetParameter("block row", Teuchos::ParameterEntry(row));
    AFact->SetParameter("block col", Teuchos::ParameterEntry(col));
    M.SetFactory("A", AFact);

    // RCP<Q2Q1PFactory>  Q2Q1Fact = rcp(new Q2Q1PFactory);
    RCP<Q2Q1uPFactory> Q2Q1Fact = rcp(new Q2Q1uPFactory);
    ParameterList q2q1ParamList = *(Q2Q1Fact->GetValidParameterList());
    q2q1ParamList.set("mode", mode);
    // q2q1ParamList.set("phase2", false);
    Q2Q1Fact->SetParameterList(q2q1ParamList);
    Q2Q1Fact->SetFactory("A", AFact);
    M.SetFactory("Ptent", Q2Q1Fact);

    RCP<PatternFactory> patternFact = rcp(new PatternFactory);
    ParameterList patternParams = *(patternFact->GetValidParameterList());
    patternParams.set("k", 0);
    patternFact->SetParameterList(patternParams);
    patternFact->SetFactory("A", AFact);
    patternFact->SetFactory("P", Q2Q1Fact);
    M.SetFactory("Ppattern", patternFact);

    RCP<ConstraintFactory> CFact = rcp(new ConstraintFactory);
    CFact->SetFactory("Ppattern", patternFact);
    M.SetFactory("Constraint", CFact);

    RCP<EminPFactory> EminPFact = rcp(new EminPFactory());
    EminPFact->SetFactory("A",          AFact);
    EminPFact->SetFactory("Constraint", CFact);
    EminPFact->SetFactory("P",          Q2Q1Fact);
    M.SetFactory("P", EminPFact);
  }

  Teuchos::RCP<MultiVector> BuildCoords(Teuchos::RCP<const Map> map, int NDim) {
    Teuchos::RCP<MultiVector> coords = MultiVectorFactory::Build(map, NDim);

    TEUCHOS_TEST_FOR_EXCEPTION(NDim != 2, MueLu::Exceptions::RuntimeError, "Need dimension 2");

    Teuchos::ArrayRCP<Teuchos::ArrayRCP<SC> > coord1D(NDim);
    coord1D[0] = coords->getDataNonConst(0);
    coord1D[1] = coords->getDataNonConst(1);

    double hx = 1.0;
    double hy = 1.0;

    int N = map->getGlobalNumElements();
    int n = sqrt(N);
    std::cout << "N = " << N << ", n = " << n << std::endl;
    if (N == n*n) {
      // pressure coords
      for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
          coord1D[0][j*n+i] = i*hx;
          coord1D[1][j*n+i] = j*hy;
        }

    } else {
      // velocity coords
      n = sqrt(N/2);

      hx *= 0.5;
      hy *= 0.5;

      for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
          coord1D[0][2*(j*n+i)+0] = coord1D[0][2*(j*n+i)+1] = i*hx;
          coord1D[1][2*(j*n+i)+0] = coord1D[1][2*(j*n+i)+1] = j*hy;
        }
    }
    return coords;
  }

}
