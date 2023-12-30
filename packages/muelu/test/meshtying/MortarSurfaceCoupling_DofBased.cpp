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

// Belos
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosStatusTestCombo.hpp>
#include <BelosXpetraStatusTestGenResSubNorm.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>

// MueLu
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_ConfigDefs.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

// Teuchos
#include <Teuchos_XMLParameterListHelpers.hpp>

// Xpetra
#include <Xpetra_BlockedMultiVector.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MapUtils.hpp>
#include <Xpetra_MatrixUtils.hpp>

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using ST  = Teuchos::ScalarTraits<Scalar>;
  using GST = Teuchos::ScalarTraits<GO>;

  RCP<const Teuchos::Comm<int>> comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out     = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  GO numGlobalDofsPrimal = -GST::one();
  clp.setOption("nPrimalDofs", &numGlobalDofsPrimal, "total number of primal DOFs");
  GO numGlobalDofsDual = -GST::one();
  clp.setOption("nDualDofs", &numGlobalDofsDual, "total number of dual DOFs");
  int numPrimalDofsPerNode = -1;
  clp.setOption("numPrimalDofsPerNode", &numPrimalDofsPerNode, "number of primal DOFs per mesh node");
  int numDualDofsPerNode = -1;
  clp.setOption("numDualDofsPerNode", &numDualDofsPerNode, "number of dual DOFs per interface node");
  std::string probName = "";
  clp.setOption("probName", &probName, "short name of the problem to be done. Used to read-in the problem from files.");
  std::string xmlFile = "";
  clp.setOption("xml", &xmlFile, "xml-file with MueLu configuration");
  int expectedNumIterations = -1;
  clp.setOption("expectedNumIts", &expectedNumIterations, "expected number of iterations for this test");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalDofsPrimal == -GST::one(), MueLu::Exceptions::InvalidArgument, "Please specify the global number of primal DOFs on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalDofsDual == -GST::one(), MueLu::Exceptions::InvalidArgument, "Please specify the global number of dual DOFs on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numPrimalDofsPerNode == -1, MueLu::Exceptions::InvalidArgument, "Please specify the number of primal DOFs per mesh node on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(numDualDofsPerNode == -1, MueLu::Exceptions::InvalidArgument, "Please specify the number of dual DOFs per interface node on the command line.");
  TEUCHOS_TEST_FOR_EXCEPTION(probName == "", MueLu::Exceptions::InvalidArgument, "Please specify a valid problem name.");
  TEUCHOS_TEST_FOR_EXCEPTION(xmlFile == "", MueLu::Exceptions::InvalidArgument, "Please specify a valid xml-file for the MueLu preconditioner.");

  const std::string matrixFileName           = probName + "_matrix.mm";
  const std::string rhsFileName              = probName + "_rhs.mm";
  const std::string nullspace1FileName       = probName + "_nullspace1.mm";
  const std::string dualInterfaceMapFileName = probName + "_interface_dof_map.mm";

  // Create maps for primal DOFs
  std::vector<size_t> stridingInfoPrimal;
  stridingInfoPrimal.push_back(numPrimalDofsPerNode);
  RCP<const StridedMap> dofRowMapPrimal = StridedMapFactory::Build(lib, numGlobalDofsPrimal, Teuchos::ScalarTraits<GO>::zero(), stridingInfoPrimal, comm, -1);

  // Create maps for dual DOFs
  std::vector<size_t> stridingInfoDual;
  stridingInfoDual.push_back(numDualDofsPerNode);
  RCP<const StridedMap> dofRowMapDual = StridedMapFactory::Build(lib, numGlobalDofsDual, Teuchos::ScalarTraits<GO>::zero(), stridingInfoDual, comm, -1, numGlobalDofsPrimal);

  // Construct the blocked map of the global system
  std::vector<RCP<const Map>> rowmaps;
  rowmaps.push_back(dofRowMapPrimal);
  rowmaps.push_back(dofRowMapDual);
  RCP<const Map> fullRowMap        = MapUtils::concatenateMaps(rowmaps);
  RCP<const BlockedMap> blockedMap = rcp(new BlockedMap(fullRowMap, rowmaps));

  // Read the matrix from file and transform it into a block matrix
  RCP<Matrix> mat = Xpetra::IO<SC, LO, GO, NO>::Read(matrixFileName, fullRowMap);
  RCP<MapExtractor> rangeMapExtractor =
      Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullRowMap, rowmaps);
  RCP<BlockedCrsMatrix> blockedMatrix =
      Xpetra::MatrixUtils<SC, LO, GO, NO>::SplitMatrix(*mat, rangeMapExtractor, rangeMapExtractor);
  blockedMatrix->fillComplete();

  // Read the right-hand side vector from file and transform it into a block vector
  RCP<MultiVector> rhs               = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(rhsFileName, fullRowMap);
  RCP<BlockedMultiVector> blockedRhs = rcp(new BlockedMultiVector(blockedMap, rhs));

  // Read the nullspace vector of the (0,0)-block from file
  RCP<MultiVector> nullspace1 = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector(nullspace1FileName, dofRowMapPrimal);

  // Read the interface slave dof row map from file
  TEUCHOS_ASSERT(!dualInterfaceMapFileName.empty());
  RCP<const Map> primalInterfaceDofMapDofs = Xpetra::IO<SC, LO, GO, NO>::ReadMap(dualInterfaceMapFileName, lib, comm);

  // Distribute interface slave dof row map across all procs
  RCP<const Map> primalInterfaceDofMap = Teuchos::null;
  {
    std::vector<int> interfaceGlobalDofsOnCurrentProc = std::vector<int>(primalInterfaceDofMapDofs->getGlobalNumElements());
    std::vector<int> interfaceGlobalDofs              = std::vector<int>(primalInterfaceDofMapDofs->getGlobalNumElements());
    if (comm->getRank() == 0) {
      for (size_t i = 0; i < primalInterfaceDofMapDofs->getGlobalNumElements(); ++i) {
        interfaceGlobalDofsOnCurrentProc[i] = primalInterfaceDofMapDofs->getLocalElementList()[i];
      }
    }
    Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_MAX, interfaceGlobalDofs.size(),
                            &interfaceGlobalDofsOnCurrentProc[0], &interfaceGlobalDofs[0]);

    Array<GlobalOrdinal> primalInterfaceDofMapDofsOnCurProc;
    for (size_t i = 0; i < interfaceGlobalDofs.size(); ++i) {
      if (dofRowMapPrimal->isNodeGlobalElement((GlobalOrdinal)interfaceGlobalDofs[i]))
        primalInterfaceDofMapDofsOnCurProc.push_back((GlobalOrdinal)interfaceGlobalDofs[i]);
    }
    primalInterfaceDofMap = MapFactory::Build(lib, primalInterfaceDofMapDofs->getGlobalNumElements(), primalInterfaceDofMapDofsOnCurProc, Teuchos::ScalarTraits<GO>::zero(), comm);
  }

  // Create the default nullspace vector of the (1,1)-block
  RCP<MultiVector> nullspace2 = MultiVectorFactory::Build(dofRowMapDual, 2, true);
  const int dimNS             = 2;
  for (int dim = 0; dim < dimNS; ++dim) {
    ArrayRCP<Scalar> nsValues = nullspace2->getDataNonConst(dim);
    const int numBlocks       = nsValues.size() / dimNS;
    for (int j = 0; j < numBlocks; ++j)
      nsValues[j * dimNS + dim] = Teuchos::ScalarTraits<Scalar>::one();
  }

  RCP<ParameterList> params = Teuchos::getParametersFromXmlFile(xmlFile);
  ParameterListInterpreter mueLuFactory(*params, comm);

  RCP<Hierarchy> hierarchy = mueLuFactory.CreateHierarchy();
  hierarchy->IsPreconditioner(true);
  hierarchy->SetDefaultVerbLevel(MueLu::Extreme);
  hierarchy->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(blockedMatrix));
  hierarchy->GetLevel(0)->Set("Nullspace1", nullspace1);
  hierarchy->GetLevel(0)->Set("Primal interface DOF map", primalInterfaceDofMap);

  mueLuFactory.SetupHierarchy(*hierarchy);

  // Create the preconditioned GMRES solver
  using OP = Belos::OperatorT<MultiVector>;

  using blockStatusTestClass = Belos::StatusTestGenResSubNorm<SC, MultiVector, OP>;
  using StatusTestComboClass = Belos::StatusTestCombo<SC, MultiVector, OP>;

  typename ST::magnitudeType tol  = 1e-4;
  typename ST::magnitudeType bTol = 1e-5;

  RCP<blockStatusTestClass> primalBlockStatusTest = rcp(new blockStatusTestClass(bTol, 0));
  RCP<blockStatusTestClass> dualBlockStatusTest   = rcp(new blockStatusTestClass(bTol, 1));

  RCP<StatusTestComboClass> statusTestCombo = rcp(new StatusTestComboClass(StatusTestComboClass::SEQ));
  statusTestCombo->addStatusTest(primalBlockStatusTest);
  statusTestCombo->addStatusTest(dualBlockStatusTest);

  RCP<ParameterList> belosParams = rcp(new ParameterList);
  belosParams->set("Flexible Gmres", false);
  belosParams->set("Num Blocks", 100);
  belosParams->set("Convergence Tolerance", tol);
  belosParams->set("Maximum Iterations", 100);
  belosParams->set("Verbosity", 33);
  belosParams->set("Output Style", 1);
  belosParams->set("Output Frequency", 1);

  using BLinProb    = Belos::LinearProblem<SC, MultiVector, OP>;
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(blockedMatrix));
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(hierarchy));

  RCP<MultiVector> X = MultiVectorFactory::Build(fullRowMap, 1, true);

  RCP<BLinProb> blinproblem = rcp(new BLinProb(belosOp, X, rhs));

  blinproblem->setRightPrec(belosPrec);
  blinproblem->setProblem();
  RCP<Belos::SolverManager<SC, MultiVector, OP>> blinsolver =
      rcp(new Belos::PseudoBlockGmresSolMgr<SC, MultiVector, OP>(blinproblem, belosParams));

  blinsolver->setUserConvStatusTest(statusTestCombo);

  Belos::ReturnType ret   = blinsolver->solve();
  const int numIterations = blinsolver->getNumIters();

  if (ret == Belos::Converged && (expectedNumIterations == -1 || numIterations == expectedNumIterations))
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

//-----------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
