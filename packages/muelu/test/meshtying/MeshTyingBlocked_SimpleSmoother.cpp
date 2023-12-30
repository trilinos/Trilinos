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

// MueLu
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_ConfigDefs.hpp>
#include <MueLu_ParameterListInterpreter.hpp>

// Teuchos
#include <Teuchos_XMLParameterListHelpers.hpp>

// Belos
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosStatusTestCombo.hpp>
#include <BelosXpetraStatusTestGenResSubNorm.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>

template <typename GlobalOrdinal>
void read_Lagr2Dof(std::string filemane, std::map<GlobalOrdinal, GlobalOrdinal> &lagr2Dof) {
  std::fstream lagr2DofFile;
  lagr2DofFile.open(filemane);
  TEUCHOS_ASSERT(lagr2DofFile.is_open())

  GlobalOrdinal key;
  GlobalOrdinal value;
  while (lagr2DofFile >> key >> value) {
    lagr2Dof[key] = value;
  }
  lagr2DofFile.close();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib &lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using SparseMatrixType    = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using tpetra_mvector_type = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using tpetra_map_type     = Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>;

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace Teuchos;

  using ST = ScalarTraits<Scalar>;

  oblackholestream blackhole;

  RCP<const Comm<int>> comm = DefaultComm<int>::getComm();
  RCP<FancyOStream> out     = fancyOStream(rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  const GO globalPrimalNumDofs    = 1599;
  const GO globalDualNumDofs      = 100;
  const GO globalNumDofs          = globalPrimalNumDofs + globalDualNumDofs;  // used for the maps
  const size_t nPrimalDofsPerNode = 3;
  const GO globalPrimalNumNodes   = globalPrimalNumDofs / nPrimalDofsPerNode;
  const size_t nDualDofsPerNode   = 1;

  std::map<GO, GO> lagr2Dof;
  std::map<LO, LO> myLagr2Dof;
  read_Lagr2Dof<GO>("Lagr2Dof.txt", lagr2Dof);

  // Construct the blocked map in Thyra mode
  RCP<const tpetra_map_type> primalNodeMap = Tpetra::createUniformContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(globalPrimalNumNodes, comm);
  const GO indexBase                       = primalNodeMap->getIndexBase();
  ArrayView<const GO> myPrimalNodes        = primalNodeMap->getLocalElementList();

  const size_t numMyPrimalNodes = primalNodeMap->getLocalNumElements();
  const size_t numMyPrimalDofs  = numMyPrimalNodes * nPrimalDofsPerNode;

  Array<GO> myPrimalDofs(numMyPrimalDofs);

  LO current_i = 0;
  for (size_t i = 0; i < numMyPrimalNodes; ++i)
    for (size_t j = 0; j < nPrimalDofsPerNode; ++j)
      myPrimalDofs[current_i++] = myPrimalNodes[i] * nPrimalDofsPerNode + j;

  RCP<const tpetra_map_type> primalMap = rcp(new tpetra_map_type(globalPrimalNumDofs, myPrimalDofs, indexBase, comm));

  size_t numMyDualDofs = 0;

  for (auto i = lagr2Dof.begin(); i != lagr2Dof.end(); ++i)
    if (primalMap->isNodeGlobalElement(nPrimalDofsPerNode * (i->second)))
      ++numMyDualDofs;

  numMyDualDofs *= nDualDofsPerNode;

  const size_t numMyDofs = numMyPrimalDofs + numMyDualDofs;

  Array<GO> myDualDofs(numMyDualDofs);
  Array<GO> myDofs(numMyDofs);

  for (size_t i = 0; i < numMyPrimalDofs; ++i)
    myDofs[i] = myPrimalDofs[i];

  current_i = 0;
  for (auto i = lagr2Dof.begin(); i != lagr2Dof.end(); ++i)
    if (primalMap->isNodeGlobalElement(nPrimalDofsPerNode * (i->second))) {
      for (size_t j = 0; j < nDualDofsPerNode; ++j) {
        myDualDofs[nDualDofsPerNode * current_i + j]               = (i->first) * nDualDofsPerNode + j;
        myDofs[numMyPrimalDofs + nDualDofsPerNode * current_i + j] = globalPrimalNumDofs + (i->first) * nDualDofsPerNode + j;
      }
      GO primalDof          = nPrimalDofsPerNode * (i->second);
      myLagr2Dof[current_i] = primalMap->getLocalElement(primalDof) / nPrimalDofsPerNode;
      ++current_i;
    }

  RCP<const tpetra_map_type> dualMap = rcp(new tpetra_map_type(globalDualNumDofs, myDualDofs, indexBase, comm));
  RCP<const tpetra_map_type> fullMap = rcp(new tpetra_map_type(globalNumDofs, myDofs, indexBase, comm));

  RCP<const Map> fullXMap   = rcp(new TpetraMap(fullMap));
  RCP<const Map> primalXMap = rcp(new TpetraMap(primalMap));
  RCP<const Map> dualXMap   = rcp(new TpetraMap(dualMap));

  std::vector<RCP<const Map>> xsubmaps = {primalXMap, dualXMap};
  RCP<BlockedMap> blockedMap           = rcp(new BlockedMap(fullXMap, xsubmaps, true));

  // Read input matrices
  typedef Tpetra::MatrixMarket::Reader<SparseMatrixType> reader_type;

  RCP<Matrix> xQ  = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("Q_mm.txt", primalXMap, null, primalXMap, primalXMap);
  RCP<Matrix> xG  = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("G_mm.txt", primalXMap, null, dualXMap, primalXMap);
  RCP<Matrix> xGT = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("GT_mm.txt", dualXMap, null, primalXMap, dualXMap);
  RCP<Matrix> xC  = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("C_mm.txt", dualXMap, null, dualXMap, dualXMap);

  RCP<CrsMatrixWrap> xwQ  = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xQ);
  RCP<CrsMatrixWrap> xwG  = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xG);
  RCP<CrsMatrixWrap> xwGT = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xGT);
  RCP<CrsMatrixWrap> xwC  = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(xC);

  // Construct the blocked saddle-point matrix
  RCP<BlockedCrsMatrix> blockedMatrix = rcp(new BlockedCrsMatrix(blockedMap, blockedMap, 8));

  blockedMatrix->setMatrix(0, 0, xwQ);
  blockedMatrix->setMatrix(0, 1, xwG);
  blockedMatrix->setMatrix(1, 0, xwGT);
  blockedMatrix->setMatrix(1, 1, xwC);
  blockedMatrix->fillComplete();

  // Create the preconditioner
  std::string xmlFile = "simple_1dof.xml";

  RCP<ParameterList> params     = Teuchos::getParametersFromXmlFile(xmlFile);
  ParameterList &userDataParams = params->sublist("user data");
  userDataParams.set<RCP<std::map<LO, LO>>>("DualNodeID2PrimalNodeID", rcpFromRef(myLagr2Dof));

  RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner(Teuchos::rcp_dynamic_cast<Matrix>(blockedMatrix, true), *params);
  H->IsPreconditioner(true);

  // Create the preconditioned GMRES solver
  typedef typename tpetra_mvector_type::dot_type belos_scalar;
  typedef Belos::OperatorT<MultiVector> OP;

  typedef Belos::StatusTestGenResSubNorm<belos_scalar, MultiVector, OP> blockStatusTestClass;
  typedef Belos::StatusTestCombo<belos_scalar, MultiVector, OP> StatusTestComboClass;

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
  belosParams->set("Convergence Tolerance", belos_scalar(tol));
  belosParams->set("Maximum Iterations", 100);
  belosParams->set("Verbosity", 33);
  belosParams->set("Output Style", 1);
  belosParams->set("Output Frequency", 1);

  typedef Belos::LinearProblem<belos_scalar, MultiVector, OP> BLinProb;
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(blockedMatrix));
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(H));

  RCP<tpetra_mvector_type> rhsMultiVector      = reader_type::readDenseFile("f_mm.txt", comm, fullMap);
  RCP<tpetra_mvector_type> solutionMultiVector = rcp(new tpetra_mvector_type(fullMap, 1));

  RCP<MultiVector> rhsXMultiVector      = rcp(new TpetraMultiVector(rhsMultiVector));
  RCP<MultiVector> solutionXMultiVector = rcp(new TpetraMultiVector(solutionMultiVector));

  RCP<BLinProb> blinproblem = rcp(new BLinProb(belosOp, solutionXMultiVector, rhsXMultiVector));

  blinproblem->setRightPrec(belosPrec);
  blinproblem->setProblem();
  RCP<Belos::SolverManager<belos_scalar, MultiVector, OP>> blinsolver =
      rcp(new Belos::PseudoBlockGmresSolMgr<belos_scalar, MultiVector, OP>(blinproblem, belosParams));

  blinsolver->setUserConvStatusTest(statusTestCombo);

  Belos::ReturnType ret = blinsolver->solve();

  if (ret == Belos::Converged)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
