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
void read_Lagr2Dof(std::string filemane, std::map<GlobalOrdinal, GlobalOrdinal> &lagr2Dof)
{
  std::fstream lagr2DofFile;
  lagr2DofFile.open(filemane);
  GlobalOrdinal key;
  GlobalOrdinal value;
  while (lagr2DofFile >> key >> value)
  {
    lagr2Dof[key] = value;
  }
  lagr2DofFile.close();
}

int main(int argc, char *argv[])
{
#ifdef HAVE_MUELU_TPETRA
  typedef double Scalar;
  typedef MueLu::DefaultLocalOrdinal LocalOrdinal;
  typedef MueLu::DefaultGlobalOrdinal GlobalOrdinal;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef MueLu::DefaultNode Node;

  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> SparseMatrixType;

  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> tpetra_mvector_type;

  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> tpetra_map_type;
  typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> xpetra_map_type;
  typedef Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> xpetra_tmap_type;
  typedef Xpetra::BlockedMap<LocalOrdinal, GlobalOrdinal, Node> xpetra_bmap_type;

  typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_crs_type;
  typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_matrix_type;
  typedef Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_wcrs_type;
  typedef Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_bcrs_type;
  typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_tcrs_type;

  typedef Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_mvector_type;
  typedef Xpetra::TpetraMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> xpetra_tmvector_type;

#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace Teuchos;

  typedef ScalarTraits<Scalar> ST;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc, &argv, &blackhole);

  RCP<const Comm<int>> comm = DefaultComm<int>::getComm();
  RCP<FancyOStream> out = fancyOStream(rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  GO globalPrimalNumDofs = 1599;
  GO globalDualNumDofs = 100;
  GO globalNumDofs = globalPrimalNumDofs + globalDualNumDofs; // used for the maps
  size_t nPrimalDofsPerNode = 3;
  const GO globalPrimalNumNodes = globalPrimalNumDofs / nPrimalDofsPerNode;
  size_t nDualDofsPerNode = 1;

  std::map<GO, GO> lagr2Dof;
  std::map<LO, LO> myLagr2Dof;
  read_Lagr2Dof<GO>("Lagr2Dof.txt", lagr2Dof);

  // Construct the blocked map in Thyra mode
  RCP<const tpetra_map_type> primalNodeMap = Tpetra::createUniformContigMapWithNode<LocalOrdinal, GlobalOrdinal, Node>(globalPrimalNumNodes, comm);
  GO indexBase = primalNodeMap->getIndexBase();
  ArrayView<const GO> myPrimalNodes = primalNodeMap->getNodeElementList();

  const size_t numMyPrimalNodes = primalNodeMap->getNodeNumElements();
  const size_t numMyPrimalDofs = numMyPrimalNodes * nPrimalDofsPerNode;

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
    if (primalMap->isNodeGlobalElement(nPrimalDofsPerNode * (i->second)))
    {
      for (size_t j = 0; j < nDualDofsPerNode; ++j)
      {
        myDualDofs[nDualDofsPerNode * current_i + j] = (i->first) * nDualDofsPerNode + j;
        myDofs[numMyPrimalDofs + nDualDofsPerNode * current_i + j] = globalPrimalNumDofs + (i->first) * nDualDofsPerNode + j;
      }
      GO primalDof = nPrimalDofsPerNode * (i->second);
      myLagr2Dof[current_i] = primalMap->getLocalElement(primalDof) / nPrimalDofsPerNode;
      ++current_i;
    }

  RCP<const tpetra_map_type> dualMap = rcp(new tpetra_map_type(globalDualNumDofs, myDualDofs, indexBase, comm));
  RCP<const tpetra_map_type> fullMap = rcp(new tpetra_map_type(globalNumDofs, myDofs, indexBase, comm));

  RCP<const xpetra_map_type> fullXMap = rcp(new xpetra_tmap_type(fullMap));
  RCP<const xpetra_map_type> primalXMap = rcp(new xpetra_tmap_type(primalMap));
  RCP<const xpetra_map_type> dualXMap = rcp(new xpetra_tmap_type(dualMap));

  std::vector<RCP<const xpetra_map_type>> xsubmaps = {primalXMap, dualXMap};
  RCP<xpetra_bmap_type> blockedMap = rcp(new xpetra_bmap_type(fullXMap, xsubmaps, true));

  // Read input matrices
  typedef Tpetra::MatrixMarket::Reader<SparseMatrixType> reader_type;

  RCP<xpetra_matrix_type> xQ = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("Q_mm.txt", primalXMap, null, primalXMap, primalXMap);
  RCP<xpetra_matrix_type> xG = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("G_mm.txt", primalXMap, null, dualXMap, primalXMap);
  RCP<xpetra_matrix_type> xGT = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("GT_mm.txt", dualXMap, null, primalXMap, dualXMap);
  RCP<xpetra_matrix_type> xC = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Read("C_mm.txt", dualXMap, null, dualXMap, dualXMap);

  RCP<xpetra_wcrs_type> xwQ = Teuchos::rcp_dynamic_cast<xpetra_wcrs_type>(xQ);
  RCP<xpetra_wcrs_type> xwG = Teuchos::rcp_dynamic_cast<xpetra_wcrs_type>(xG);
  RCP<xpetra_wcrs_type> xwGT = Teuchos::rcp_dynamic_cast<xpetra_wcrs_type>(xGT);
  RCP<xpetra_wcrs_type> xwC = Teuchos::rcp_dynamic_cast<xpetra_wcrs_type>(xC);

  // Construct the blocked saddle-point matrix
  RCP<xpetra_bcrs_type> blockedMatrix = rcp(new xpetra_bcrs_type(blockedMap, blockedMap, 8));

  blockedMatrix->setMatrix(0, 0, xwQ);
  blockedMatrix->setMatrix(0, 1, xwG);
  blockedMatrix->setMatrix(1, 0, xwGT);
  blockedMatrix->setMatrix(1, 1, xwC);
  blockedMatrix->fillComplete();

  // Create the preconditioner
  std::string xmlFile = "myXML.xml";

  RCP<ParameterList> params = Teuchos::getParametersFromXmlFile(xmlFile);
  ParameterList& userDataParams = params->sublist("user data");
  userDataParams.set<RCP<std::map<LO, LO>>>("DualNodeID2PrimalNodeID", rcpFromRef(myLagr2Dof));

  RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner(Teuchos::rcp_dynamic_cast<Matrix>(blockedMatrix, true), *params);
  H->IsPreconditioner(true);

  // Create the preconditioned GMRES solver
  typedef xpetra_mvector_type MV;
  typedef typename tpetra_mvector_type::dot_type belos_scalar;
  typedef Belos::OperatorT<MV> OP;

  typedef Belos::StatusTestGenResSubNorm<belos_scalar, MV, OP> blockStatusTestClass;
  typedef Belos::StatusTestCombo<belos_scalar, MV, OP> StatusTestComboClass;

  typename ST::magnitudeType tol = 1e-4;
  typename ST::magnitudeType bTol = 1e-5;

  RCP<blockStatusTestClass> primalBlockStatusTest = rcp(new blockStatusTestClass(bTol, 0));
  RCP<blockStatusTestClass> dualBlockStatusTest = rcp(new blockStatusTestClass(bTol, 1));

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

  typedef Belos::LinearProblem<belos_scalar, MV, OP> BLinProb;
  RCP<OP> belosOp = rcp(new Belos::XpetraOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(blockedMatrix));
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(H));

  RCP<tpetra_mvector_type> rhsMultiVector = reader_type::readDenseFile("f_mm.txt", comm, fullMap);
  RCP<tpetra_mvector_type> solutionMultiVector = rcp(new tpetra_mvector_type(fullMap, 1));

  RCP<xpetra_mvector_type> rhsXMultiVector = rcp(new xpetra_tmvector_type(rhsMultiVector));
  RCP<xpetra_mvector_type> solutionXMultiVector = rcp(new xpetra_tmvector_type(solutionMultiVector));

  RCP<BLinProb> blinproblem = rcp(new BLinProb(belosOp, solutionXMultiVector, rhsXMultiVector));

  blinproblem->setRightPrec(belosPrec);
  blinproblem->setProblem();
  RCP<Belos::SolverManager<belos_scalar, MV, OP>> blinsolver =
      rcp(new Belos::PseudoBlockGmresSolMgr<belos_scalar, MV, OP>(blinproblem, belosParams));

  blinsolver->setUserConvStatusTest(statusTestCombo);

  Belos::ReturnType ret = blinsolver->solve();

  if (ret == Belos::Converged)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;

#else
  std::cout << "Tpetra not enabled. Skip test." << std::endl;
  return EXIT_SUCCESS;
#endif
}
