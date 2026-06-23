// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * Author: Zhen Wang
 * Email: wangz@ornl.gov
 *        zhen.wang@alum.emory.edu
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>

// Teuchos
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp_decl.hpp"

// Tpetra
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Import.hpp"
#include "MatrixMarket_Tpetra.hpp"

// Teko
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_StridedTpetraOperator.hpp"
#include "Teko_TpetraBlockPreconditioner.hpp"
#include "Teko_TpetraBlockedMappingStrategy.hpp"
#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#include "Teko_ConfigDefs.hpp"

// AL / ModAL
#include "Teko_ALOperator.hpp"
#include "Teko_InvModALStrategy.hpp"
#include "Teko_ModALPreconditionerFactory.hpp"

// Belos
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosTpetraAdapter.hpp"

using Teuchos::null;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using namespace Thyra;
using namespace Teko;
using namespace std;

TEUCHOS_UNIT_TEST(tModALPreconditioner, test_tpetra) {
  using ST = Teko::ST;
  using LO = Teko::LO;
  using GO = Teko::GO;
  using NT = Teko::NT;

  using map_t = Tpetra::Map<LO, GO, NT>;
  using crs_t = Tpetra::CrsMatrix<ST, LO, GO, NT>;
  using vec_t = Tpetra::Vector<ST, LO, GO, NT>;
  using mv_t  = Tpetra::MultiVector<ST, LO, GO, NT>;
  using op_t  = Tpetra::Operator<ST, LO, GO, NT>;

  auto comm = Tpetra::getDefaultComm();

  int myPID = comm->getRank();
  out << "MPI_PID = " << myPID << ", UNIX_PID = " << getpid() << std::endl;

  // Maps
  int dim   = 2;
  GO numVel = 4225, numPre = 4225;
  int errCode = 0;

  RCP<const map_t> mapVel = rcp(new map_t(numVel, 0, comm));
  RCP<const map_t> mapPre = rcp(new map_t(numPre, 0, comm));
  RCP<const map_t> mapAll = rcp(new map_t(numVel * dim + numPre, 0, comm));

  // Reorder
  std::vector<GO> reorderedVec;

  for (LO lid = 0; lid < static_cast<LO>(mapVel->getLocalNumElements()); ++lid) {
    GO gid = mapVel->getGlobalElement(lid);
    for (int i = 0; i < dim; i++) {
      reorderedVec.push_back(gid + numVel * i);
    }
  }

  for (LO lid = 0; lid < static_cast<LO>(mapPre->getLocalNumElements()); ++lid) {
    GO gid = mapPre->getGlobalElement(lid);
    reorderedVec.push_back(gid + numVel * dim);
  }

  RCP<const map_t> mapReorder =
      rcp(new map_t(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
                    Teuchos::ArrayView<const GO>(reorderedVec), 0, comm));

  RCP<Tpetra::Import<LO, GO, NT> > importReorder =
      rcp(new Tpetra::Import<LO, GO, NT>(mapAll, mapReorder));

  std::vector<std::vector<GO> > blockedVec;
  for (int i = 0; i < dim; i++) {
    std::vector<GO> blk;
    for (LO lid = 0; lid < static_cast<LO>(mapVel->getLocalNumElements()); ++lid) {
      GO gid = mapVel->getGlobalElement(lid);
      blk.push_back(gid + numVel * i);
    }
    blockedVec.push_back(blk);
  }
  {
    std::vector<GO> blk;
    for (LO lid = 0; lid < static_cast<LO>(mapPre->getLocalNumElements()); ++lid) {
      GO gid = mapPre->getGlobalElement(lid);
      blk.push_back(gid + numVel * dim);
    }
    blockedVec.push_back(blk);
  }

  // Read matrices and RHS
  RCP<crs_t> ptrMat = Tpetra::MatrixMarket::Reader<crs_t>::readSparseFile("data/szdMat.mm", comm);
  TEST_ASSERT(!ptrMat.is_null());

  RCP<crs_t> ptrMp = Tpetra::MatrixMarket::Reader<crs_t>::readSparseFile("data/szdMp.mm", comm);
  TEST_ASSERT(!ptrMp.is_null());

  LinearOp lpMp = Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrMp->getRangeMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrMp->getDomainMap()), ptrMp);

  RCP<const map_t> mapAllL = mapAll;
  RCP<vec_t> ptrb = Tpetra::MatrixMarket::Reader<crs_t>::readVectorFile("data/szdRHS.mm", comm,
                                                                        mapAllL, false, false);
  TEST_ASSERT(!ptrb.is_null());

  // Reorder matrix
  RCP<crs_t> mat = rcp(new crs_t(mapReorder, ptrMat->getGlobalMaxNumRowEntries()));
  mat->doImport(*ptrMat, *importReorder, Tpetra::INSERT);
  mat->fillComplete();

  // Build AL operator
  double gamma = 0.05;
  Teuchos::RCP<Teko::NS::ALOperator> al =
      Teuchos::rcp(new Teko::NS::ALOperator(blockedVec, mat, lpMp, gamma));

  vec_t x(mapReorder), b(mapReorder);
  b.doImport(*ptrb, *importReorder, Tpetra::INSERT);
  x.putScalar(0.0);

  // Build augmented RHS
  vec_t bAugmented(b);
  al->augmentRHS(b, bAugmented);

  // Build inverse factory
  RCP<Teko::InverseLibrary> invLib           = Teko::InverseLibrary::buildFromStratimikos();
  Teuchos::RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("Amesos2");

  // Build preconditioner factory
  Teuchos::RCP<Teko::NS::ModALPreconditionerFactory> precFact =
      Teuchos::rcp(new Teko::NS::ModALPreconditionerFactory(inverse, lpMp));

  // Build modified AL preconditioner
  precFact->setGamma(gamma);
  Teko::TpetraHelpers::TpetraBlockPreconditioner prec(precFact);
  prec.buildPreconditioner(al);

  // Solve with Belos
  RCP<mv_t> X = rcp(new mv_t(mapReorder, 1));
  RCP<mv_t> B = rcp(new mv_t(mapReorder, 1));

  X->getVectorNonConst(0)->assign(x);
  B->getVectorNonConst(0)->assign(b);

  using problem_t        = Belos::LinearProblem<ST, mv_t, op_t>;
  RCP<problem_t> problem = rcp(new problem_t(al, X, B));

  RCP<op_t> precOp = Teuchos::rcp(&prec, false);
  problem->setRightPrec(precOp);

  const bool set = problem->setProblem();
  TEST_ASSERT(set);

  Teuchos::ParameterList belosList;
  belosList.set("Maximum Iterations", 100);
  belosList.set("Convergence Tolerance", 1e-6);
  belosList.set("Num Blocks", 50);
  belosList.set("Verbosity",
                Belos::Errors + Belos::Warnings + Belos::IterationDetails + Belos::FinalSummary);
  belosList.set("Output Frequency", 10);

  using solver_t = Belos::BlockGmresSolMgr<ST, mv_t, op_t>;
  solver_t solver(problem, rcpFromRef(belosList));

  Belos::ReturnType result = solver.solve();

  int iters = solver.getNumIters();
  out << iters << std::endl;

  if (result == Belos::Converged && iters < 100) {
    out << "GMRES with modified AL preconditioner has converged." << std::endl;
    errCode = 0;
  } else {
    out << "GMRES with modified AL preconditioner has NOT converged." << std::endl;
    errCode = -1;
  }

  TEST_ASSERT(errCode == 0);
}