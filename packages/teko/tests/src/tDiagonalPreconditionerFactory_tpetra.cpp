// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tDiagonalPreconditionerFactory_tpetra.hpp"
#include "Teko_DiagonalPreconditionerFactory.hpp"
#include "Teko_DiagonalPreconditionerOp.hpp"
#include "EpetraExt_PointToBlockDiagPermute.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "Teko_Utilities.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include <vector>

#include <math.h>

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;

tDiagonalPreconditionerFactory_tpetra::~tDiagonalPreconditionerFactory_tpetra() {
  delete fact;
  delete pstate;
  delete[] block_starts;
  delete[] block_gids;
}

void tDiagonalPreconditionerFactory_tpetra::initializeTest() {
  const Epetra_Comm &comm_epetra                   = *GetComm();
  const RCP<const Teuchos::Comm<int> > comm_tpetra = GetComm_tpetra();

  tolerance_ = 1.0e-14;

  int nx = 39;  // essentially random values
  int ny = 53;

  // create some big blocks to play with
  Trilinos_Util::CrsMatrixGallery FGallery("laplace_2d", comm_epetra, false);
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  RCP<Epetra_CrsMatrix> epetraF = rcp(new Epetra_CrsMatrix(*FGallery.GetMatrix()));
  epetraF->FillComplete(true);
  tpetraF = Teko::TpetraHelpers::epetraCrsMatrixToTpetra(epetraF, comm_tpetra);
  F_      = Thyra::constTpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(tpetraF->getRangeMap()), tpetraF);
}

void tDiagonalPreconditionerFactory_tpetra::buildParameterList(int blocksize) {
  const Tpetra::CrsMatrix<ST, LO, GO, NT> *F = &*tpetraF;
  TEUCHOS_ASSERT(F);

  if (blocksize > 0) {
    // Cleanup of last run
    delete[] block_starts;
    delete[] block_gids;

    // Allocs
    LO Nr        = F->getLocalNumRows();
    int Nb       = (int)ceil(((double)Nr) / ((double)blocksize));
    block_starts = new GO[Nb + 1];
    block_gids   = new GO[Nr];

    // Fill out block data structures
    block_starts[0] = 0;
    for (int i = 0; i < Nb; i++) block_starts[i + 1] = block_starts[i] + blocksize;
    block_starts[Nb] = Nr;

    for (int i = 0; i < Nr; i++) block_gids[i] = F->getRowMap()->getGlobalElement(i);

    // Set the list
    Teuchos::ParameterList sublist;
    List_.set("number of local blocks", Nb);
    List_.set("block start index", block_starts);
    List_.set("block entry gids", block_gids);
    sublist.set("apply mode", "invert");
    List_.set("blockdiagmatrix: list", sublist);
  } else {
    List_.set("Diagonal Type", "Diagonal");
  }
}

int tDiagonalPreconditionerFactory_tpetra::runTest(int verbosity, std::ostream &stdstrm,
                                                   std::ostream &failstrm, int &totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tDiagonalPreconditionerFactory_tpetra";

  status = test_createPrec(verbosity, failstrm, 0);
  Teko_TEST_MSG(stdstrm, 1, "   \"createPrec\" ... PASSED", "   \"createPrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_initializePrec(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"initializePrec\" ... PASSED", "   \"initializePrec\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_canApply(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"canApply\" ... PASSED", "   \"canApply\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tDiagonalPreconditionedFactory...PASSED",
                  "tDiagonalPreconditionedFactory...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tDiagonalPreconditionedFactory...FAILED");
  }

  return failcount;
}

bool tDiagonalPreconditionerFactory_tpetra::test_initializePrec(int verbosity, std::ostream &os) {
  delete pstate;

  pstate = new DiagonalPrecondState();
  pop    = fact->buildPreconditionerOperator(F_, *pstate);

  // Check that a diagonal linear op was produced
  RCP<const Thyra::DiagonalLinearOpBase<ST> > dop =
      rcp_dynamic_cast<const Thyra::DiagonalLinearOpBase<ST> >(pop);
  if (dop.is_null()) return false;

  return true;
}

bool tDiagonalPreconditionerFactory_tpetra::test_createPrec(int verbosity, std::ostream &os,
                                                            int blocksize) {
  buildParameterList(blocksize);

  // Cleanup
  delete fact;

  // New preconditioner
  fact = new DiagonalPreconditionerFactory();
  fact->initializeFromParameterList(List_);
  if (!fact) return false;

  return true;
}

bool tDiagonalPreconditionerFactory_tpetra::test_canApply(int verbosity, std::ostream &os) {
  RCP<const Tpetra::Map<LO, GO, NT> > domain_ = tpetraF->getDomainMap();
  RCP<const Tpetra::Map<LO, GO, NT> > range_  = tpetraF->getRangeMap();

  RCP<Tpetra::Vector<ST, LO, GO, NT> > X = Tpetra::createVector<ST, LO, GO, NT>(domain_);
  RCP<Tpetra::Vector<ST, LO, GO, NT> > Y = Tpetra::createVector<ST, LO, GO, NT>(range_);
  RCP<Tpetra::Vector<ST, LO, GO, NT> > Z = Tpetra::createVector<ST, LO, GO, NT>(range_);
  Y->putScalar(0.0);
  Z->putScalar(1.0);

  // Let X = diag(F). Then applying the preconditioner to X should yield a vector of ones
  tpetraF->getLocalDiagCopy(*X);

  // Build Thyra wrappers
  MultiVector tX =
      Thyra::createVector<ST, LO, GO, NT>(X, Thyra::createVectorSpace<ST, LO, GO, NT>(domain_));
  MultiVector tY =
      Thyra::createVector<ST, LO, GO, NT>(Y, Thyra::createVectorSpace<ST, LO, GO, NT>(range_));

  // Do the apply via thyra
  Teko::applyOp(pop, tX, tY, 1.0, 0.0);

  // Compare solutions
  double znrm, dnrm;
  znrm = Z->norm2();
  Z->update(-1.0, *Y, 1.0);
  dnrm = Z->norm2();

  if (!tpetraF->getComm()->getRank()) std::cout << "||Z-Y||/||Z|| = " << dnrm / znrm << std::endl;
  if (dnrm / znrm > 1e-12)
    return false;
  else
    return true;
}

}  // end namespace Test
}  // end namespace Teko
