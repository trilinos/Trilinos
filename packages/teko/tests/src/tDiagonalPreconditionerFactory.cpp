// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tDiagonalPreconditionerFactory.hpp"
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

#include <vector>

#include <math.h>

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;

tDiagonalPreconditionerFactory::~tDiagonalPreconditionerFactory() {
  delete fact;
  delete pstate;
  delete[] block_starts;
  delete[] block_gids;
}

void tDiagonalPreconditionerFactory::initializeTest() {
  const Epetra_Comm &comm = *GetComm();

  tolerance_ = 1.0e-14;

  int nx = 39;  // essentially random values
  int ny = 53;

  // create some big blocks to play with
  Trilinos_Util::CrsMatrixGallery FGallery("laplace_2d", comm,
                                           false);  // CJ TODO FIXME: change for Epetra64
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  epetraF = rcp(new Epetra_CrsMatrix(*FGallery.GetMatrix()));
  epetraF->FillComplete(true);
  F_ = Thyra::epetraLinearOp(epetraF);
}

void tDiagonalPreconditionerFactory::buildParameterList(int blocksize) {
  const Epetra_CrsMatrix *F = &*epetraF;
  TEUCHOS_ASSERT(F);

  // Cleanup of last run
  delete[] block_starts;
  delete[] block_gids;

  // Allocs
  int Nr       = F->NumMyRows();
  int Nb       = (int)ceil(((double)Nr) / ((double)blocksize));
  block_starts = new int[Nb + 1];
  block_gids   = new int[Nr];

  // Fill out block data structures
  block_starts[0] = 0;
  for (int i = 0; i < Nb; i++) block_starts[i + 1] = block_starts[i] + blocksize;
  block_starts[Nb] = Nr;

  for (int i = 0; i < Nr; i++) block_gids[i] = F->GRID(i);

  // Set the list
  Teuchos::ParameterList sublist;
  List_.set("number of local blocks", Nb);
  List_.set("block start index", block_starts);
  List_.set("block entry gids", block_gids);
  sublist.set("apply mode", "invert");
  List_.set("blockdiagmatrix: list", sublist);
}

int tDiagonalPreconditionerFactory::runTest(int verbosity, std::ostream &stdstrm,
                                            std::ostream &failstrm, int &totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tDiagonalPreconditionerFactory";

  status = test_createPrec(verbosity, failstrm, 2);
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

bool tDiagonalPreconditionerFactory::test_initializePrec(int verbosity, std::ostream &os) {
  delete pstate;

  pstate = new DiagonalPrecondState();
  pop    = fact->buildPreconditionerOperator(F_, *pstate);

  // const DiagonalPreconditionerOp *dop=dynamic_cast<const DiagonalPreconditionerOp*>(&*pop);
  // if(!dop) return false;
  if (Thyra::get_Epetra_Operator(*pop) == Teuchos::null) return false;

  //  pstate->BDP_->Print(std::cout);
  return true;
}

bool tDiagonalPreconditionerFactory::test_createPrec(int verbosity, std::ostream &os,
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

bool tDiagonalPreconditionerFactory::test_canApply(int verbosity, std::ostream &os) {
  const Epetra_Map &domain_ = epetraF->DomainMap(), &range_ = epetraF->RangeMap();

  RCP<Epetra_Vector> X = rcp(new Epetra_Vector(domain_)), Y = rcp(new Epetra_Vector(range_));
  Epetra_Vector Z(range_);
  X->PutScalar(1.0);
  Y->PutScalar(0.0);
  Z.PutScalar(-5.0);

  // Build Thyra wrappers
  MultiVector tX = Thyra::create_Vector(X, F_->domain());
  MultiVector tY = Thyra::create_Vector(Y, F_->range());

  // Do the apply via thrya
  //   const DiagonalPreconditionerOp *dop=dynamic_cast<const DiagonalPreconditionerOp*>(&*pop);
  //   dop->implicitApply(tX,tY,1.0,0.0);
  Teko::applyOp(pop, tX, tY, 1.0, 0.0);

  // Do the apply via epetra
  pstate->BDP_->ApplyInverse(*X, Z);

  // Compare solutions
  double znrm, ynrm, dnrm;
  Y->Norm2(&ynrm);
  Z.Norm2(&znrm);
  Z.Update(-1.0, *Y, 1.0);
  Z.Norm2(&dnrm);

  if (!epetraF->Comm().MyPID()) std::cout << "||Z-Y||/||Z|| = " << dnrm / znrm << std::endl;
  if (dnrm / znrm > 1e-12)
    return false;
  else
    return true;
}

}  // end namespace Test
}  // end namespace Teko
