// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tParallelInverse.hpp"

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Config.h"
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_StridedEpetraOperator.hpp"

#include "Thyra_EpetraLinearOp.hpp"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Thyra::epetraLinearOp;

void tParallelInverse::initializeTest() { tolerance_ = 1.0e-7; }

void tParallelInverse::loadMatrix() {
  // Read in the matrix, store pointer as an RCP
  Epetra_CrsMatrix* ptrA = 0;
  TEUCHOS_TEST_FOR_EXCEPT(
      EpetraExt::MatrixMarketFileToCrsMatrix("data/lsc_F_2.mm", *GetComm(), ptrA));
  F_ = Thyra::epetraLinearOp(rcp(ptrA));
}

void tParallelInverse::loadStridedMatrix() {
  // Read in the matrix, store pointer as an RCP
  Epetra_CrsMatrix* ptrA = 0;
  std::vector<int> vec(2);
  vec[0] = 1;
  vec[1] = 2;
  TEUCHOS_TEST_FOR_EXCEPT(
      EpetraExt::MatrixMarketFileToCrsMatrix("data/nsjac.mm", *GetComm(), ptrA));
  RCP<Epetra_CrsMatrix> A = rcp(ptrA);

  // Block the linear system using a strided epetra operator
  Teko::Epetra::StridedEpetraOperator sA(vec, A);

  // get 0,0 block
  F_ = Thyra::epetraLinearOp(sA.GetBlock(0, 0));
}

int tParallelInverse::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                              int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tParallelInverse";

  status = test_inverse(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"inverse\" ... PASSED", "   \"inverse\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

#ifdef Teko_ENABLE_DEV_MODE  // so the file nsjac.mm isn't required for release mode
  status = test_stridedInverse(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"stridedInverse\" ... PASSED", "   \"stridedInverse\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;
#endif

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tParallelInverse...PASSED", "tParallelInverse...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tParallelInverse...FAILED");
  }

  return failcount;
}

bool tParallelInverse::test_inverse(int verbosity, std::ostream& os) {
  // bool status = false;
  bool allPassed = true;

  loadMatrix();

  // build an InverseLibrary
  RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

  // build the inverse factory needed by the example preconditioner
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("Amesos");

  Teko::LinearOp inv = invFact->buildInverse(F_);

  return allPassed;
}

bool tParallelInverse::test_stridedInverse(int verbosity, std::ostream& os) {
  // bool status = false;
  bool allPassed = true;

  loadStridedMatrix();

  // build an InverseLibrary
  TEST_MSG("\n   Building inverse library");
  RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();

  // build the inverse factory needed by the example preconditioner
  TEST_MSG("\n   Building inverse factory");
  RCP<const Teko::InverseFactory> invFact = invLib->getInverseFactory("Amesos");

  TEST_MSG("\n   Building inverse");
  Teko::LinearOp inv = invFact->buildInverse(F_);

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
