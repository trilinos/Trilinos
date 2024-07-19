// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Thyra testing tools
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"

// Thyra includes
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraOperatorWrapper.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#include "tBlockedEpetraOperator.hpp"

#include "Teko_BlockedEpetraOperator.hpp"

#include "EpetraExt_RowMatrixOut.h"

namespace Teko {
namespace Test {

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Thyra::createMember;
using Thyra::LinearOpBase;
using Thyra::LinearOpTester;
using Thyra::VectorBase;

void tBlockedEpetraOperator::buildBlockGIDs(std::vector<std::vector<int> >& gids,
                                            const Epetra_Map& map) const {
  int numLocal = map.NumMyElements();
  int numHalf  = numLocal / 2;
  numHalf += ((numHalf % 2 == 0) ? 0 : 1);

  gids.clear();
  gids.resize(3);

  std::vector<int>& blk0 = gids[0];
  std::vector<int>& blk1 = gids[1];
  std::vector<int>& blk2 = gids[2];

  // loop over global IDs: treat first block as strided
  int gid = -1;
  for (int i = 0; i < numHalf; i += 2) {
    gid = map.GID(i);
    blk0.push_back(gid);

    gid = map.GID(i + 1);
    blk1.push_back(gid);
  }

  // loop over global IDs: treat remainder as contiguous
  for (int i = numHalf; i < numLocal; i++) {
    gid = map.GID(i);
    blk2.push_back(gid);
  }

  TEUCHOS_ASSERT(int(blk0.size() + blk1.size() + blk2.size()) == numLocal);
}

void tBlockedEpetraOperator::initializeTest() { tolerance_ = 1e-14; }

int tBlockedEpetraOperator::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                                    int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tBlockedEpetraOperator";

  status = test_vector_constr(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"vector_constr\" ... PASSED", "   \"vector_constr\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 0);
  Teko_TEST_MSG(stdstrm, 1, "   \"reorder(flat reorder)\" ... PASSED",
                "   \"reorder(flat reorder)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 1);
  Teko_TEST_MSG(stdstrm, 1, "   \"reorder(composite reorder = " << 1 << ")\" ... PASSED",
                "   \"reorder(composite reorder)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_reorder(verbosity, failstrm, 2);
  Teko_TEST_MSG(stdstrm, 1, "   \"reorder(composite reorder = " << 2 << ")\" ... PASSED",
                "   \"reorder(composite reorder)\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tBlockedEpetraOperator...PASSED",
                  "tBlockedEpetraOperator...FAILED");
  } else {  // Normal Operating Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tBlockedEpetraOperator...FAILED");
  }

  return failcount;
}

bool tBlockedEpetraOperator::test_vector_constr(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;

  const Epetra_Comm& comm = *GetComm();

  TEST_MSG("\n   tBlockedEpetraOperator::test_vector_constr: "
           << "Running on " << comm.NumProc() << " processors");

  // pick
  int nx = 5 * comm.NumProc();
  int ny = 5 * comm.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm,
                                           false);  // CJ TODO FIXME: change for Epetra64
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  RCP<Epetra_CrsMatrix> A = rcp(FGallery.GetMatrix(), false);
  double beforeNorm       = A->NormOne();

  int width = 3;
  Epetra_MultiVector x(A->OperatorDomainMap(), width);
  Epetra_MultiVector ys(A->OperatorRangeMap(), width);
  Epetra_MultiVector y(A->OperatorRangeMap(), width);

  // EpetraExt::RowMatrixToMatrixMarketFile("fullop.mm",*A);

  std::vector<std::vector<int> > vars;
  buildBlockGIDs(vars, A->RowMap());

  Teko::Epetra::BlockedEpetraOperator shell(vars, A);

  // test the operator against a lot of random vectors
  int numtests = 50;
  double max   = 0.0;
  double min   = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<double> norm(width);
    std::vector<double> rel(width);
    x.Random();

    shell.Apply(x, y);
    A->Apply(x, ys);

    Epetra_MultiVector e(y);
    e.Update(-1.0, ys, 1.0);
    e.Norm2(&norm[0]);

    // compute relative error
    ys.Norm2(&rel[0]);
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "\n   tBlockedEpetraOperator::test_vector_constr: "
                              << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "\n   tBlockedEpetraOperator::test_vector_constr: "
                                     << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  int *indexOffset, *indicies;
  double* values;
  A->ExtractCrsDataPointers(indexOffset, indicies, values);
  for (int i = 0; i < A->NumMyNonzeros(); i++) values[i] *= 2.0;  // square everything!

  double afterNorm = A->NormOne();
  TEST_ASSERT(beforeNorm != afterNorm, "\n   tBlockedEpetraOperator::test_vector_constr "
                                           << toString(status) << ": "
                                           << "verify matrix has been modified");

  shell.RebuildOps();

  // test the operator against a lot of random vectors
  numtests = 50;
  max      = 0.0;
  min      = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<double> norm(width);
    std::vector<double> rel(width);
    x.Random();

    shell.Apply(x, y);
    A->Apply(x, ys);

    Epetra_MultiVector e(y);
    e.Update(-1.0, ys, 1.0);
    e.Norm2(&norm[0]);

    // compute relative error
    ys.Norm2(&rel[0]);
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "\n   tBlockedEpetraOperator::test_vector_constr (rebuild): "
                              << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "\n   tBlockedEpetraOperator::test_vector_constr (rebuild): "
                                     << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  return allPassed;
}

bool tBlockedEpetraOperator::test_reorder(int verbosity, std::ostream& os, int total) {
  bool status    = false;
  bool allPassed = true;

  const Epetra_Comm& comm = *GetComm();

  std::string tstr = total ? "(composite reorder)" : "(flat reorder)";

  TEST_MSG("\n   tBlockedEpetraOperator::test_reorder" << tstr << ": "
                                                       << "Running on " << comm.NumProc()
                                                       << " processors");

  // pick
  int nx = 3 * 25 * comm.NumProc();
  int ny = 3 * 50 * comm.NumProc();

  // create a big matrix to play with
  // note: this matrix is not really strided
  //       however, I just need a nontrivial
  //       matrix to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm,
                                           false);  // CJ TODO FIXME: change for Epetra64
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  RCP<Epetra_CrsMatrix> A = rcp(FGallery.GetMatrix(), false);

  int width = 3;
  Epetra_MultiVector x(A->OperatorDomainMap(), width);
  Epetra_MultiVector yf(A->OperatorRangeMap(), width);
  Epetra_MultiVector yr(A->OperatorRangeMap(), width);

  std::vector<std::vector<int> > vars;
  buildBlockGIDs(vars, A->RowMap());

  Teko::Epetra::BlockedEpetraOperator flatShell(vars, A, "Af");
  Teko::Epetra::BlockedEpetraOperator reorderShell(vars, A, "Ar");

  Teko::BlockReorderManager brm;
  switch (total) {
    case 0:
      brm.SetNumBlocks(3);
      brm.SetBlock(0, 1);
      brm.SetBlock(1, 0);
      brm.SetBlock(2, 2);
      break;
    case 1:
      brm.SetNumBlocks(2);
      brm.SetBlock(0, 1);
      brm.GetBlock(1)->SetNumBlocks(2);
      brm.GetBlock(1)->SetBlock(0, 0);
      brm.GetBlock(1)->SetBlock(1, 2);
      break;
    case 2:
      brm.SetNumBlocks(2);
      brm.GetBlock(0)->SetNumBlocks(2);
      brm.GetBlock(0)->SetBlock(0, 0);
      brm.GetBlock(0)->SetBlock(1, 2);
      brm.SetBlock(1, 1);
      break;
  }
  reorderShell.Reorder(brm);
  TEST_MSG("\n   tBlockedEpetraOperator::test_reorder" << tstr << ": patern = " << brm.toString());

  TEST_MSG("\n   tBlockedEpetraOperator::test_reorder" << tstr << ":\n");
  TEST_MSG("\n      " << Teuchos::describe(*reorderShell.getThyraOp(), Teuchos::VERB_HIGH)
                      << std::endl);

  // test the operator against a lot of random vectors
  int numtests = 10;
  double max   = 0.0;
  double min   = 1.0;
  for (int i = 0; i < numtests; i++) {
    std::vector<double> norm(width);
    std::vector<double> rel(width);
    x.Random();

    flatShell.Apply(x, yf);
    reorderShell.Apply(x, yr);

    Epetra_MultiVector e(yf);
    e.Update(-1.0, yr, 1.0);
    e.Norm2(&norm[0]);

    // compute relative error
    yf.Norm2(&rel[0]);
    for (int j = 0; j < width; j++) {
      max = max > norm[j] / rel[j] ? max : norm[j] / rel[j];
      min = min < norm[j] / rel[j] ? min : norm[j] / rel[j];
    }
  }
  TEST_ASSERT(max >= min, "   tBlockedEpetraOperator::test_reorder"
                              << tstr << ": " << toString(status) << ": "
                              << "sanity checked - " << max << " >= " << min);
  TEST_ASSERT(max <= tolerance_, "   tBlockedEpetraOperator::test_reorder"
                                     << tstr << ": " << toString(status) << ": "
                                     << "testing tolerance over many matrix vector multiplies ( "
                                     << max << " <= " << tolerance_ << " )");

  return allPassed;
}

}  // namespace Test
}  // namespace Teko
