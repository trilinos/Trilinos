// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

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

#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_Utilities.hpp"

#define SS_ECHO(ops)      \
  {                       \
    std::stringstream ss; \
    ss << ops;            \
    ECHO(ss.str())        \
  };

using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Thyra::createMember;
using Thyra::LinearOpBase;
using Thyra::LinearOpTester;
using Thyra::VectorBase;

double tolerance = 1e-14;

const Teuchos::RCP<Epetra_Comm> GetComm() {
#ifdef HAVE_MPI
  static Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  static Epetra_SerialComm comm;
#endif

  return Teuchos::rcpFromRef(comm);
}

TEUCHOS_UNIT_TEST(tStridedEpetraOperator, test_numvars_constr) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  SS_ECHO("tStridedEpetraOperator::test_numvars: "
          << "Running on " << comm.NumProc() << " processors");

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
  double beforeNorm       = A->NormOne();

  int vars  = 3;
  int width = 3;
  Epetra_MultiVector x(A->OperatorDomainMap(), width);
  Epetra_MultiVector ys(A->OperatorRangeMap(), width);
  Epetra_MultiVector y(A->OperatorRangeMap(), width);

  Teko::Epetra::StridedEpetraOperator shell(vars, A);

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
  TEST_ASSERT(max >= min);
  TEST_ASSERT(max <= tolerance);

  int *indexOffset, *indicies;
  double* values;
  A->ExtractCrsDataPointers(indexOffset, indicies, values);
  for (int i = 0; i < A->NumMyNonzeros(); i++) values[i] *= 2.0;  // square everything!

  double afterNorm = A->NormOne();
  TEST_ASSERT(beforeNorm != afterNorm);

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
  TEST_ASSERT(max >= min);
  TEST_ASSERT(max <= tolerance);
}

TEUCHOS_UNIT_TEST(tStridedEpetraOperator, test_vector_constr) {
  const Epetra_Comm& comm = *GetComm();

  SS_ECHO("\n   tStridedEpetraOperator::test_vector_constr: "
          << "Running on " << comm.NumProc() << " processors");

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

  double beforeNorm = A->NormOne();

  int width = 3;
  Epetra_MultiVector x(A->OperatorDomainMap(), width);
  Epetra_MultiVector ys(A->OperatorRangeMap(), width);
  Epetra_MultiVector y(A->OperatorRangeMap(), width);

  std::vector<int> vars;
  vars.push_back(2);
  vars.push_back(1);
  Teko::Epetra::StridedEpetraOperator shell(vars, A);

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
  TEST_ASSERT(max >= min);
  TEST_ASSERT(max <= tolerance)

  int *indexOffset, *indicies;
  double* values;
  A->ExtractCrsDataPointers(indexOffset, indicies, values);
  for (int i = 0; i < A->NumMyNonzeros(); i++) values[i] *= 2.0;  // square everything!

  double afterNorm = A->NormOne();
  TEST_ASSERT(beforeNorm != afterNorm);

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
  TEST_ASSERT(max >= min);
  TEST_ASSERT(max <= tolerance);
}

TEUCHOS_UNIT_TEST(tStridedEpetraOperator, test_reorder) {
  const Epetra_Comm& comm = *GetComm();

  for (int total = 0; total < 3; total++) {
    std::string tstr = total ? "(composite reorder)" : "(flat reorder)";

    SS_ECHO("\n   tStridedEpetraOperator::test_reorder" << tstr << ": "
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

    Teko::Epetra::StridedEpetraOperator flatShell(3, A, "Af");
    Teko::Epetra::StridedEpetraOperator reorderShell(3, A, "Ar");

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
    SS_ECHO("\n   tStridedEpetraOperator::test_reorder" << tstr << ": patern = " << brm.toString());

    SS_ECHO("\n   tStridedEpetraOperator::test_reorder" << tstr << ":\n");
    SS_ECHO("\n      " << Teuchos::describe(*reorderShell.getThyraOp(), Teuchos::VERB_HIGH)
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
    TEST_ASSERT(max >= min);
    TEST_ASSERT(max <= tolerance);
  }
}

TEUCHOS_UNIT_TEST(tStridedEpetraOperator, test_print_norm) {
  const Epetra_Comm& comm = *GetComm();

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

  std::vector<int> vars;
  vars.push_back(1);
  vars.push_back(1);
  vars.push_back(1);
  Teko::Epetra::StridedEpetraOperator shell(vars, A);

  std::string normString = shell.PrintNorm();
  *Teko::getOutputStream() << std::endl << normString << std::endl;
}
