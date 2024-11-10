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

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

TEUCHOS_UNIT_TEST(tExplicitOps, transpose) {
// build global (or serial communicator)
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int nx = 39;  // essentially random values
  int ny = 53;

  // create some big blocks to play with
  Trilinos_Util::CrsMatrixGallery FGallery("recirc_2d", comm,
                                           false);  // CJ TODO FIXME: change for Epetra64
  FGallery.Set("nx", nx);
  FGallery.Set("ny", ny);
  Epetra_CrsMatrix& epetraF = FGallery.GetMatrixRef();
  Teko::LinearOp F          = Thyra::epetraLinearOp(rcp(new Epetra_CrsMatrix(epetraF)));
  Teko::LinearOp F_T        = Teko::explicitTranspose(F);
  Teko::LinearOp aF         = Thyra::adjoint(F);

  Teuchos::RCP<const Epetra_Operator> eOp = Thyra::get_Epetra_Operator(*F_T);
  TEST_ASSERT(eOp != Teuchos::null);
  Teuchos::RCP<const Epetra_CrsMatrix> eCrs =
      Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(eOp);
  TEST_ASSERT(eCrs != Teuchos::null);

  Thyra::LinearOpTester<double> tester;
  tester.set_all_error_tol(1e-14);
  tester.show_all_tests(true);

  {
    std::stringstream ss;
    const bool result = tester.compare(*aF, *F_T, Teuchos::ptrFromRef(out));
    TEST_ASSERT(result);
  }
}
