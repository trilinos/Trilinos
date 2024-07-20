// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "ml_include.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Epetra_Map.h"
#include "Epetra_LinearProblem.h"

#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorIn.h>
#include <EpetraExt_MultiVectorOut.h>
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Teko_EpetraBlockPreconditioner.hpp"
#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_SmootherPreconditionerFactory.hpp"
#include "Teko_mlutils.hpp"

#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"

using namespace Teuchos;

int main(int argc, char* argv[]) {
  int MyPID = 0;
  Epetra_CrsMatrix* A;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  MyPID = Comm.MyPID();
#else
  Epetra_SerialComm Comm;
#endif

  EpetraExt::MatrixMarketFileToCrsMatrix("data/F.mm", Comm, A);
  Teuchos::RCP<Epetra_CrsMatrix> ptrA = Teuchos::rcp(A);

  ParameterList MLMainList;

  ML_Epetra::SetDefaults("SA", MLMainList);

#define ML_OUT_LEVEL 0

  Teuchos::RCP<Teuchos::ParameterList> tekoPL = Teuchos::getParametersFromXmlFile("ml_teko.xml");
  Teuchos::RCP<Teuchos::ParameterList> invLibPL =
      Teuchos::rcpFromRef(tekoPL->sublist("Inverse Library"));

  // big outer parameter list
  MLMainList.set("ML output", 8);
  MLMainList.set("ML print initial list", -2);
  MLMainList.set("ML print final list", -2);
  MLMainList.set("ML label", "1x1 system");
  MLMainList.set("max levels", 3);
  MLMainList.set("PDE equations", 1);
  MLMainList.set("smoother: pre or post", "both");
  MLMainList.set("smoother: sweeps", 2);
  if (argc == 2) {
    MLMainList.set("smoother: type", "Gauss-Seidel");
    MLMainList.set("coarse: type", "Amesos-KLU");
    std::cout << " *** ML GS SMOOTHER *** " << std::endl;
  } else {
    MLMainList.set("smoother: type", "teko");
    MLMainList.set("coarse: type", "teko");
    MLMainList.set("coarse: teko inverse", std::string("Amesos"));
    MLMainList.set("smoother: teko parameter list", invLibPL);
    MLMainList.set("smoother: teko inverse", std::string("GS-Single"));
    MLMainList.set("smoother: teko is blocked", 0);
    std::cout << " *** Teko GS SMOOTHER *** " << std::endl;
  }

  /*********************************************************/
  /* Constructor for composite AMG. Does AMG on A[k] using */
  /* corresponding parameter lists. Then, builds a new AMG */
  /* hierarchy with block diagonal grid transfers. The     */
  /* (k,k)th diagonal block is obtained by extracting the  */
  /* grid transfers associated with AMG on A[k].           */
  /*********************************************************/
  ML_Epetra::MultiLevelPreconditioner* MLPre =
      new ML_Epetra::MultiLevelPreconditioner(*A, MLMainList, true);

  /* read in rhs */
  Epetra_Vector* RHS = new Epetra_Vector(A->OperatorRangeMap());
  RHS->PutScalar(7.0);
  /* set initial guess */
  Epetra_Vector LHS(A->OperatorDomainMap());
  LHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(A, &LHS, RHS);
  AztecOO solver(Problem);

  solver.SetPrecOperator(MLPre);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_conv, AZ_noscaled);
  solver.SetAztecOption(AZ_output, 1);
  solver.Iterate(20, 1e-8);

  delete MLPre;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (EXIT_SUCCESS);
}
