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
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_VectorIn.h>
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

#define CMP2RAY
// #undef CMP2RAY

EpetraExt::CrsMatrix_SolverMap RowMatrixColMapTrans_;

ML_Operator *readNS(ML_Comm *ml_comm, Epetra_Comm &Comm);

int main(int argc, char *argv[]) {
  int *idummy1, *idummy2;
  double *ddummy;
  int MyPID = 0;
  Epetra_RowMatrix *A[2];

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  MyPID = Comm.MyPID();
#else
  Epetra_SerialComm Comm;
#endif
  ML_Comm *ml_comm;
  ML_Comm_Create(&ml_comm);

  Teuchos::RCP<Teuchos::ParameterList> tekoPL = Teuchos::getParametersFromXmlFile("ml_teko.xml");
  bool cForPressure                           = tekoPL->get<bool>("Use C For Pressure");

  /***********************************************/
  /* Read in BlkMatrix, set A[0] to BlkMat(0,0), */
  /* and then read BBt into A[1].                */
  /***********************************************/

  ML_Operator *BlkMat = readNS(ml_comm, Comm);
  ML_Operator *tmpF   = ML_Operator_BlkMatExtract(BlkMat, 0, 0);
  ML_Operator *tmpC   = ML_Operator_BlkMatExtract(BlkMat, 1, 1);

  // CRS or ROW matrix?
  if (Epetra_ML_GetCrsDataptrs(tmpF, &ddummy, &idummy1, &idummy2))
    A[0] = (Epetra_RowMatrix *)tmpF->data;
  else
    A[0] = dynamic_cast<Epetra_RowMatrix *>((Epetra_CrsMatrix *)tmpF->data);

  Teko::ModifiableLinearOp BBt;
  if (cForPressure) {
    std::cout << "Using C for SA" << std::endl;
    if (Epetra_ML_GetCrsDataptrs(tmpC, &ddummy, &idummy1, &idummy2))
      A[1] = (Epetra_RowMatrix *)tmpC->data;
    else
      A[1] = dynamic_cast<Epetra_RowMatrix *>((Epetra_CrsMatrix *)tmpC->data);
  } else {
    std::cout << "Using approximate Schur Complement for SA" << std::endl;
    ML_Operator *tmpBt = ML_Operator_BlkMatExtract(BlkMat, 0, 1);
    ML_Operator *tmpB  = ML_Operator_BlkMatExtract(BlkMat, 1, 0);

    Teko::LinearOp Bt = Thyra::epetraLinearOp(Teuchos::rcp((Epetra_CrsMatrix *)tmpBt->data, false));
    Teko::LinearOp B  = Thyra::epetraLinearOp(Teuchos::rcp((Epetra_CrsMatrix *)tmpB->data, false));
    Teko::LinearOp F  = Thyra::epetraLinearOp(Teuchos::rcp((Epetra_CrsMatrix *)tmpF->data, false));
    Teko::LinearOp C  = Thyra::epetraLinearOp(Teuchos::rcp((Epetra_CrsMatrix *)tmpC->data, false));
    Teko::LinearOp idF = Teko::getInvDiagonalOp(F, Teko::AbsRowSum);

    BBt  = Teko::explicitAdd(C, Teko::scale(-1.0, Teko::explicitMultiply(B, idF, Bt)), BBt);
    A[1] = &*Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*BBt));
  }

  /**************************************************/
  /* Read ML parameter lists for subblocks and for  */
  /* block 2x2 system.                              */
  /**************************************************/

  Teuchos::RCP<Teuchos::ParameterList> invLibPL =
      Teuchos::rcpFromRef(tekoPL->sublist("Inverse Library"));
  int maxAztecIters = tekoPL->get<int>("Max Aztec Iters", 100);
  double aztecTols  = tekoPL->get<double>("Aztec Tolerance", 1e-8);

  Teuchos::ParameterList MLMainList;
  Teuchos::ParameterList MLList[2];
  ML_Epetra::SetDefaults("SA", MLMainList);
  ML_Epetra::SetDefaults("SA", MLList[0]);
  ML_Epetra::SetDefaults("SA", MLList[1]);

  MLMainList.setParameters(tekoPL->sublist("Main ML List"));
  MLList[0].setParameters(tekoPL->sublist("u ML List"));
  MLList[1].setParameters(tekoPL->sublist("p ML List"));

  MLMainList.set("smoother: teko parameter list", invLibPL);

  /*********************************************************/
  /* Constructor for composite AMG. Does AMG on A[k] using */
  /* corresponding parameter lists. Then, builds a new AMG */
  /* hierarchy with block diagonal grid transfers. The     */
  /* (k,k)th diagonal block is obtained by extracting the  */
  /* grid transfers associated with AMG on A[k].           */
  /*********************************************************/
  ML_Epetra::MultiLevelPreconditioner *MLPre =
      new ML_Epetra::MultiLevelPreconditioner(BlkMat, MLMainList, A, MLList, 2);
  // ML * mlptr = const_cast<ML*>(MLPre->GetML());
  // ML_Operator_Print(&(mlptr->Pmat[1]),"CompP1");
  // ML_Operator_Print(&(mlptr->Rmat[0]),"CompR0");

  /* convert BlkMat to an Epetra BlkMat */
  ML_Epetra::RowMatrix EBlkMat(BlkMat, &Comm);

  /* read in rhs */
  Teuchos::RCP<Epetra_Vector> RHS = Teuchos::rcp(new Epetra_Vector(EBlkMat.OperatorRangeMap()));
  RHS->PutScalar(7.0);
  Teuchos::RCP<Epetra_Vector> srcRHS = Teuchos::rcp(new Epetra_Vector(*RHS));
  EBlkMat.Apply(*srcRHS, *RHS);

  /* set initial guess */
  Epetra_Vector LHS(EBlkMat.OperatorDomainMap());
  LHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(&EBlkMat, &LHS, &*RHS);
  AztecOO solver(Problem);

  solver.SetPrecOperator(MLPre);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 1);
  solver.Iterate(maxAztecIters, aztecTols);

  delete MLPre;
  ML_Operator_Destroy(&BlkMat);
  ML_Comm_Destroy(&ml_comm);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return (EXIT_SUCCESS);
}

/***********************************************************/
/* Little temporary utility to read in a bunch of Epetra   */
/* CRS matrices and store the whole thing within a block   */
/* ML matrix.                                              */
/***********************************************************/

ML_Operator *convertToSubOp(ML_Comm *ml_comm, Epetra_Comm &Comm, const std::string &fileName) {
  ML_Operator *tmp         = 0;
  const char *str          = fileName.c_str();
  Epetra_RowMatrix *EMat   = 0;
  Epetra_CrsMatrix *crsMat = 0;

  int finfo = EpetraExt::MatrixMarketFileToCrsMatrix(str, Comm, crsMat);

  if (finfo == 0) {
    EMat = ML_Epetra::ModifyEpetraMatrixColMap(*crsMat, RowMatrixColMapTrans_, str);
    tmp  = ML_Operator_Create(ml_comm);
    ML_Operator_WrapEpetraMatrix(EMat, tmp);
  } else
    TEUCHOS_ASSERT(false);

  return tmp;
}

ML_Operator *readNS(ML_Comm *ml_comm, Epetra_Comm &Comm) {
  ML_Operator *BlkMat = ML_Operator_Create(ml_comm);
  ML_Operator_BlkMatInit(BlkMat, ml_comm, 2, 2, ML_DESTROY_EVERYTHING);
  ML_Operator_BlkMatInsert(BlkMat, convertToSubOp(ml_comm, Comm, "data/F.mm"), 0, 0);
  ML_Operator_BlkMatInsert(BlkMat, convertToSubOp(ml_comm, Comm, "data/C.mm"), 1, 1);
  ML_Operator_BlkMatInsert(BlkMat, convertToSubOp(ml_comm, Comm, "data/Bt.mm"), 0, 1);
  ML_Operator_BlkMatInsert(BlkMat, convertToSubOp(ml_comm, Comm, "data/B.mm"), 1, 0);
  ML_Operator_BlkMatFinalize(BlkMat);

  return BlkMat;
}
