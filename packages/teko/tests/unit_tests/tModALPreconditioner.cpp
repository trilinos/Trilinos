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

// Teuchos
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

// Epetra
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Export.h"

// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"

#include "ml_epetra_preconditioner.h"

// Thyra
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp_decl.hpp"

// Teko
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_SIMPLEPreconditionerFactory.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_StridedEpetraOperator.hpp"
#include "Teko_EpetraBlockPreconditioner.hpp"
#include "Teko_BlockedMappingStrategy.hpp"
#include "Teko_BlockedEpetraOperator.hpp"
#include "Teko_BlockUpperTriInverseOp.hpp"
#include "Teko_BlockedMappingStrategy.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

// AL.
#include "Teko_ALOperator.hpp"
#include "Teko_InvModALStrategy.hpp"
#include "Teko_ModALPreconditionerFactory.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

using Teuchos::null;
using Teuchos::ParameterList;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using namespace Thyra;
using namespace Teko;
using namespace Teko::Epetra;
using namespace std;

TEUCHOS_UNIT_TEST(tModALPreconditioner, test) {
  // build global (or serial communicator
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // get process information
  int myPID = comm.MyPID();
  out << "MPI_PID = " << myPID << ", UNIX_PID = " << getpid() << std::endl;

  // Maps.
  // Stable finite elements.
  // int dim = 2, numVel = 4225, numPre = 1089, errCode;
  // Stabilized finite elements.
  int dim = 2, numVel = 4225, numPre = 4225, errCode;
  Epetra_Map mapVel(numVel, 0, comm), mapPre(numPre, 0, comm),
      mapAll(numVel * dim + numPre, 0, comm);

  // Reorder.
  std::vector<int> reorderedVec;
  int numMyLen = mapVel.NumMyElements();
  int *myGlb;
  myGlb = mapVel.MyGlobalElements();
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < numMyLen; j++) reorderedVec.push_back(myGlb[j] + numVel * i);
  numMyLen = mapPre.NumMyElements();
  myGlb    = mapPre.MyGlobalElements();
  for (int j = 0; j < numMyLen; j++) reorderedVec.push_back(myGlb[j] + numVel * dim);

  Teuchos::RCP<Epetra_Map> mapReorder =
      Teuchos::rcp(new Epetra_Map(-1, reorderedVec.size(), &reorderedVec[0], 0, comm));
  Teuchos::RCP<Epetra_Import> importReorder = Teuchos::rcp(new Epetra_Import(*mapReorder, mapAll));

  std::vector<std::vector<int> > blockedVec;
  numMyLen = mapVel.NumMyElements();
  myGlb    = mapVel.MyGlobalElements();
  for (int i = 0; i < dim; i++) {
    reorderedVec.clear();
    for (int j = 0; j < numMyLen; j++) reorderedVec.push_back(myGlb[j] + numVel * i);
    blockedVec.push_back(reorderedVec);
  }
  numMyLen = mapPre.NumMyElements();
  myGlb    = mapPre.MyGlobalElements();
  reorderedVec.clear();
  for (int j = 0; j < numMyLen; j++) reorderedVec.push_back(myGlb[j] + numVel * dim);
  blockedVec.push_back(reorderedVec);

  // Read matrices and right-hand side.
  Epetra_CrsMatrix *ptrMat = 0, *ptrMp = 0;
  Epetra_Vector *ptrb = 0;
  TEUCHOS_ASSERT(EpetraExt::MatrixMarketFileToCrsMatrix("data/szdMat.mm", mapAll, ptrMat) == 0);
  TEUCHOS_ASSERT(EpetraExt::MatrixMarketFileToCrsMatrix("data/szdMp.mm", mapPre, ptrMp) == 0);
  LinearOp lpMp = Thyra::epetraLinearOp(Teuchos::rcpFromRef(*ptrMp));
  TEUCHOS_ASSERT(EpetraExt::MatrixMarketFileToVector("data/szdRHS.mm", mapAll, ptrb) == 0);

  // Use Epetra_Import to reorder.
  Teuchos::RCP<Epetra_CrsMatrix> mat =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *mapReorder, ptrMat->GlobalMaxNumEntries()));
  errCode = mat->Import(*ptrMat, *importReorder, Insert);
  errCode = mat->FillComplete();

  // Build AL operator.
  double gamma = 0.05;
  Teuchos::RCP<Teko::NS::ALOperator> al =
      Teuchos::rcp(new Teko::NS::ALOperator(blockedVec, mat, lpMp, gamma));

  Epetra_Vector x(*mapReorder, false), b(*mapReorder, false);
  b.Import(*ptrb, *importReorder, Insert);
  x.PutScalar(0.0);

  // Build augmented right-hand side.
  Epetra_Vector bAugmented(b);
  al->augmentRHS(b, bAugmented);

  // Build an InverseLibrary
  RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
  // RCP<Teuchos::ParameterList> ml = Teuchos::getParametersFromXmlFile ("ml.xml");
  // RCP<InverseLibrary> invLib = rcp(new InverseLibrary());
  // Teuchos::ParameterList pl;
  // pl.set("Type","ML");
  // Teuchos::ParameterList & pls = pl.sublist("ML Settings");
  // pls.setParameters(*ml);
  // invLib->addInverse("ML",pl);
  Teuchos::RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("Amesos");

  // Build the preconditioner factory
  // RCP<Teko::NS::InvModALStrategy> strategy = rcp(new Teko::NS::InvModALStrategy(inverse));
  // strategy->setParameters(lpMp, gamma);
  Teuchos::RCP<Teko::NS::ModALPreconditionerFactory> precFact =
      Teuchos::rcp(new Teko::NS::ModALPreconditionerFactory(inverse, lpMp));

  // Build modified AL preconditioner.
  precFact->setGamma(gamma);
  Teko::Epetra::EpetraBlockPreconditioner prec(precFact);
  prec.buildPreconditioner(al);  // fill epetra preconditioner using the strided operator

  // Solve.
  AztecOO solver;
  solver.SetUserOperator(&*al);
  solver.SetRHS(&b);
  solver.SetLHS(&x);
  solver.SetPrecOperator(&prec);
  solver.Iterate(100, 1e-6);
  out << solver.NumIters() << std::endl;

  if (solver.NumIters() < 100) {
    out << "GMRES with modified AL preconditioner has converged." << std::endl;
    errCode = 0;
  } else {
    out << "GMRES with modified AL preconditioner has NOT converged." << std::endl;
    errCode = -1;
  }

  TEST_ASSERT(errCode == 0);
}
