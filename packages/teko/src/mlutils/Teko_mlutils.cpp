// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_mlutils.hpp"

#include <vector>

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"

#include "ml_epetra_utils.h"
#include "ml_op_utils.h"

#include "Thyra_EpetraLinearOp.hpp"

#include "EpetraExt_RowMatrixOut.h"

#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerFactory.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"
#include "Teko_SmootherPreconditionerFactory.hpp"

namespace Teko {
namespace mlutils {

// build a very simple row map from the ML_Operator
Teuchos::RCP<Epetra_Map> buildRowMap(ML_Operator *mlOp) {
  ML_Comm *comm = mlOp->comm;
#ifdef ML_MPI
  MPI_Comm mpi_comm;
  mpi_comm = comm->USR_comm;
  Epetra_MpiComm eComm(mpi_comm);
#else
  Epetra_SerialComm eComm;
#endif

  // get operators row count, build map...who cares what GIDs are
  int rowCnt = mlOp->outvec_leng;
  return Teuchos::rcp(new Epetra_Map(-1, rowCnt, 0, eComm));
}

/** convert to an Epetra_CrsMatrix, using a specified row map
 * or the default one build from <code>buildRowMap</code>.
 */
Teuchos::RCP<Epetra_CrsMatrix> convertToCrsMatrix(ML_Operator *mlOp,
                                                  const Teuchos::RCP<Epetra_Map> &rowMapArg) {
  ML_Comm *comm = mlOp->comm;
#ifdef ML_MPI
  MPI_Comm mpi_comm;
  mpi_comm = comm->USR_comm;
  Epetra_MpiComm eComm(mpi_comm);
#else
  Epetra_SerialComm eComm;
#endif

  // build a default row map
  Teuchos::RCP<Epetra_Map> rowMap = rowMapArg;
  if (rowMapArg == Teuchos::null) rowMap = buildRowMap(mlOp);

  // build lightweight CrsMatrix wrapper
  Epetra_CrsMatrix *crsMatrixPtr = 0;
  int maxNumNonzeros             = 0;
  double cpuTime                 = 0.0;
  bool verbose                   = false;
  ML_Operator2EpetraCrsMatrix(mlOp, crsMatrixPtr, maxNumNonzeros, false, cpuTime, verbose);

  return Teuchos::rcp(crsMatrixPtr);
}

Teko::LinearOp buildTekoBlockOp(ML_Operator *mlOp, int level) {
  Teko_DEBUG_SCOPE("Teko::mlutils::buildTekoBlockOp", 0);

  using Teuchos::RCP;

  int numRows = ML_Operator_BlkMatNumBlockRows(mlOp);
  int numCols = ML_Operator_BlkMatNumBlockCols(mlOp);

  Teko::BlockedLinearOp tekoOp = Teko::createBlockedOp();
  Teko::beginBlockFill(tekoOp, numRows, numCols);
  for (int i = 0; i < numRows; i++) {
    for (int j = 0; j < numCols; j++) {
      ML_Operator *subBlock = ML_Operator_BlkMatExtract(mlOp, i, j);

      // if required construct and add a block to tekoOp
      if (subBlock != 0) {
        // create a CRS matrix from ML operator
        RCP<const Epetra_CrsMatrix> eCrsMat = convertToCrsMatrix(subBlock);

#if 0
               std::stringstream ss;
               ss << "A(" << level << ")_" << i << "_" << j << ".mm";
               EpetraExt::RowMatrixToMatrixMarketFile(ss.str().c_str(),*eCrsMat);
#endif

        // build Teko operator, add to Teko blocked operator
        Teko::LinearOp tekoSubBlock = Thyra::epetraLinearOp(eCrsMat);
        Teko::setBlock(i, j, tekoOp, tekoSubBlock);
      }
    }
  }
  Teko::endBlockFill(tekoOp);

  return tekoOp.getConst();
}

int smoother(ML_Smoother *mydata, int leng1, double x[], int leng2, double rhs[]) {
  Teko_DEBUG_SCOPE("Teko::mlutils::smoother", 10);

  // std::cout << "init guess = " << mydata->init_guess << std::endl;

  // grab data object
  SmootherData *smootherData = (struct SmootherData *)ML_Get_MySmootherData(mydata);

  Epetra_Vector X(View, smootherData->Amat->OperatorDomainMap(), x);
  Epetra_Vector Y(View, smootherData->Amat->OperatorRangeMap(), rhs);

  smootherData->smootherOperator->Apply(Y, X);

  return 0;
}

Teuchos::RCP<Teko::InverseFactory> ML_Gen_Init_Teko_Prec(const std::string smoothName,
                                                         const Teko::InverseLibrary &invLib) {
  Teko_DEBUG_SCOPE("Teko::mlutils::ML_Gen_Init_Teko_Prec", 10);

  // Teuchos::RCP<Teko::RequestHandler> rh = Teuchos::rcp(new Teko::RequestHandler());
  // Teuchos::RCP<Teko::InverseLibrary> invLib
  //       = Teko::InverseLibrary::buildFromParameterList(pl);
  // invLib->setRequestHandler(rh);

  Teuchos::RCP<Teko::PreconditionerFactory> precFact;
  Teuchos::RCP<Teko::InverseFactory> invFact = invLib.getInverseFactory(smoothName);

  return invFact;
}

extern "C" void ML_Destroy_Smoother_Teko(void *data) {
  Teko::mlutils::SmootherData *sData = (Teko::mlutils::SmootherData *)data;
  delete sData;
}

extern "C" int ML_Gen_Smoother_Teko(ML *ml, int level, int pre_or_post, int ntimes,
                                    const Teuchos::RCP<const Teuchos::ParameterList> &tekoPL,
                                    const Teuchos::RCP<const Teko::InverseLibrary> &invLib_in,
                                    const std::string &inverse, bool isBlocked) {
  Teko_DEBUG_SCOPE("Teko::mlutils::ML_Gen_Smoother_Teko", 10);
  ML_Operator *BlkMat = &(ml->Amat[level]);

  Teuchos::RCP<const Teko::InverseLibrary> invLib;
  if (invLib_in == Teuchos::null) {
    Teuchos::RCP<Teko::RequestHandler> rh = Teuchos::rcp(new Teko::RequestHandler());
    Teuchos::RCP<Teko::InverseLibrary> ncInvLib =
        Teko::InverseLibrary::buildFromParameterList(*tekoPL);
    ncInvLib->setRequestHandler(rh);

    invLib = ncInvLib.getConst();
  } else {
    // this assumes a request handler has already been set
    invLib = invLib_in;
  }

  Teko::LinearOp tekoAmat;
  if (isBlocked) {
    // build teko blocked operator
    tekoAmat = Teko::mlutils::buildTekoBlockOp(BlkMat, level);
  } else {
    // build unblocked operator
    Teuchos::RCP<const Epetra_CrsMatrix> eCrsMat = convertToCrsMatrix(BlkMat);
    tekoAmat                                     = Thyra::epetraLinearOp(eCrsMat);
  }

  // Teuchos::RCP<Teko::mlutils::SmootherData> data = Teuchos::rcp(new Teko::mlutils::SmootherData);
  Teko::mlutils::SmootherData *data = new Teko::mlutils::SmootherData;
  data->Amat = Teuchos::rcp(new Teko::Epetra::EpetraOperatorWrapper(tekoAmat));

  // Teuchos::ParameterList pl = *Teuchos::getParametersFromXmlFile(filename);
  Teuchos::RCP<Teko::InverseFactory> invFact = ML_Gen_Init_Teko_Prec(inverse, *invLib);
  Teuchos::RCP<Teko::RequestHandler> rh      = invFact->getRequestHandler();

  // build smoother operator
  Teko::LinearOp precOp     = Teko::buildInverse(*invFact, tekoAmat);
  Teko::LinearOp smootherOp = Teko::buildSmootherLinearOp(tekoAmat, precOp, ntimes, true);

  // get an epetra operator wrapper
  data->smootherOperator = Teuchos::rcp(new Teko::Epetra::EpetraOperatorWrapper(smootherOp));

  int ret_val =
      ML_Set_Smoother(ml, level, pre_or_post, (void *)data, Teko::mlutils::smoother, NULL);
  ml->post_smoother[level].data_destroy = ML_Destroy_Smoother_Teko;

  return ret_val;
}

}  // namespace mlutils
}  // namespace Teko
