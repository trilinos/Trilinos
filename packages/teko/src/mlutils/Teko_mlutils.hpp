// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_mlutils_hpp__
#define __Teko_mlutils_hpp__

#include "ml_operator.h"
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

#include "Teuchos_RCP.hpp"

#include "Teko_Utilities.hpp"

namespace Teko {

class InverseLibrary;

namespace mlutils {

//! build a very simple row map from the ML_Operator
Teuchos::RCP<Epetra_Map> buildRowMap(ML_Operator *mlOp);

/** convert to an Epetra_CrsMatrix, using a specified row map
 * or the default one build from <code>buildRowMap</code>.
 */
Teuchos::RCP<Epetra_CrsMatrix> convertToCrsMatrix(
    ML_Operator *mlOp, const Teuchos::RCP<Epetra_Map> &rowMap = Teuchos::null);

Teko::LinearOp buildTekoBlockOp(ML_Operator *mlOp, int level);

/** Data structure for teko smoothing information
 */
struct SmootherData {
  Teuchos::RCP<Epetra_Operator> Amat;
  Teuchos::RCP<Epetra_Operator> smootherOperator;
};

/** Smoother function to be send to ML.
 */
int smoother(ML_Smoother *mydata, int leng1, double x[], int leng2, double rhs[]);

extern "C" int ML_Gen_Smoother_Teko(ML *ml, int level, int pre_or_post, int ntimes,
                                    const Teuchos::RCP<const Teuchos::ParameterList> &tekoPL,
                                    const Teuchos::RCP<const Teko::InverseLibrary> &invLib,
                                    const std::string &inverse, bool isBlocked);

}  // namespace mlutils
}  // namespace Teko

#endif
