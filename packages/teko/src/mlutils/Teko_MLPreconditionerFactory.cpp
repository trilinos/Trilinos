// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_MLPreconditionerFactory.hpp"

#include "Teko_MLLinearOp.hpp"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"

#include "Thyra_get_Epetra_Operator.hpp"

namespace Teko {

/////////////////////////////////////////////////////////////
// MLPreconditionerState
/////////////////////////////////////////////////////////////

MLPreconditionerState::MLPreconditionerState() : isFilled_(false) {}

void MLPreconditionerState::setMLComm(ML_Comm *comm) {
  mlComm_ = Teuchos::rcpWithDealloc(comm, Teuchos::deallocFunctorDelete<ML_Comm>(cleanup_ML_Comm));
}

void MLPreconditionerState::setMLOperator(ML_Operator *op) {
  mlOp_ =
      Teuchos::rcpWithDealloc(op, Teuchos::deallocFunctorDelete<ML_Operator>(cleanup_ML_Operator));
}

void MLPreconditionerState::setIsFilled(bool value) { isFilled_ = value; }

bool MLPreconditionerState::isFilled() const { return isFilled_; }

void MLPreconditionerState::cleanup_ML_Comm(ML_Comm *mlComm) { ML_Comm_Destroy(&mlComm); }

void MLPreconditionerState::cleanup_ML_Operator(ML_Operator *mlOp) { ML_Operator_Destroy(&mlOp); }

void MLPreconditionerState::setAggregationMatrices(const std::vector<Epetra_RowMatrix *> &diags) {
  diagonalOps_ = diags;
}

Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> MLPreconditionerState::constructMLPreconditioner(
    const Teuchos::ParameterList &mainList,
    const std::vector<Teuchos::RCP<const Teuchos::ParameterList> > &coarseningParams) {
  TEUCHOS_ASSERT(isFilled());
  std::vector<Teuchos::ParameterList> cpls(coarseningParams.size());
  for (std::size_t i = 0; i < coarseningParams.size(); i++) cpls[i] = *coarseningParams[i];

  mlPreconditioner_ = rcp(new ML_Epetra::MultiLevelPreconditioner(
      &*mlOp_, mainList, &diagonalOps_[0], &cpls[0], diagonalOps_.size()));

  return mlPreconditioner_;
}

/////////////////////////////////////////////////////////////
// MLPreconditionerFactory
/////////////////////////////////////////////////////////////

MLPreconditionerFactory::MLPreconditionerFactory() {}

LinearOp MLPreconditionerFactory::buildPreconditionerOperator(
    BlockedLinearOp &blo, BlockPreconditionerState &state) const {
  MLPreconditionerState &mlState = Teuchos::dyn_cast<MLPreconditionerState>(state);

  if (not mlState.isFilled()) fillMLPreconditionerState(blo, mlState);
  TEUCHOS_ASSERT(mlState.isFilled());

  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> mlPrecOp =
      mlState.constructMLPreconditioner(mainParams_, coarseningParams_);

  // return Thyra::epetraLinearOp(mlPrecOp);
  return Teuchos::rcp(new MLLinearOp(mlPrecOp));
}

Teuchos::RCP<PreconditionerState> MLPreconditionerFactory::buildPreconditionerState() const {
  return Teuchos::rcp(new MLPreconditionerState());
}

void MLPreconditionerFactory::fillMLPreconditionerState(const BlockedLinearOp &blo,
                                                        MLPreconditionerState &mlState) const {
  TEUCHOS_ASSERT(not mlState.isFilled());

  EpetraExt::CrsMatrix_SolverMap RowMatrixColMapTrans;

  int rowCnt = blockRowCount(blo);
  int colCnt = blockRowCount(blo);

  // construct comm object...add to state
  ML_Comm *mlComm;
  ML_Comm_Create(&mlComm);
  mlState.setMLComm(mlComm);

  ML_Operator *mlBlkMat = ML_Operator_Create(mlComm);
  ML_Operator_BlkMatInit(mlBlkMat, mlComm, rowCnt, colCnt, ML_DESTROY_EVERYTHING);

  std::vector<Epetra_RowMatrix *> aggMats;
  for (int r = 0; r < rowCnt; r++) {
    for (int c = 0; c < colCnt; c++) {
      ML_Operator *tmp       = 0;
      Epetra_RowMatrix *EMat = 0;
      std::stringstream ss;

      ss << "fine_" << r << "_" << c;

      // extract crs matrix
      Teuchos::RCP<Epetra_CrsMatrix> crsMat = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(
          Thyra::get_Epetra_Operator(*blo->getNonconstBlock(r, c)));
      EMat = ML_Epetra::ModifyEpetraMatrixColMap(*crsMat, RowMatrixColMapTrans, ss.str().c_str());
      if (r == c)  // setup diagonal scaling matrices
        aggMats.push_back(EMat);

      // extract ml sub operator
      tmp = ML_Operator_Create(mlComm);
      ML_Operator_WrapEpetraMatrix(EMat, tmp);

      // add it to the block
      ML_Operator_BlkMatInsert(mlBlkMat, tmp, r, c);
    }
  }
  ML_Operator_BlkMatFinalize(mlBlkMat);

  // finish setting up state object
  mlState.setMLOperator(mlBlkMat);

  // now set aggregation matrices
  mlState.setAggregationMatrices(aggMats);

  mlState.setIsFilled(true);  // register state object as filled
}

/** \brief This function builds the internals of the preconditioner factory
 *        from a parameter list.
 */
void MLPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList &settings) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teko_DEBUG_SCOPE("MLPreconditionerFactory::initializeFromParameterList", 10);

  // clear initial state
  coarseningParams_.clear();
  blockRowCount_ = 0;

  blockRowCount_ = settings.get<int>("Block Row Count");

  // read in main parameter list: with smoothing information
  ////////////////////////////////////////////
  mainParams_ = settings.sublist("Smoothing Parameters");
  mainParams_.set<Teuchos::RCP<const Teko::InverseLibrary> >("smoother: teko inverse library",
                                                             getInverseLibrary());

  // read in aggregation sub lists
  ////////////////////////////////////////////
  const Teuchos::ParameterList &aggSublist = settings.sublist("Block Aggregation");

  for (int block = 0; block < blockRowCount_; block++) {
    // write out sub list name: "Block #"
    std::stringstream ss;
    ss << "Block " << block;
    std::string sublistName = ss.str();

    // grab sublist
    RCP<Teuchos::ParameterList> userSpec =
        rcp(new Teuchos::ParameterList(aggSublist.sublist(sublistName)));
    coarseningParams_.push_back(userSpec);
  }
}

}  // end namespace Teko
