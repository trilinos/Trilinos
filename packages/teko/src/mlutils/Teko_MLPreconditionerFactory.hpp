// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_MLPrecondtionerFactory_hpp__
#define __Teko_MLPrecondtionerFactory_hpp__

#include "Teko_BlockPreconditionerFactory.hpp"

#include "ml_include.h"
#include "ml_operator.h"

namespace ML_Epetra {
class MultiLevelPreconditioner;
}

namespace Teko {

//! Contains operator internals need for ML
class MLPreconditionerState : public BlockPreconditionerState {
 public:
  MLPreconditionerState();

  //! Build an ML preconditioner object using the set of coarsening parameters
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> constructMLPreconditioner(
      const Teuchos::ParameterList &mainList,
      const std::vector<Teuchos::RCP<const Teuchos::ParameterList> > &coarseningParams);

  // Set functions
  /////////////////////////////////////////////

  //! set ML Comm pointer...it will be automatically cleaned up
  void setMLComm(ML_Comm *comm);

  //! set ML Operator pointer...it will be automatically cleaned up
  void setMLOperator(ML_Operator *op);

  //! Set if ML internals been constructed yet
  void setIsFilled(bool value);

  //! Set matrices to build multigrid hierarcy from
  void setAggregationMatrices(const std::vector<Epetra_RowMatrix *> &diags);

  // Get functions
  /////////////////////////////////////////////

  //! Has this object been filled yet
  bool isFilled() const;

 protected:
  static void cleanup_ML_Comm(ML_Comm *mlComm);
  static void cleanup_ML_Operator(ML_Operator *mlComm);

  bool isFilled_;
  Teuchos::RCP<ML_Comm> mlComm_;    // note that this has to be properly clean up!
  Teuchos::RCP<ML_Operator> mlOp_;  // note that this has to be properly clean up!

  std::vector<Epetra_RowMatrix *> diagonalOps_;  // patterns for setting up aggregation
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> mlPreconditioner_;
};

/** \brief Class that constructs and returns an ML preconditioner object
 *        that is capable of doing block smoothing.
 */
class MLPreconditionerFactory : public BlockPreconditionerFactory {
 public:
  MLPreconditionerFactory();

  /** \brief Function that is called to build the preconditioner
   *        for the linear operator that is passed in.
   */
  virtual LinearOp buildPreconditionerOperator(BlockedLinearOp &blo,
                                               BlockPreconditionerState &state) const;

  /** \brief Function that permits the construction of an arbitrary
   *        PreconditionerState object.
   */
  virtual Teuchos::RCP<PreconditionerState> buildPreconditionerState() const;

  /** \brief This function builds the internals of the preconditioner factory
   *        from a parameter list.
   */
  void initializeFromParameterList(const Teuchos::ParameterList &settings);

 protected:
  /** \brief Fills an ML preconditioner state object
   *
   * \note This function assumes the blocked linear operator is flat
   *       and no nesting has occured. (Each operator is actually a Epetra_CrsMatrix)
   */
  void fillMLPreconditionerState(const BlockedLinearOp &blo, MLPreconditionerState &mlState) const;

  int blockRowCount_;
  std::vector<Teuchos::RCP<const Teuchos::ParameterList> > coarseningParams_;
  Teuchos::ParameterList mainParams_;
};

}  // end namespace Teko

#endif
