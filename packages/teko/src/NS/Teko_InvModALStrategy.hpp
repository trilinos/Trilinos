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

#ifndef __Teko_ModALStrategy_hpp__
#define __Teko_ModALStrategy_hpp__

#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_BlockPreconditionerFactory.hpp"

namespace Teko {

namespace NS {

class ModALPrecondState;

class InvModALStrategy {
 public:
  //! Empty constructor.
  InvModALStrategy();

  InvModALStrategy(const Teuchos::RCP<InverseFactory>& factory);

  InvModALStrategy(const Teuchos::RCP<InverseFactory>& factory, LinearOp& pressureMassMatrix);

  InvModALStrategy(const Teuchos::RCP<InverseFactory>& invFactA,
                   const Teuchos::RCP<InverseFactory>& invFactS);

  InvModALStrategy(const Teuchos::RCP<InverseFactory>& invFactA,
                   const Teuchos::RCP<InverseFactory>& invFactS, LinearOp& pressureMassMatrix);

  //! Destructor.
  virtual ~InvModALStrategy() {}

  /** Get the inverse of the \f$A_{11}p = A_{11} + \gamma B^T_1 W^{-1} B_1 \f$ block.
   *
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$A_{11}p\f$.
   */
  virtual LinearOp getInvA11p(BlockPreconditionerState& state) const;

  /** Get the inverse of the \f$A_{22}p = A_{22} + \gamma B^T_2 W^{-1} B_2 \f$ block.
   *
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$A_{22}p\f$.
   */
  virtual LinearOp getInvA22p(BlockPreconditionerState& state) const;

  /** Get the inverse of the \f$A_{33}p = A_{33} + \gamma B^T_3 W^{-1} B_3 \f$ block.
   *
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$A_{33}p\f$.
   */
  virtual LinearOp getInvA33p(BlockPreconditionerState& state) const;

  /** Get the inverse of the pressure Schur complement \f$ S \f$.
   *
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$S\f$.
   */
  virtual LinearOp getInvS(BlockPreconditionerState& state) const;

  /** This informs the strategy object to build the state associated
   * with this operator.
   *
   * \param[in] A The linear operator to be preconditioned by modified AL.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   */
  virtual void buildState(const BlockedLinearOp& A, BlockPreconditionerState& state) const;

  /**
   * Initialize the state object using this blocked linear operator.
   */
  virtual void initializeState(const BlockedLinearOp& A, ModALPrecondState* state) const;

  /** Compute the inverses.
   *
   * \param[in] A ALOperator.
   *
   * \note This method assumes that the operators required have been constructed.
   */
  virtual void computeInverses(const BlockedLinearOp& A, ModALPrecondState* state) const;

  /** Set pressure mass matrix.
   *
   * \param[in] pressureMassMatrix
   *            Pressure mass matrix.
   */
  void setPressureMassMatrix(const LinearOp& pressureMassMatrix);

  /** Set the augmentation parameter gamma.
   *
   * \param[in] gamma
   *            Augmentation paramter.
   */
  void setGamma(double gamma);

  /** Tell strategy that this operator is supposed to be symmetric.
   *
   * \param[in] isSymmetric Is this operator symmetric?
   */
  virtual void setSymmetric(bool isSymmetric) { isSymmetric_ = isSymmetric; }

 protected:
  // In the modified AL preconditioner, we need to two methods,
  // one for solving \f$ A_{ii}, i = 1, 2(, 3) \f$,
  // the other for solving \f$ S \f$.
  Teuchos::RCP<InverseFactory> invFactoryA_;
  Teuchos::RCP<InverseFactory> invFactoryS_;
  LinearOp pressureMassMatrix_;
  double gamma_;

  DiagonalType scaleType_;
  bool isSymmetric_;
  int dim_;
};

}  // end namespace NS

}  // end namespace Teko

#endif /* __Teko_ModALStrategy_hpp__ */
