// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_StaticLSCStrategy_hpp__
#define __Teko_StaticLSCStrategy_hpp__

#include "Teko_LSCStrategy.hpp"

namespace Teko {
namespace NS {

class LSCPrecondState;  // forward declaration

// constant, not very flexible strategy for driving LSCPreconditioenrFactory
class StaticLSCStrategy : public LSCStrategy {
 public:
  // Staiblized constructor
  StaticLSCStrategy(const LinearOp& invF, const LinearOp& invBQBtmC, const LinearOp& invD,
                    const LinearOp& invMass);

  // Stable constructor
  StaticLSCStrategy(const LinearOp& invF, const LinearOp& invBQBtmC, const LinearOp& invMass);

  /** This informs the strategy object to build the state associated
   * with this operator.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   */
  virtual void buildState(BlockedLinearOp& /* A */, BlockPreconditionerState& /* state */) const {}

  /** Get the inverse of the \f$F\f$ block.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$F\f$.
   */
  virtual LinearOp getInvF(const BlockedLinearOp& /* A */,
                           BlockPreconditionerState& /* state */) const {
    return invF_;
  }

  /** Get the inverse of \f$B Q_u^{-1} B^T\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$B Q_u^{-1} B^T\f$.
   */
  virtual LinearOp getInvBQBt(const BlockedLinearOp& /* A */,
                              BlockPreconditionerState& /* state */) const {
    return invBQBtmC_;
  }

  /** Get the inverse of \f$B H B^T - \gamma C\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$B H B^T - \gamma C\f$.
   */
  virtual LinearOp getInvBHBt(const BlockedLinearOp& /* A */,
                              BlockPreconditionerState& /* state */) const {
    return invBQBtmC_;
  }

  /** Get the inverse for stabilizing the whole Schur complement approximation.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The operator to stabilize the whole Schur complement (\f$\alpha D^{-1} \f$).
   */
  virtual LinearOp getOuterStabilization(const BlockedLinearOp& /* A */,
                                         BlockPreconditionerState& /* state */) const {
    return invD_;
  }

  virtual LinearOp getInnerStabilization(const BlockedLinearOp& /* A */,
                                         BlockPreconditionerState& /* state */) const {
    return Teuchos::null;
  }

  /** Get the inverse mass matrix.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The inverse of the mass matrix \f$Q_u\f$.
   */
  virtual LinearOp getInvMass(const BlockedLinearOp& /* A */,
                              BlockPreconditionerState& /* state */) const {
    return invMass_;
  }

  /** Get the \f$H\f$ scaling matrix.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The \f$H\f$ scaling matrix.
   */
  virtual LinearOp getHScaling(const BlockedLinearOp& /* A */,
                               BlockPreconditionerState& /* state */) const {
    return invMass_;
  }

  /** Should the approximation of the inverse use a full LDU decomposition, or
   * is a upper triangular approximation sufficient.
   *
   * \returns True if the full LDU decomposition should be used, otherwise
   *          only an upper triangular version is used.
   */
  virtual bool useFullLDU() const { return false; }

  /** Tell strategy that this operator is supposed to be symmetric.
   * Behavior of LSC is slightly different for non-symmetric case.
   *
   * \param[in] isSymmetric Is this operator symmetric?
   */
  virtual void setSymmetric(bool /* isSymmetric */) {}

 protected:
  // protected memebers
  LinearOp invF_;
  LinearOp invBQBtmC_;
  LinearOp invD_;
  LinearOp invMass_;
};

}  // end namespace NS
}  // end namespace Teko

#endif
