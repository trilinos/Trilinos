// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_LSCSIMPLECStrategy_hpp__
#define __Teko_LSCSIMPLECStrategy_hpp__

#include "Teko_LSCStrategy.hpp"

namespace Teko {
namespace NS {

class LSCPrecondState;  // forward declration

/** \brief A strategy that takes a single inverse factory and
 *        uses that for all inverses. If no mass matrix is
 *        passed in the diagonal of the 1,1 block is used.
 *
 * A strategy that takes a single inverse factory and uses that
 * for all inverses. Optionally the mass matrix can be passed
 * in, if it is the diagonal is extracted and that is used to
 * form the inverse approximation.
 */
class LSCSIMPLECStrategy : public LSCStrategy {
 public:
  //! \name Constructors
  //@{
  LSCSIMPLECStrategy();
  //@}

  virtual ~LSCSIMPLECStrategy() {}

  //! Functions inherited from LSCStrategy
  //@{

  /** This informs the strategy object to build the state associated
   * with this operator.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   */
  virtual void buildState(BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the inverse of \f$B Q_u^{-1} B^T\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$B Q_u^{-1} B^T\f$.
   */
  virtual LinearOp getInvBQBt(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the inverse of \f$B H B^T - \gamma C\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$B H B^T - \gamma C\f$.
   */
  virtual LinearOp getInvBHBt(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the inverse of the \f$F\f$ block.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns An (approximate) inverse of \f$F\f$.
   */
  virtual LinearOp getInvF(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the inverse for stabilizing the whole schur complement approximation.
   * Returns \f$\alpha D^{-1}\f$.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The operator to stabilize the whole Schur complement.
   */
  virtual LinearOp getOuterStabilization(const BlockedLinearOp & /* A */,
                                         BlockPreconditionerState & /* state */) const {
    return Teuchos::null;
  }

  virtual LinearOp getInnerStabilization(const BlockedLinearOp &A,
                                         BlockPreconditionerState &state) const;

  /** Get the inverse mass matrix.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The inverse of the mass matrix \f$Q_u\f$.
   */
  virtual LinearOp getInvMass(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Get the \f$H\f$ scaling matrix.
   *
   * \param[in] A The linear operator to be preconditioned by LSC.
   * \param[in] state State object for storying reusable information about
   *                  the operator A.
   *
   * \returns The \f$H\f$ scaling matrix.
   */
  virtual LinearOp getHScaling(const BlockedLinearOp &A, BlockPreconditionerState &state) const;

  /** Should the approximation of the inverse use a full LDU decomposition, or
   * is a upper triangular approximation sufficient.
   *
   * \returns True if the full LDU decomposition should be used, otherwise
   *          only an upper triangular version is used.
   */
  virtual bool useFullLDU() const { return useFullLDU_; }

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList &pl,
                                           const InverseLibrary &invLib);

  //! Initialize the state object using this blocked linear operator
  virtual void initializeState(const BlockedLinearOp &A, LSCPrecondState *state) const;

  /** Compute the inverses required for the LSC Schur complement
   *
   * \note This method assumes that the BQBt and BHBt operators have
   *       been constructed.
   */
  void computeInverses(const BlockedLinearOp &A, LSCPrecondState *state) const;

  //! Set to true to use the Full LDU decomposition, false otherwise
  virtual void setUseFullLDU(bool val) { useFullLDU_ = val; }

  virtual void setSymmetric(bool /* isSymmetric */) {}

 protected:
  // how to invert the matrices
  Teuchos::RCP<InverseFactory> invFactoryF_;
  Teuchos::RCP<InverseFactory> invFactoryS_;

  // flags for handling various options
  bool useFullLDU_;
  DiagonalType scaleType_;

 private:
  LSCSIMPLECStrategy(const LSCSIMPLECStrategy &);
};

}  // end namespace NS
}  // end namespace Teko

#endif
