// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_LU2x2DiagonalStrategy_hpp__
#define __Teko_LU2x2DiagonalStrategy_hpp__

#include "Teko_LU2x2Strategy.hpp"

// Teuchos includes
#include "Teuchos_Time.hpp"

namespace Teko {

/** @brief Strategy for computing \f$A_00^{-1}\f$ and \f$S^{-1}\f$ in the
 *         LU2x2PreconditionerFactory. Uses the diagonal of \f$A_00\f$ to
 *         build \f$S\f$.
 *
 * This should be paired with a LU2x2PreconditionerFactory,
 * it build the \f$A_{00}^{-1}\f$ and \f$S^{-1}\f$ operators.
 * where the Shur complement is
 * \f$ S = -A_{11} + A_{10} \mbox{diag}(A_{00})^{-1} A_{01} \f$.
 */
class LU2x2DiagonalStrategy : public LU2x2Strategy {
 public:
  //! default Constructor
  LU2x2DiagonalStrategy();

  //! Constructor to set the inverse factories.
  LU2x2DiagonalStrategy(const Teuchos::RCP<InverseFactory>& invFA,
                        const Teuchos::RCP<InverseFactory>& invS);

  //! Destructor (does nothing)
  virtual ~LU2x2DiagonalStrategy() {}

  /** returns the first (approximate) inverse of \f$A_{00}\f$ */
  virtual const Teko::LinearOp getHatInvA00(const Teko::BlockedLinearOp& A,
                                            BlockPreconditionerState& state) const;

  /** returns the second (approximate) inverse of \f$A_{00}\f$ */
  virtual const Teko::LinearOp getTildeInvA00(const Teko::BlockedLinearOp& A,
                                              BlockPreconditionerState& state) const;

  /** returns an (approximate) inverse of \f$S = -A_{11} + A_{10} \mbox{diag}(A_{00})^{-1} A_{01}\f$
   */
  virtual const Teko::LinearOp getInvS(const Teko::BlockedLinearOp& A,
                                       BlockPreconditionerState& state) const;

  /** \brief This function builds the internals of the state from a parameter list.
   *
   * This function builds the internals of the LU 2x2 state
   * from a parameter list. Furthermore, it allows a
   * developer to easily add a factory to the build system.
   *
   * \param[in] settings Parameter list to use as the internal settings
   * \param[in] invLib Inverse library to use for building inverse factory objects
   *
   * \note The default implementation does nothing.
   */
  virtual void initializeFromParameterList(const Teuchos::ParameterList& settings,
                                           const InverseLibrary& invLib);

 protected:
  /** Build timers for this type of object.
   */
  static void buildTimers();

  /** Initialize the operator's state. This builds the Schur complement and the inverse
   * operators. If the state has already been initialized this method does nothing.
   *
   * \param[in] A Operator to intialize with.
   * \param[in] state Storage object for this operator.
   */
  void initializeState(const Teko::BlockedLinearOp& A, BlockPreconditionerState& state) const;

  // how to invert the matrices
  Teuchos::RCP<InverseFactory> invFactoryA00_;  // for \tilde{A_00}\f$
  Teuchos::RCP<InverseFactory> invFactoryS_;

  DiagonalType a00InverseType_;

  static Teuchos::RCP<Teuchos::Time> initTimer_;
  static Teuchos::RCP<Teuchos::Time> invSTimer_;
  static Teuchos::RCP<Teuchos::Time> invA00Timer_;
  static Teuchos::RCP<Teuchos::Time> opsTimer_;
};

}  // end namespace Teko

#endif
