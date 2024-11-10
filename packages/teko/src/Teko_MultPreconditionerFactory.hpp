// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_MultiPreconditionerFactory_hpp__
#define __Teko_MultiPreconditionerFactory_hpp__

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"
#include "Teko_BlockImplicitLinearOp.hpp"

namespace Teko {

/** Preconditioning factories must supply a 'State' class which
 * is where data specific to the preconditioner construction
 * is stored. The constructor will be invoked elsewhere.
 */
class MultPrecondState : public Teko::BlockPreconditionerState {
 public:
  MultPrecondState() {}

  Teuchos::RCP<BlockPreconditionerState> StateOne_;
  Teuchos::RCP<BlockPreconditionerState> StateTwo_;
};

/** Declaration of a Teko::BlockImplicitLinearOp.
 * BlockImplicitLinearOp's are useful
 * when there is no simple or cheap matrix representation
 * of something like f(r). This particular class
 * corresponds to
 *         f(r) = (M2 + M1 - M2 * A * M1) r
 *
 * which is an application of a multiplicative preconditioner.
 * Notice that when M1 = inv(A) or when M2 = inv(A), we get
 * f(r) = inv(A)*r.
 *
 * While the above matrix represenation could be used
 * instead of writing an implicit function, it requires
 * an additional matrix-vector product than an efficient
 * implementation. It should be noted (see comments below)
 * that there is an efficient matrix represenation of this
 * f(r). Namely,
 *
 *             f(r) =  [M2 I] [I  -A] [ I]
 *                            [0   I] [ M1]
 *
 * so we didn't really need to create an implicit operator.
 */
class MultPrecsLinearOp : public Teko::BlockImplicitLinearOp {
 public:
  //! Constructor
  MultPrecsLinearOp(const Teko::LinearOp &A, const Teko::LinearOp &M1, const Teko::LinearOp &M2)
      : A_(A), M1_(M1), M2_(M2) {}

  virtual Teko::VectorSpace range() const { return M1_->range(); }
  virtual Teko::VectorSpace domain() const { return M1_->domain(); }
  virtual void implicitApply(const Teko::BlockedMultiVector &r, Teko::BlockedMultiVector &y,
                             const double alpha = 1.0, const double beta = 0.0) const;

 protected:
  using Teko::BlockImplicitLinearOp::implicitApply;

  Teko::LinearOp A_, M1_, M2_;

 private:
  // hide me!
  MultPrecsLinearOp();
  MultPrecsLinearOp(const MultPrecsLinearOp &);
};

/** Declaration of preconditioner factory that creates
 * a preconditioner which is the multiplicative combination
 * of two other preconditioners.
 */
class MultPreconditionerFactory : public Teko::BlockPreconditionerFactory {
 public:
  //! Constructor
  MultPreconditionerFactory(
      const Teuchos::RCP<const Teko::BlockPreconditionerFactory> &FirstFactory,
      const Teuchos::RCP<const Teko::BlockPreconditionerFactory> &SecondFactory);

  MultPreconditionerFactory();

  //! Function inherited from Teko::BlockPreconditionerFactory
  Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp &blo,
                                             Teko::BlockPreconditionerState &state) const;

  //! Build the MultPrecondState object
  virtual Teuchos::RCP<Teko::PreconditionerState> buildPreconditionerState() const;

 protected:
  using Teko::BlockPreconditionerFactory::buildPreconditionerOperator;

  // class members
  Teuchos::RCP<const Teko::BlockPreconditionerFactory> FirstFactory_;
  Teuchos::RCP<const Teko::BlockPreconditionerFactory> SecondFactory_;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList &pl);
};

}  // end namespace Teko

#endif
