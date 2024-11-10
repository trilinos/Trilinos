// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_AddPreconditionerFactory_hpp__
#define __Teko_AddPreconditionerFactory_hpp__

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_Utilities.hpp"

namespace Teko {

/** Preconditioning factories must supply a 'State' class which
 * is where data specific to the preconditioner construction
 * is stored. The constructor will be invoked elsewhere.
 */
class AddPrecondState : public Teko::BlockPreconditionerState {
 public:
  AddPrecondState() {}

  Teuchos::RCP<BlockPreconditionerState> StateOne_;
  Teuchos::RCP<BlockPreconditionerState> StateTwo_;
};

/** Declaration of preconditioner factory that creates
 * a preconditioner which is the sum (additive) of two
 * other preconditioners.
 */
class AddPreconditionerFactory : public Teko::BlockPreconditionerFactory {
 public:
  //! Constructor
  AddPreconditionerFactory(
      const Teuchos::RCP<const Teko::BlockPreconditionerFactory>& FirstFactory,
      const Teuchos::RCP<const Teko::BlockPreconditionerFactory>& SecondFactory);

  AddPreconditionerFactory();

  //! Function inherited from Teko::BlockPreconditionerFactory
  Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp& blo,
                                             Teko::BlockPreconditionerState& state) const;

  //! Build the AddPrecondState object
  virtual Teuchos::RCP<Teko::PreconditionerState> buildPreconditionerState() const;

 protected:
  using Teko::BlockPreconditionerFactory::buildPreconditionerOperator;

  // class members
  Teuchos::RCP<const Teko::BlockPreconditionerFactory> FirstFactory_;
  Teuchos::RCP<const Teko::BlockPreconditionerFactory> SecondFactory_;

  //! Initialize from a parameter list
  virtual void initializeFromParameterList(const Teuchos::ParameterList& pl);
};

}  // end namespace Teko

#endif
