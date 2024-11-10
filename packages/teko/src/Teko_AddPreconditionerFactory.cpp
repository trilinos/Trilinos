// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_AddPreconditionerFactory.hpp"

namespace Teko {

using Teuchos::RCP;

AddPreconditionerFactory::AddPreconditionerFactory(
    const RCP<const BlockPreconditionerFactory> &FirstFactory,
    const RCP<const BlockPreconditionerFactory> &SecondFactory)
    : FirstFactory_(FirstFactory), SecondFactory_(SecondFactory) {}

AddPreconditionerFactory::AddPreconditionerFactory() {}

//! Build the AddPrecondState object
RCP<PreconditionerState> AddPreconditionerFactory::buildPreconditionerState() const {
  AddPrecondState *mystate = new AddPrecondState();
  mystate->StateOne_       = Teuchos::rcp_dynamic_cast<BlockPreconditionerState>(
      FirstFactory_->buildPreconditionerState());
  mystate->StateTwo_ = Teuchos::rcp_dynamic_cast<BlockPreconditionerState>(
      SecondFactory_->buildPreconditionerState());
  return rcp(mystate);
}

// Use the factory to build the preconditioner (this is where the work goes)
LinearOp AddPreconditionerFactory ::buildPreconditionerOperator(
    BlockedLinearOp &blockOp, BlockPreconditionerState &state) const {
  // The main tricky thing here is that we have to take the 'state' object
  // associated with AddPreconditionerFactory(), pull out the states for
  // the individual preconditioners, and pass these on to
  // buildPreconditionerOperator() for each subpreconditioner.

  AddPrecondState *MyState = dynamic_cast<AddPrecondState *>(&state);
  TEUCHOS_ASSERT(MyState != 0);

  LinearOp M1 = FirstFactory_->buildPreconditionerOperator(blockOp, *MyState->StateOne_);
  LinearOp M2 = SecondFactory_->buildPreconditionerOperator(blockOp, *MyState->StateTwo_);

  LinearOp invA = add(M1, M2);

  // return fully constructed preconditioner
  return invA;
}

//! Initialize from a parameter list
void AddPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList &pl) {
  RCP<const InverseLibrary> invLib = getInverseLibrary();

  // get string specifying inverse
  std::string aStr = "", bStr = "";

  // "parse" the parameter list
  aStr = pl.get<std::string>("Preconditioner A");
  bStr = pl.get<std::string>("Preconditioner B");

  RCP<const Teuchos::ParameterList> aSettings = invLib->getParameterList(aStr);
  RCP<const Teuchos::ParameterList> bSettings = invLib->getParameterList(bStr);

  // build preconditioner from the parameters
  std::string aType                      = aSettings->get<std::string>("Preconditioner Type");
  RCP<Teko::PreconditionerFactory> precA = Teko::PreconditionerFactory::buildPreconditionerFactory(
      aType, aSettings->sublist("Preconditioner Settings"), invLib);

  // build preconditioner from the parameters
  std::string bType                      = bSettings->get<std::string>("Preconditioner Type");
  RCP<Teko::PreconditionerFactory> precB = Teko::PreconditionerFactory::buildPreconditionerFactory(
      bType, bSettings->sublist("Preconditioner Settings"), invLib);

  // set preconditioners
  FirstFactory_  = Teuchos::rcp_dynamic_cast<const Teko::BlockPreconditionerFactory>(precA);
  SecondFactory_ = Teuchos::rcp_dynamic_cast<const Teko::BlockPreconditionerFactory>(precB);
}

}  // end namespace Teko
