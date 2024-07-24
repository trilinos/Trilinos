// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_MultPreconditionerFactory.hpp"

namespace Teko {

using Teuchos::RCP;

void MultPrecsLinearOp::implicitApply(const Teko::BlockedMultiVector &r,
                                      Teko::BlockedMultiVector &y, const double /* alpha */,
                                      const double /* beta */) const {
  // Casting is a bit delicate. We basically use
  //
  //  1) deepcopy      to copy & cast BlockedMultiVectors to MultiVectors.
  //
  //  2) toMultiVector to cast BlockedMultiVectors to MultiVectors.
  //
  Teko::MultiVector MOne_r = Teko::deepcopy(r);
  Teko::MultiVector t      = Teko::deepcopy(r);
  Teko::MultiVector w      = Teko::toMultiVector(y);

  Teko::applyOp(M1_, r, MOne_r);
  Teko::applyOp(A_, MOne_r, t);
  Teko::update(1., r, -1., t);
  Teko::applyOp(M2_, t, w);
  Teko::update(1., MOne_r, 1., w);
}

//! Constructor
MultPreconditionerFactory ::MultPreconditionerFactory(
    const RCP<const Teko::BlockPreconditionerFactory> &FirstFactory,
    const RCP<const Teko::BlockPreconditionerFactory> &SecondFactory)
    : FirstFactory_(FirstFactory), SecondFactory_(SecondFactory) {}

MultPreconditionerFactory::MultPreconditionerFactory() {}

//! Build the MultPrecondState object
RCP<Teko::PreconditionerState> MultPreconditionerFactory::buildPreconditionerState() const {
  MultPrecondState *mystate = new MultPrecondState();
  mystate->StateOne_        = Teuchos::rcp_dynamic_cast<BlockPreconditionerState>(
      FirstFactory_->buildPreconditionerState());
  mystate->StateTwo_ = Teuchos::rcp_dynamic_cast<BlockPreconditionerState>(
      SecondFactory_->buildPreconditionerState());
  return rcp(mystate);
}

//! Use the factory to build the preconditioner (this is where the work goes)
Teko::LinearOp MultPreconditionerFactory ::buildPreconditionerOperator(
    Teko::BlockedLinearOp &blockOp, Teko::BlockPreconditionerState &state) const {
  MultPrecondState *MyState = dynamic_cast<MultPrecondState *>(&state);

  TEUCHOS_ASSERT(MyState != 0);

  Teko::LinearOp M1 = FirstFactory_->buildPreconditionerOperator(blockOp, *MyState->StateOne_);
  Teko::LinearOp M2 = SecondFactory_->buildPreconditionerOperator(blockOp, *MyState->StateTwo_);

  /*************************************************************************
     A different way to create the same preconditioner using the funky
     matrix representation discussed above. At the present time, there
     appears to be some kind of bug in Thrya so this doesn't work.

  const RCP<const Thyra::LinearOpBase<double>> Mat1=
  Thyra::block2x1(Teko::identity(Teko::rangeSpace(M1)) ,M1); const RCP<const
  Thyra::LinearOpBase<double>> Mat3= Thyra::block1x2(M2,Teko::identity(Teko::rangeSpace(M1))); const
  RCP<const Thyra::LinearOpBase<double>> Mat2= Thyra::block2x2(
                    Teko::identity(Teko::rangeSpace(M1)),
  Teko::scale(-1.,Teko::toLinearOp(blockOp)),
                    Thyra::zero<double>(Teko::rangeSpace(M1),Teko::domainSpace(M1)),
  Teko::identity(Teko::rangeSpace(M1))); Teko::LinearOp invA = Teko::multiply(Mat3,Mat2,Mat1);

  return invA;
   *************************************************************************/

  // construct an implicit operator corresponding to multiplicative
  // preconditioning, wrap it in an rcp pointer and return.

  return Teuchos::rcp(new MultPrecsLinearOp(blockOp, M1, M2));
}

//! Initialize from a parameter list
void MultPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList &pl) {
  RCP<const InverseLibrary> invLib = getInverseLibrary();

  // get string specifying inverse
  std::string aStr = "", bStr = "";

  // "parse" the parameter list
  aStr = pl.get<std::string>("Preconditioner A");
  bStr = pl.get<std::string>("Preconditioner B");

  RCP<const Teuchos::ParameterList> aSettings = invLib->getParameterList(aStr);
  RCP<const Teuchos::ParameterList> bSettings = invLib->getParameterList(bStr);

  // build preconditioner from the parameters
  std::string aType = aSettings->get<std::string>("Preconditioner Type");
  RCP<Teko::PreconditionerFactory> precA =
      Teko::BlockPreconditionerFactory::buildPreconditionerFactory(
          aType, aSettings->sublist("Preconditioner Settings"), invLib);

  // build preconditioner from the parameters
  std::string bType = bSettings->get<std::string>("Preconditioner Type");
  RCP<Teko::PreconditionerFactory> precB =
      Teko::BlockPreconditionerFactory::buildPreconditionerFactory(
          bType, bSettings->sublist("Preconditioner Settings"), invLib);

  // set preconditioners
  FirstFactory_  = Teuchos::rcp_dynamic_cast<const Teko::BlockPreconditionerFactory>(precA);
  SecondFactory_ = Teuchos::rcp_dynamic_cast<const Teko::BlockPreconditionerFactory>(precB);
}

}  // end namespace Teko
