// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Teko includes
#include "Teko_DiagonallyScaledPreconditionerFactory.hpp"

#include "Teko_PreconditionerInverseFactory.hpp"

namespace Teko {

//! Default constructor, for use with the AutoClone class.
DiagonallyScaledPreconditionerFactory::DiagonallyScaledPreconditionerFactory()
    : invFactory_(Teuchos::null), scalingType_(ROW_SCALING), diagonalType_(AbsRowSum) {}

/** Construct a preconditioner factory that applies a specified
 * preconditioner, a fixed number of times.
 */
DiagonallyScaledPreconditionerFactory::DiagonallyScaledPreconditionerFactory(
    const Teuchos::RCP<Teko::InverseFactory>& invFactory, ScalingType scalingType,
    DiagonalType diagonalType)
    : invFactory_(invFactory), scalingType_(scalingType), diagonalType_(diagonalType) {}

DiagonallyScaledPreconditionerFactory::~DiagonallyScaledPreconditionerFactory() {}

/** \brief Function that is called to build the preconditioner
 *        for the linear operator that is passed in.
 */
LinearOp DiagonallyScaledPreconditionerFactory::buildPreconditionerOperator(
    LinearOp& lo, PreconditionerState& state) const {
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  Teko_DEBUG_SCOPE("DiagonallyScaledPreconditionerFactory::buildPreconditionerOperator", 10);

  TEUCHOS_TEST_FOR_EXCEPTION(
      invFactory_ == Teuchos::null, std::runtime_error,
      "ERROR: Teko::DiagonallyScaledPreconditionerFactory::buildPreconditionerOperator requires "
      "that an "
          << "inverse factory has been set. Currently it is null!");

  // get diagonal matrix
  LinearOp invD = getInvDiagonalOp(lo, diagonalType_);

  // M = A * invD
  ModifiableLinearOp& M = state.getModifiableOp("op_M");
  if (scalingType_ == COLUMN_SCALING)
    M = explicitMultiply(lo, invD, M);
  else
    M = explicitMultiply(invD, lo, M);

  // build inverse operator
  ModifiableLinearOp& invM = state.getModifiableOp("op_invM");
  if (invM == Teuchos::null)
    invM = buildInverse(*invFactory_, M);
  else
    rebuildInverse(*invFactory_, M, invM);

  // return invD * invM
  if (scalingType_ == COLUMN_SCALING)
    return multiply(invD, invM.getConst());
  else
    return multiply(invM.getConst(), invD);
}

/** \brief This function builds the internals of the preconditioner factory
 *        from a parameter list.
 */
void DiagonallyScaledPreconditionerFactory::initializeFromParameterList(
    const Teuchos::ParameterList& settings) {
  TEUCHOS_TEST_FOR_EXCEPTION(
      not settings.isParameter("Inverse Factory"), std::runtime_error,
      "Parameter \"Inverse Factory\" is required by a Teko::DiagonallyScaledPreconditionerFactory");

  // grab library and preconditioner name
  std::string invName = settings.get<std::string>("Inverse Factory");

  // build preconditioner factory
  Teuchos::RCP<const InverseLibrary> il = getInverseLibrary();
  invFactory_                           = il->getInverseFactory(invName);
  TEUCHOS_TEST_FOR_EXCEPTION(invFactory_ == Teuchos::null, std::runtime_error,
                             "ERROR: \"Inverse Factory\" = " << invName << " could not be found");

  // get scaling type specified by XML file
  const std::string defaultScaleType  = "Row";
  const std::string scalingTypeString = "Scaling Type";
  std::string scaleType               = defaultScaleType;
  if (settings.isParameter(scalingTypeString))
    scaleType = settings.get<std::string>(scalingTypeString);

  if (defaultScaleType == scaleType)
    scalingType_ = ROW_SCALING;
  else
    scalingType_ = COLUMN_SCALING;

  if (settings.isParameter("Diagonal Type"))
    diagonalType_ = Teko::getDiagonalType(settings.get<std::string>("Diagonal Type"));
}

}  // end namespace Teko
