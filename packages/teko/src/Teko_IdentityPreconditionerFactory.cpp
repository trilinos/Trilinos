// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_IdentityPreconditionerFactory.hpp"

#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"

using Teuchos::rcp;

namespace Teko {

/** Build a Identity preconditioner factory from a parameter list
 */
IdentityPreconditionerFactory::IdentityPreconditionerFactory() : scaling_(1.0) {}

LinearOp IdentityPreconditionerFactory::buildPreconditionerOperator(
    LinearOp& lo, PreconditionerState& /* state */) const {
  return Thyra::scale(scaling_, Thyra::identity(rangeSpace(lo)));
}

//! Initialize from a parameter list
void IdentityPreconditionerFactory::initializeFromParameterList(const Teuchos::ParameterList& pl) {
  Teko_DEBUG_SCOPE("IdentityPreconditionerFactory::initializeFromParameterList", 10);
  Teko_DEBUG_MSG_BEGIN(9);
  DEBUG_STREAM << "Parameter list: " << std::endl;
  pl.print(DEBUG_STREAM);
  Teko_DEBUG_MSG_END();

  // get string specifying default inverse
  std::string scaleStr = "Scaling";
  if (pl.isParameter(scaleStr)) scaling_ = pl.get<double>(scaleStr);
}

}  // namespace Teko
