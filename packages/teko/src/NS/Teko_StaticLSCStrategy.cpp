// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NS/Teko_StaticLSCStrategy.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"

#include "Teuchos_Time.hpp"

// Teko includes
#include "Teko_Utilities.hpp"

using Teuchos::RCP;
using Teuchos::rcp_const_cast;
using Teuchos::rcp_dynamic_cast;

namespace Teko {
namespace NS {

// Staiblized constructor
StaticLSCStrategy::StaticLSCStrategy(const LinearOp& invF, const LinearOp& invBQBtmC,
                                     const LinearOp& invD, const LinearOp& invMass)
    : invF_(invF), invBQBtmC_(invBQBtmC), invD_(invD), invMass_(invMass) {}

// Stable constructor
StaticLSCStrategy::StaticLSCStrategy(const LinearOp& invF, const LinearOp& invBQBtmC,
                                     const LinearOp& invMass)
    : invF_(invF), invBQBtmC_(invBQBtmC), invD_(Teuchos::null), invMass_(invMass) {}

}  // end namespace NS
}  // end namespace Teko
