// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_VALIDPIROPARAMETERS
#define PIRO_VALIDPIROPARAMETERS

#include "Teuchos_ParameterList.hpp"

namespace Piro {

Teuchos::RCP<const Teuchos::ParameterList>
getValidPiroParameters();

}

#endif //PIRO_VALIDPIROPARAMETERS
