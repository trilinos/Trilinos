// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "BelosConfigDefs.hpp"

/*! \file BelosVersion.cpp
    \brief Simple function for returning the current version number [necessary for portability]
*/

namespace Belos {

	std::string Belos_Version() { 
		return("Belos Version 1.3d - 9/17/2008"); 
	}

} // namespace Belos 
