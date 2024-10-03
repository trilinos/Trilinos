// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziVersion.cpp
    \brief Version function that will return the current version of Anasazi being utilized
*/

#include "AnasaziConfigDefs.hpp"
#include "Trilinos_version.h"

namespace Anasazi {

   std::string Anasazi_Version() { 
		return("Anasazi in Trilinos " TRILINOS_VERSION_STRING); 
	}

} // namespace Anasazi
