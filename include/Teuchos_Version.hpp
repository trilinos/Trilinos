// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_VERSION_HPP
#define TEUCHOS_VERSION_HPP

#include "Teuchos_ConfigDefs.hpp"
#ifdef TEUCHOS_STANDALONE_PACKAGE
#  include "Teuchos_version.h"
#else
#  include "Trilinos_version.h"
#endif

namespace Teuchos {

	std::string Teuchos_Version() {
#ifdef TEUCHOS_STANDALONE_PACKAGE
		return ("Teuchos standalone package " TEUCHOS_VERSION_STRING);
#else
		return ("Teuchos in Trilinos " TRILINOS_VERSION_STRING);
#endif
	}

} // namespace Teuchos

#endif // TEUCHOS_VERSION_HPP

