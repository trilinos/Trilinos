// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOSCORE_CONFIGDEFS_HPP
#define TEUCHOSCORE_CONFIGDEFS_HPP

/*! \file Teuchos_ConfigDefs.hpp
    \brief Teuchos header file which uses auto-configuration information
	to include necessary C++ headers.
*/

#include "TeuchosCore_config.h"

#ifdef HAVE_TEUCHOSCORE_CXX11
#  define TEUCHOS_NOEXCEPT_FALSE noexcept(false)
#else
#  define TEUCHOS_NOEXCEPT_FALSE
#endif

#endif /* TEUCHOSCORE_CONFIGDEFS_HPP */
