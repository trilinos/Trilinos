// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _THYRA_CONFIGDEFS_H_
#define _THYRA_CONFIGDEFS_H_

// Let Teuchos' configure process do all of the work!
#include <Teuchos_ConfigDefs.hpp>

#include <Thyra_Config.h>

// We need these macros in a lot of files and this is small so let's include
// it here.
#include "Teuchos_CompilerCodeTweakMacros.hpp"

#if defined(HAVE_THYRA_DEBUG) && !defined(THYRA_DEBUG)
#  define THYRA_DEBUG
#endif

#endif /*_THYRA_CONFIGDEFS_H_*/
