// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


/** \file  Intrepid_ConfigDefs.hpp
    \brief Intrepid header file which uses auto-configuration information
           to include necessary C++ headers.
*/
#ifndef INTREPID2_CONFIGDEFS_HPP
#define INTREPID2_CONFIGDEFS_HPP

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 */

#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include "Intrepid2_config.h"

/******************************************************************************
 *   Choose header file flavor: either ANSI-style (no .h, e.g. <iostream>) or
 *   old-style (with .h, e.g., <iostream.h>).
 *****************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <iomanip>

/*
 * Intrepid_Version() method
*/
namespace Intrepid2 {
  std::string Intrepid_Version();
}


#endif /* INTREPID_CONFIGDEFS_HPP */
