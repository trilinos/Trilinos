// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER


/** \file  Intrepid_ConfigDefs.hpp
    \brief Intrepid header file which uses auto-configuration information
           to include necessary C++ headers.
*/

#ifndef INTREPID_CONFIGDEFS_HPP
#define INTREPID_CONFIGDEFS_HPP

#ifndef TRILINOS_NO_CONFIG_H

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

#include "Intrepid_config.h"



#ifdef __cplusplus

/******************************************************************************
 *   Choose header file flavor: either ANSI-style (no .h, e.g. <iostream>) or
 *   old-style (with .h, e.g., <iostream.h>).
 *****************************************************************************/

#if HAVE_CSTDIO
#include <cstdio>
#elif HAVE_STDIO_H
#include <stdio.h>
#else
#error "Found neither cstdio nor stdio.h"
#endif

#if HAVE_CSTDLIB
#include <cstdlib>
#elif HAVE_STDLIB_H
#include <stdlib.h>
#else
#error "Found neither cstdlib nor stdlib.h"
#endif

#if HAVE_STRING
#include <string>
#elif HAVE_STRING_H
#include <string.h>
#else
#error "Found neither string nor string.h"
#endif

#if HAVE_IOSTREAM
#include <iostream>
#elif HAVE_IOSTREAM_H
#include <iostream.h>
#else
#error "Found neither iostream nor iostream.h"
#endif

#if HAVE_IOSTREAM
#include <fstream>
#elif HAVE_IOSTREAM_H
#include <fstream.h>
#else
#error "Found neither fstream nor fstream.h"
#endif

#if HAVE_SSTREAM
#include <sstream>
#elif HAVE_SSTREAM_H
#include <sstream.h>
#elif HAVE_STRSTREAM
#include <strstream>
#elif HAVE_STRSTREAM_H
#include <strstream.h>
#else
#error "Found neither sstream, sstream.h, strstream.h, nor strstream"
#endif

#ifdef HAVE_CASSERT
#include <cassert>
#else
#include <assert.h>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef HAVE_VECTOR
#include <vector>
#elif defined(HAVE_VECTOR_H)
#include <vector.h>
#endif

#ifdef HAVE_MAP
#include <map>
#elif defined(HAVE_MAP_H)
#include <map.h>
#endif

#ifdef HAVE_STRING
#include <string>
#elif defined(HAVE_STRING_H)
#include <string.h>
#elif defined(HAVE_STRINGS_H)
#include <strings.h>
#endif

#ifdef HAVE_IOMANIP
#include <iomanip>
#else
#include <iomanip.h>
#endif

#else /* __cplusplus not defined */

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#endif /* __cplusplus */

#endif /* end ndef TRILINOS_NO_CONFIG_H */

/*
 * Intrepid_Version() method
*/
namespace Intrepid {
  std::string Intrepid_Version();
}


#endif /* INTREPID_CONFIGDEFS_HPP */
