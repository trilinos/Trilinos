// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziConfigDefs.hpp
  \brief Anasazi header file which uses auto-configuration information to include
  necessary C++ headers
*/

#ifndef ANASAZI_CONFIGDEFS_HPP
#define ANASAZI_CONFIGDEFS_HPP

#include "Teuchos_ConfigDefs.hpp"

#ifndef __cplusplus
#  define __cplusplus
#endif

#ifndef TRILINOS_NO_CONFIG_H

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 * KL 11/25/02
 */
#  ifdef PACKAGE
#    undef PACKAGE
#  endif

#  ifdef PACKAGE_NAME
#    undef PACKAGE_NAME
#  endif

#  ifdef PACKAGE_BUGREPORT
#    undef PACKAGE_BUGREPORT
#  endif

#  ifdef PACKAGE_STRING
#    undef PACKAGE_STRING
#  endif

#  ifdef PACKAGE_TARNAME
#    undef PACKAGE_TARNAME
#  endif

#  ifdef PACKAGE_VERSION
#    undef PACKAGE_VERSION
#  endif

#  ifdef VERSION
#    undef VERSION
#  endif

#  include <Anasazi_config.h>

#  ifdef HAVE_MPI
#    ifndef EPETRA_MPI
#      define EPETRA_MPI
#    endif
#  endif

#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <cctype>
#include <numeric>
#include <complex>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <cmath>
#include <functional>

#else /*TRILINOS_NO_CONFIG_H is defined*/

#  include <iterator>
#  include <iostream>
#  include <string>

#  if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(CPLANT) || defined (TFLOP)
#    include <stdlib.h>
#    include <stdio.h>
#    include <math.h>
#  else
#    include <cstdlib>
#    include <cstdio>
#    include <cmath>
#  endif

#  include <vector>
#  include <map>
#  include <deque>
#  include <algorithm>
#  include <numeric>
#  include <functional>

#endif /*ndef TRILINOS_NO_CONFIG_H*/

/* Define some macros */
#define ANASAZI_MAX(x,y) (( (x) > (y) ) ? (x)  : (y) )     /* max function  */
#define ANASAZI_MIN(x,y) (( (x) < (y) ) ? (x)  : (y) )     /* min function  */
#define ANASAZI_SGN(x)   (( (x) < 0.0 ) ? -1.0 : 1.0 )     /* sign function */

#include "Anasazi_DLLExportMacro.h"

/*
 * Anasazi_Version() method
 */
namespace Anasazi {
  ANASAZI_LIB_DLL_EXPORT std::string Anasazi_Version();
}

#endif /*ANASAZI_CONFIGDEFS_HPP*/
