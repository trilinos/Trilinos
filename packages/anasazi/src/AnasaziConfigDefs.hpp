// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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
