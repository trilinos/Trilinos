
//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_CONFIGDEFS_H
#define EPETRAEXT_CONFIGDEFS_H

#ifndef __cplusplus
#define __cplusplus
#endif

#ifdef HAVE_CONFIG_H

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and need to
 * be undef'd here to avoid warnings when this file is included from another package.
 * KL 11/25/02
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

#include <EpetraExt_config.h>

#ifdef HAVE_CSTDLIB
#include <cstdlib>
#elif defined(HAVE_STDLIB_H)
#include <stdlib.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#elif defined(HAVE_STDIO_H)
#include <stdio.h>
#endif

#ifdef HAVE_CASSERT
#include <cassert>
#elif defined(HAVE_ASSERT_H)
#include <assert.h>
#endif

#ifdef HAVE_STRING
#include <string>
#elif defined(HAVE_STRING_H)
#include <string.h>
#endif

#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#endif

#ifdef HAVE_ALGORITHM
#include <algorithm>
#elif defined(HAVE_ALGORITHM_H)
#include <algorithm.h>
#endif

/* Every line that begins with 'using' should eventually be dependent
   on some check within the configure script */

#ifndef TFLOP

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

using namespace std;

#else /* TFLOP defined */

#ifdef HAVE_IOMANIP
#include <iomanip>
#else
#include <iomanip.h>
#endif

#ifdef HAVE_STRING
using std::string;
#endif

#ifdef HAVE_IOSTREAM
using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;
#endif

#endif

/*-----------------------------------------------------------------------
  Must refine the following up to #else HAVE_CONFIG_H is not defined
  -----------------------------------------------------------------------*/

#else /*HAVE_CONFIG_H is not defined*/

#ifndef __cplusplus
#define __cplusplus
#endif

#if defined(SGI) || defined(SGI64) || defined(SGI32) || defined(CPLANT)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <math.h>
#include <string>
using namespace std;

#elif defined(TFLOP)

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string>
using std::string;
#include <iostream>
#include <iomanip>
using std::istream;
using std::ostream;
using std::cerr;
using std::cout;
using std::endl;

#else

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <cmath>
#include <string>
using namespace std;

#endif

#endif

#endif /* EPETRAEXT_CONFIGDEFS_H */
