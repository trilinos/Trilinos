/* $Id$ */
/* $Source$ */
/* 
@HEADER 
*************************************************************************

                         Sacado Package
               Copyright (2006) Sandia Corporation

Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
the U.S. Government retains certain rights in this software.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
(etphipp@sandia.gov).

************************************************************************
@HEADER 
*/

#ifndef SACADO_CONFIGDEFS_H
#define SACADO_CONFIGDEFS_H

#ifndef __cplusplus
#define __cplusplus
#endif

/*
 * The macros PACKAGE, PACKAGE_NAME, etc, get defined for each package and 
 * need to be undef'd here to avoid warnings when this file is included from 
 * another package.
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

#include <Sacado_config.h>

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

#ifdef HAVE_CMATH
#include <cmath>
#elif defined(HAVE_MATH_H)
#include <math.h>
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

#ifdef HAVE_IOMANIP
#include <iomanip>
#elif defined(HAVE_IOMANIP_H)
#include <iomanip.h>
#endif

#ifdef HAVE_ALGORITHM
#include <algorithm>
#elif defined(HAVE_ALGORITHM_H)
#include <algorithm.h>
#elif defined(HAVE_ALGO_H)
#include <algo.h>
#endif

#ifdef HAVE_NEW
#include <new>
#elif defined(HAVE_NEW_H)
#include <new.h>
#endif

#ifdef HAVE_VECTOR
#include <vector>
#elif defined(HAVE_VECTOR_H)
#include <vector.h>
#endif

#ifdef HAVE_VALARRAY
#include <valarray>
#elif defined(HAVE_VALARRAY_H)
#include <valarray.h>
#endif

#ifdef HAVE_MAP
#include <map>
#elif defined(HAVE_MAP_H)
#include <map.h>
#endif

#ifdef HAVE_ITERATOR
#include <iterator>
#elif defined(HAVE_ITERATOR_H)
#include <iterator.h>
#endif

#ifdef HAVE_TYPEINFO
#include <typeinfo>
#elif defined(HAVE_TYPEINFO_H)
#include <typeinfo.h>
#endif

#endif /* SACADO_CONFIGDEFS_H */
