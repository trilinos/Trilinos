// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_COMMON_H
#define NOX_COMMON_H

// The following undef's are to avoid conflicts with other packages

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

#ifdef HAVE_CONFIG_H
#include "NOX_Config.h"
#endif

#ifdef HAVE_CSTDLIB
#include <cstdlib>
#elif defined(HAVE_STDLIB_H)
#include <stdlib.h>
#else
// Try to include it anyway
#include <stdlib.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#elif defined(HAVE_STDIO_H)
#include <stdio.h>
#else
// Try to include it anyway
#include <stdio.h>
#endif

#ifdef HAVE_CASSERT
#include <cassert>
#elif defined(HAVE_ASSERT_H)
#include <assert.h>
#else
// Try to include it anyway
#include <assert.h>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#elif defined(HAVE_MATH_H)
#include <math.h>
#else
// Try to include it anyway
#include <math.h>
#endif

#ifdef HAVE_IOMANIP
#include <iomanip>
#elif defined(HAVE_IOMANIP_H)
#include <iomanip.h>
#else
// Try to include it anyway
#include <iomanip>
#endif

#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#else
// Try to include it anyway
#include <iostream>
#endif

#ifdef HAVE_STRING
#include <string>
#elif defined(HAVE_STRING_H)
#include <string.h>
#else
// Try to include it anyway
#include <string>
#endif

#ifdef HAVE_MAP
#include <map>
#elif defined(HAVE_MAP_H)
#include <map.h>
#else
// Try to include it anyway
#include <map>
#endif

// Having vector implies having deque
#ifdef HAVE_VECTOR
#include <vector>
#include <deque>
#elif defined(HAVE_VECTOR_H)
#include <vector.h>
#include <deque.h>
#else
// Try to include it anyway
#include <vector>
#include <deque>
#endif

#ifdef HAVE_ALGORITHM
#include <algorithm>
#elif defined(HAVE_ALGO_H)
#include <algo.h>
#elif defined(HAVE_ALGORITHM_H)
#include <algorithm.h>
#else
// Try to include it anyway
#include <algorithm>
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
// Try to include it anyway
#include <sstream>
#endif

#ifdef HAVE_FSTREAM
#include <fstream>
#elif defined(HAVE_FSTREAM_H)
#include <fstream.h>
#else
// Try to include it anyway
#include <fstream>
#endif

using namespace std;

// Some compilers don't like this - particularly DEC!
// Others require it - particularily AIX!
#ifndef NO_USING_STD_OSTREAM
using std::ostream;
#endif

using std::map;
using std::vector;

#define NOX_MAX(a,b) ((a)>(b)) ? (a) : (b)
#define NOX_MIN(a,b) ((a)<(b)) ? (a) : (b)

#endif
