// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_COMMON_H
#define NOX_COMMON_H

#include "NOX_Config.h"

#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#else
*** MISSING HEADER FILE ***
#endif

#ifdef HAVE_IOMANIP
#include <iomanip>
#elif defined(HAVE_IOMANIP_H)
#include <iomanip.h>
#else
*** MISSING HEADER FILE ***
#endif

#ifdef HAVE_CMATH
#include <cmath>
#elif defined(HAVE_MATH_H)
#include <math.h>
#else
*** MISSING HEADER FILE ***
#endif

#ifdef HAVE_STRING
#include <string>
#else
*** MISSING HEADER FILE ***
#endif


#include <vector>
#include <map>
using std::map;
using std::vector;


// ------------------------------------------------------------------------
#ifdef ROGERS_OLD_STUFF

#if defined(SGI) || defined(SGI64) || defined(SGI32)
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream.h>
#include <strstream.h>
#include <string> // Need SGI's special string class
#include <math.h>
#include <iomanip.h>
#include <string>
using std::string;

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
#include <strstream>
#include <cstring>
#include <string>
#include <cmath>
#include <iomanip>
using namespace std;

#endif

#include <vector>
#include <map>
using std::map;
using std::vector;

#endif // ROGERS_OLD_STUFF
// ------------------------------------------------------------------------

#endif
