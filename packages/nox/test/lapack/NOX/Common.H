
//@HEADER

/*!
  \file APPSPACK_Common.hpp
  \brief All standard C/C++ and STL header files are included via this
  file, and all "using" declarations are declared here.

  The purpose of this file is to handle all the special cases for the
  way the various files are included and also to enforse general
  conventions. For example, on some systems it is vital to declare
  <code>using namespace std;</code>, while on other systems it is just
  the opposite. All these situations are handled in this file.
*/

#ifndef APPSPACK_COMMON_H
#define APPSPACK_COMMON_H

#ifdef APSPACK_HAVE_CONFIG_H
#include "APPSPACK_Config.hpp"
#endif

#ifdef HAVE_IOSTREAM
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#else
#include <iostream>
#endif

#ifdef HAVE_IOMANIP
#include <iomanip>
#elif defined(HAVE_IOMANIP_H)
#include <iomanip.h>
#else
#include <iomanip>
#endif

#ifdef HAVE_CMATH
#include <cmath>
#elif defined(HAVE_MATH_H)
#include <math.h>
#else
#include <math.h>
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>	
#elif defined(HAVE_STDIO_H)
#include <stdio.h>
#else		
#include <stdio.h>
#endif

#ifdef HAVE_CSTDLIB
#include <cstdlib>	
#elif defined(HAVE_STDLIB_H)
#include <stdlib.h>
#else		
#include <stdlib.h>
#endif

#ifdef HAVE_CTIME
#include <ctime>	
#elif defined(HAVE_TIME_H)
#include <time.h>
#else		
#include <time.h>
#endif

#ifdef HAVE_CSTRING
#include <cstring>	
#elif defined(HAVE_STRING_H)
#include <string.h>
#else		
#include <string.h>
#endif

#include <unistd.h>		
#include <string>
#include <vector>
#include <list>
#include <fstream>
#include <map>
#include <sstream>

#ifdef REQUIRE_USING_NAMESPACE_STD
using namespace std;
#endif

using std::vector;

#endif
