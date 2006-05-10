
//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

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
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_configdefs_hpp_
#define _Isorropia_configdefs_hpp_

/*
   The macros PACKAGE, PACKAGE_NAME, etc, get defined in the automatically-
   generated header Isorropia_autoheader.h. So we need to undefine them before
   including that header, in order to avoid warnings in cases where another
   package's header is also included and has already defined them.
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

#include <Isorropia_autoheader.h>

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#error "Isorropia must have <iostream>"
#endif

#ifdef HAVE_FSTREAM
#include <fstream>
#else
#error "Isorropia must have <fstream>"
#endif

#ifdef HAVE_EXCEPTION
#include <exception>
#else
#error "Isorropia must have <exception>"
#endif

#ifdef HAVE_VECTOR
#include <vector>
#else
#error "Isorropia must have <vector>"
#endif

#ifdef HAVE_SET
#include <set>
#else
#error "Isorropia must have <set>"
#endif

#ifdef HAVE_TIME_H
#include <time.h>
#endif

#endif //_Isorropia_configdefs_hpp_

