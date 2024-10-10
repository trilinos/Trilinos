// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef AMESOS_CONFIGDEFS
#define AMESOS_CONFIGDEFS

#if defined(Amesos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Amesos package is deprecated"
#endif
#endif

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

#ifdef F77_FUNC
#undef F77_FUNC
#endif

#include "Amesos_config.h"
#include "Epetra_ConfigDefs.h"
#include <vector>

#define AMESOS_PRINT(variable) { { \
                      if ( debug_ != 0) { std::cerr << "AMESOS_PRINT " << # variable << "= " << variable << ", " \
                           << __FILE__ << ", line " << __LINE__ << std::endl; }  }\
                   }

// prints out an error message if variable is not zero,
// and return this value. This is a copy of macro EPETRA_CHK_ERR,
// here modified so that the user sees an "AMESOS ERROR" instead
// of a possibly misleading "Epetra ERROR".

#define AMESOS_CHK_ERR(a) { { int amesos_err = a; \
                              if ((amesos_err < 0 && Epetra_Object::GetTracebackMode() > 0) || \
                                  (amesos_err > 0 && Epetra_Object::GetTracebackMode() > 1)) { \
                      std::cerr << "AMESOS ERROR " << amesos_err << ", " \
                           << __FILE__ << ", line " << __LINE__ << std::endl; }\
                      if (amesos_err != 0) return(amesos_err);  }\
                   }
//#define AMESOS_CHK_ERR(amesos_err) { EPETRA_CHK_ERR( amesos_err ) }

// prints out an error message if variable is not zero,
// returns void
#define AMESOS_CHK_ERRV(amesos_err) \
{ if (amesos_err != 0) { \
  std::cerr << "AMESOS ERROR " << amesos_err << ", " \
    << __FILE__ << ", line " << __LINE__ << std::endl; \
    return;  } }

// prints out an error message if variable is not zero,
// returns void
#define AMESOS_RETURN(amesos_err) \
{ \
  if (amesos_err != 0) \
    std::cerr << "AMESOS ERROR " << amesos_err << ", " \
      << __FILE__ << ", line " << __LINE__ << std::endl; \
  return(amesos_err);  }


#endif
