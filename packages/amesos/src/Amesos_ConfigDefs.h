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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef AMESOS_CONFIGDEFS
#define AMESOS_CONFIGDEFS

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

#include "Amesos_config.h"
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_VECTOR
#include <vector>
#else
  Amesos requires STL vector class
#endif

  //  Disable Kundert for now (we test the KundertOO interface, 
  //  but support the Epetra_CrsKundertSparse interface
#undef HAVE_AMESOS_KUNDERT

#define AMESOS_PRINT(variable) { { \
                      if ( debug_ != 0) { cerr << "AMESOS_PRINT " << # variable << "= " << variable << ", " \
                           << __FILE__ << ", line " << __LINE__ << endl; }  }\
                   }

// prints out an error message if variable is not zero,
// and return this value.
#define AMESOS_CHK_ERR(amesos_err) \
{ if (amesos_err != 0) { \
  std::cerr << "AMESOS ERROR " << amesos_err << ", " \
    << __FILE__ << ", line " << __LINE__ << std::endl; \
    return(amesos_err);  } }

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
  std::cerr << "AMESOS ERROR " << amesos_err << ", " \
    << __FILE__ << ", line " << __LINE__ << std::endl; \
    return(amesos_err);  }


#endif 
