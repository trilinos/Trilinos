// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#ifndef TEUCHOS_LOCAL_TESTING_HELPERS_HPP
#define TEUCHOS_LOCAL_TESTING_HELPERS_HPP


/*! \file Teuchos_LocalTestingHelpers.hpp
 *
 * \brief Utilities to make writing tests easier.
 *
 * <b>WARNING!</b> These macros are not namespaced, so you must
 * only include it in *.cpp files for unit testing only!
 *
 */


#include "Teuchos_TestingHelpers.hpp"


/** \brief . */
#define ECHO( statement ) \
  TEUCHOS_ECHO( statement, out )


/** \brief . */
#define TEST_ASSERT( v1 ) \
  TEUCHOS_TEST_ASSERT( v1, out, success )


/** \brief . */
#define TEST_EQUALITY_CONST( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY_CONST( v1, v2, out, success )


/** \brief . */
#define TEST_EQUALITY( v1, v2 ) \
  TEUCHOS_TEST_EQUALITY( v1, v2, out, success )


/** \brief . */
#define TEST_INEQUALITY_CONST( v1, v2 ) \
  TEUCHOS_TEST_INEQUALITY_CONST( v1, v2, out, success )


/** \brief . */
#define TEST_INEQUALITY( v1, v2 ) \
  TEUCHOS_TEST_INEQUALITY( v1, v2, out, success )


/** \brief . */
#define TEST_FLOATING_EQUALITY( v1, v2, tol ) \
  TEUCHOS_TEST_FLOATING_EQUALITY( v1, v2, tol, out, success )


/** \brief . */
#define TEST_ITER_EQUALITY( iter1, iter2 ) \
  TEUCHOS_TEST_ITER_EQUALITY( iter1, iter2, out, success )


/** \brief . */
#define TEST_ARRAY_ELE_EQUALITY( a, i, val ) \
   TEUCHOS_TEST_ARRAY_ELE_EQUALITY( a, i, val, false, out, local_success )


/** \brief . */
#define TEST_ARRAY_ELE_INEQUALITY( a, i, val ) \
   TEUCHOS_TEST_ARRAY_ELE_INEQUALITY( a, i, val, false, out, local_success )


/** \brief . */
#define TEST_COMPARE( v1, comp, v2 ) \
  TEUCHOS_TEST_COMPARE( v1, comp, v2, out, success )


/** \brief . */
#define TEST_COMPARE_ARRAYS( a1, a2 ) \
  { \
    const bool l_result = compareArrays(a1,#a1,a2,#a2,out); \
    if (!l_result) success = false; \
  }


/** \brief . */
#define TEST_COMPARE_FLOATING_ARRAYS( a1, a2, tol ) \
  { \
    const bool result = compareFloatingArrays(a1,#a1,a2,#a2,tol,out); \
    if (!result) success = false; \
  }


/** \brief . */
#define TEST_THROW( code, ExceptType  ) \
  TEUCHOS_TEST_THROW( code, ExceptType, out, success  )


/** \brief . */
#define TEST_NOTHROW( code  ) \
  TEUCHOS_TEST_NOTHROW( code, out, success  )


#endif  // TEUCHOS_LOCAL_TESTING_HELPERS_HPP
