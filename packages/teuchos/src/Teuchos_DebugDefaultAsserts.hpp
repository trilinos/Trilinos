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

#ifndef TEUCHOS_DEBUG_DEFAULT_ASSERTS_HPP
#define TEUCHOS_DEBUG_DEFAULT_ASSERTS_HPP


#include "Teuchos_TestForException.hpp"


#ifdef TEUCHOS_DEBUG
#define TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT() default: TEST_FOR_EXCEPT(true); break
#else
/** \brief Macro to insert switch default that throws in a debug build.
 *
 * In a non-debug build, however, this does nothing.  This is also helpful for
 * removing code that would otherwise show as not being covered.
 *
 * \ingroup teuchos_language_support_grp
 */
#define TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT() break
#endif


#ifdef TEUCHOS_DEBUG
#define TEUCHOS_IF_ELSE_DEBUG_ASSERT() else { TEST_FOR_EXCEPT(true); }
#else
/** \brief Macro to insert else block that throws in a debug build.
 *
 * In a non-debug build, however, this does nothing.  This is also helpful for
 * removing code that would otherwise show as not being covered.
 *
 * \ingroup teuchos_language_support_grp
 */
#define TEUCHOS_IF_ELSE_DEBUG_ASSERT() else {}
#endif


#endif // TEUCHOS_DEBUG_DEFAULT_ASSERTS_HPP
