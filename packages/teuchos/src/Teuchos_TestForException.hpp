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

#ifndef TEUCHOS_TEST_FOR_EXCEPTION_H
#define TEUCHOS_TEST_FOR_EXCEPTION_H

/*! \file Teuchos_TestForException.hpp
    \brief Macro for throwing an exception with breakpointing to ease debugging.
*/

#include "Teuchos_ConfigDefs.hpp"

/// The only purpose for this function is to set a breakpoint in.
void TestForException_break();

/// \defgroup TestForException_grp  Macro for throwing an exception with breakpointing 
///@{
/** 
 * \brief Macro for throwing an exception with breakpointing to ease debugging
 *
 * @param  throw_exception_test
 *               [in] Test for when to throw the exception.  This can and
 *               should be an expression that may mean something to the user.
 *               The text verbatim of this expression is included in the
 *               formed error string.
 * @param  Exception
 *               [in] This should be the name of an exception class.  The
 *               only requirement for this class is that it have a constructor
 *               that accepts and null terminated C string (i.e. <tt>const char*</tt>).
 * @param  msg   [in] This is any expression that can be included in an
 *               output stream operation.  This is useful when buinding
 *               error messages on the fly.  Note that the code in this
 *               argument only gets evaluated if <tt>throw_exception_test</tt>
 *               evaluates to <tt>true</tt> when an exception is throw.
 *
 * The way that this macro is intended to be used is to 
 * call it in the source code like a function.  For example,
 * suppose that in a piece of code in the file <tt>my_source_file.cpp</tt>
 * that the exception <tt>std::out_of_range</tt> is thrown if <tt>n > 100</tt>.
 * To use the macro, the source code would contain (at line 225
 * for instance):
 \verbatim

 TEST_FOR_EXCEPTION( n > 100, std::out_of_range
    , "Error, n = " << n << is bad" );
 \endverbatim
 * When the program runs and with <tt>n = 125 > 100</tt> for instance,
 * the <tt>std::out_of_range</tt> exception would be thrown with the
 * error message:
 \verbatim

 /home/bob/project/src/my_source_file.cpp:225: n > 100: Error, n = 125 is bad
 \endverbatim
 *
 * In order to debug this, simply open your debugger (gdb for instance),
 * set a break point at <tt>my_soure_file.cpp:225</tt> and then set the condition
 * to break for <tt>n > 100</tt> (e.g. in gdb the command
 * is <tt>cond break_point_number n > 100</tt> and then run the
 * program.  The program should stop a the point in the source file
 * right where the exception will be thrown at but before the exception
 * is thrown.  Try not to use expression for <tt>throw_exception_test</tt> that
 * includes virtual function calls, etc. as most debuggers will not be able to check
 * these types of conditions in order to stop at a breakpoint.  For example,
 * instead of:
 \verbatim

 TEST_FOR_EXCEPTION( obj1->val() > obj2->val(), std::logic_error, "Oh no!" );
 \endverbatim
 * try:
 \verbatim

 double obj1_val = obj1->val(), obj2_val = obj2->val();
 TEST_FOR_EXCEPTION( obj1_val > obj2_val, std::logic_error, "Oh no!" );
 \endverbatim
 * If the developer goes to the line in the source file that is contained
 * in the error message of the exception thrown, he/she will see the
 * underlying condition.
 *
 * As an alternative, you can set a breakpoint for any exception thrown
 * by setting a breakpoint in the function <tt>ThrowException_break()</tt>.
 */
#define TEST_FOR_EXCEPTION(throw_exception_test,Exception,msg) \
{ \
    const bool throw_exception = (throw_exception_test); \
    if(throw_exception) { \
        TestForException_break(); \
	    TeuchosOStringStream omsg; \
	    omsg << __FILE__ << ":" << __LINE__ << ": " \
             << #throw_exception_test << ": " << msg; \
	    throw Exception(TEUCHOS_OSTRINGSTREAM_GET_C_STR(omsg)); \
    } \
}

///
/** \ingroup TestForException_grp
 * \brief This macro is designed to be a short version of
 * <tt>TEST_FOR_EXCEPTION()</tt> that is easier to call.
 *
 * @param  throw_exception_test
 *               [in] Test for when to throw the exception.  This can and
 *               should be an expression that may mean something to the user.
 *               The text verbatim of this expression is included in the
 *               formed error string.
 *
 * \note The exception thrown is <tt>std::logic_error</tt>.
 */
#define TEST_FOR_EXCEPT(throw_exception_test) TEST_FOR_EXCEPTION(throw_exception_test,std::logic_error,"Error!")

///@}

#endif // TEUCHOS_TEST_FOR_EXCEPTION_H
