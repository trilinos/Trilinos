// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_TEST_FOR_EXCEPTION_H
#define TEUCHOS_TEST_FOR_EXCEPTION_H

/*! \file Teuchos_TestForException.hpp
\brief DEPRECATED (Use Teuchos_Assert.hpp).
*/

#include "Teuchos_Assert.hpp"


#ifdef __GNUC__

#  warning The non-namespaced macros in the included file \
 Teuchos_TestForException.hpp are deprecated. \
 Please use the namespaced macros in the file Teuchos_Assert.hpp.

#endif


/** \brief Deprecated. */
TEUCHOS_DEPRECATED inline
void TestForException_incrThrowNumber()
{
  Teuchos::TestForException_incrThrowNumber();
}


/** \brief Deprecated. */
TEUCHOS_DEPRECATED inline
int TestForException_getThrowNumber()
{
  return Teuchos::TestForException_getThrowNumber();
}


/** \brief Deprecated. */
TEUCHOS_DEPRECATED inline
void TestForException_break( const std::string &msg )
{
  Teuchos::TestForException_break(msg);
}


/** \brief Deprecated. */
TEUCHOS_DEPRECATED inline
void TestForException_setEnableStacktrace(bool enableStrackTrace)
{
  Teuchos::TestForException_setEnableStacktrace(enableStrackTrace);
}


/** \brief Deprecated. */
TEUCHOS_DEPRECATED inline
bool TestForException_getEnableStacktrace()
{
  return Teuchos::TestForException_getEnableStacktrace();
}


/** \brief Deprecated. */
#define TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg) \
  TEUCHOS_TEST_FOR_EXCEPTION(throw_exception_test, Exception, msg)


/** \brief Deprecated. */
#define TEST_FOR_EXCEPTION_CLASS_FUNC(throw_exception_test, Exception, msg) \
  TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(throw_exception_test, Exception, msg)


/** \brief Deprecated. */
#define TEST_FOR_EXCEPTION_PURE_MSG(throw_exception_test, Exception, msg) \
  TEUCHOS_TEST_FOR_EXCEPTION_PURE_MSG(throw_exception_test, Exception, msg)


/** \brief Deprecated. */
#define TEST_FOR_EXCEPT(throw_exception_test) \
  TEUCHOS_TEST_FOR_EXCEPT(throw_exception_test)


/** \brief Deprecated. */
#define TEST_FOR_EXCEPT_MSG(throw_exception_test, msg) \
  TEUCHOS_TEST_FOR_EXCEPT_MSG(throw_exception_test, msg)


/** \brief Deprecated. */
#define TEST_FOR_EXCEPTION_PRINT(throw_exception_test, Exception, msg, out_ptr) \
  TEUCHOS_TEST_FOR_EXCEPTION_PRINT(throw_exception_test, Exception, msg, out_ptr)


/** \brief Deprecated. */
#define TEST_FOR_EXCEPT_PRINT(throw_exception_test, out_ptr) \
  TEUCHOS_TEST_FOR_EXCEPT_PRINT(throw_exception_test, out_ptr)


#endif // TEUCHOS_TEST_FOR_EXCEPTION_H
