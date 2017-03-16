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

#ifndef TEUCHOS_COMPILER_CODE_TWEAK_MACROS_HPP
#define TEUCHOS_COMPILER_CODE_TWEAK_MACROS_HPP

#include "TeuchosCore_config.h"


/** \file Teuchos_CompilerCodeTweakMacros.hpp
 *
 * This file contains macros used to tweak code to make compilers happy and
 * avoid warnings of various types.  See the macros defined here for the
 * motivation and usage of these macros.
 */


/** \brief Avoid warning about unreachable or missing return from function.
 *
 * Consider a function like:
 *
 * \begin{code}
  int func(const ESomeEnum val)
  {
    switch (val) {
      case VAL1: return 1;
      case VAL2: return 2;
      default: TEUCHOS_TEST_FOR_EXCEPT(true);
    }
  }
 * \end{code} 
 *
 * That code will never execute out of the switch statement.  However, some
 * comilers will provide a warning that the function may not return a value.
 * Therefore, one can remove this warning by adding a dummy return value like:
 *
 * \begin{code}
  int func(const ESomeEnum val)
  {
    switch (val) {
      case VAL1: return 1;
      case VAL2: return 2;
      default: TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    return -1; // Will never get called!
  }
 * \end{code} 
 *
 * That removes the "may not return value" warning on those compilers.  But
 * other compilers will correclty warn that <tt>return -1;</tt> will never be
 * executed with a warning like "statement is reachable".  Therefore, to
 * address warnings like this, this macro is used like:
 *
 * \begin{code}
  int func(const ESomeEnum val)
  {
    switch (val) {
      case VAL1: return 1;
      case VAL2: return 2;
      default: TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    TEUCHOS_UNREACHABLE_RETURN(-1);
  }
 * \end{code} 
 *
 * On compilers that warn about a missing return, the statement <tt>return
 * -1;</tt> is provided.  However, on compilers that warn about unreachable
 * code, it will simply provide <tt>(void);</tt>.
 */
#define TEUCHOS_UNREACHABLE_RETURN(dummyReturnVal)


#endif // TEUCHOS_COMPILER_CODE_TWEAK_MACROS_HPP
