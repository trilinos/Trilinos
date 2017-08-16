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


//
// Implementations of below macros
//


#ifdef __NVCC__
  #define TEUCHOS_UNREACHABLE_RETURN_IMPL(dummyReturnVal)
#else
  #define TEUCHOS_UNREACHABLE_RETURN_IMPL(dummyReturnVal) \
    return dummyReturnVal
#endif
// Above, the unreachable return statement is only removed for __NVCC__.  This
// is because only __NVCC__ warns about the unreachable return.  Having that
// return statement added for other compilers that don't warn leads to safer
// code by avoiding undefined behavior in case a function is modified in such
// a way that this return is no longer unreachable.  In that case, we want it
// to return something to avoid undefined behavior in case it ever gets
// executed.  That is just good defensive programming.  Also, by leaving this
// return statement by default, we avoid the (incorrect) warning by some
// compilers that the function may not return a value.


//
// User interface macros
//


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
 * compilers will provide a warning that the function may not return a value.
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
 * other compilers will correctly warn that <tt>return -1;</tt> will never be
 * executed with a warning like "statement is unreachable".  Therefore, to
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
 * On compilers that warn about the return being unreachable the return
 * statement is skipped.  On every other compiler, the return statement is
 * kept which results in safer code under refactoring (by avoiding undefined
 * behavior when returning from a function by fall-through without returning
 * an explicit value.
 *
 * \ingroup teuchos_language_support_grp
 */
#define TEUCHOS_UNREACHABLE_RETURN(dummyReturnVal) \
  TEUCHOS_UNREACHABLE_RETURN_IMPL(dummyReturnVal)


#endif // TEUCHOS_COMPILER_CODE_TWEAK_MACROS_HPP
