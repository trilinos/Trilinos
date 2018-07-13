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

#ifndef TEUCHOS_DEBUG_DEFAULT_ASSERTS_HPP
#define TEUCHOS_DEBUG_DEFAULT_ASSERTS_HPP


#include "Teuchos_Assert.hpp"


#ifdef TEUCHOS_DEBUG
#define TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT() default: TEUCHOS_TEST_FOR_EXCEPT(true)
#else
/** \brief Macro to insert switch default that throws in a debug build.
 *
 * In a non-debug build, however, this does nothing.  This is also helpful for
 * removing code that would otherwise show as not being covered.
 *
 * \ingroup teuchos_language_support_grp
 */
#define TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT()
#endif

// NOTE: Some explaination for the above TEUCHOS_SWITCH_DEFAULT_DEBUG_ASSERT()
// macro:
//
// In a debug build where default: throws always, we can't have a follow-up
// break statement or some compilers (e.g. NVCC) will (correctly) complain
// that the 'break' statement is unreachable.  Older compilers did not do
// this.  However, in a non-debug build, you want to not even put in a
// 'default:' block.  That way, if all of the enum values are not covered,
// then must compilers will issue a warning about that and we want to see that
// warning in a non-debug build.

#ifdef TEUCHOS_DEBUG
#define TEUCHOS_IF_ELSE_DEBUG_ASSERT() else { TEUCHOS_TEST_FOR_EXCEPT(true); }
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
