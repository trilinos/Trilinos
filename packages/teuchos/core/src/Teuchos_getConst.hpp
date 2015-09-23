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

#ifndef TEUCHOS_GET_CONST_HPP
#define TEUCHOS_GET_CONST_HPP

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

/** \brief Return a constant reference to an object given a non-const reference.
 * \tparam T The type (without const qualifier) to make const.
 *
 * This template function provides a shorthand for
 * \verbatim const_cast<const T&>(t) \endverbatim
 * as
 * \verbatim getConst(t) \endverbatim
 * This saves you the trouble of typing the name of the class,
 * and also saves readers of your code from looking at the
 * (possibly long) class name.
 *
 * Here is an example showing why you might want to use this function.
 * \code
 * // A class with a very long name.
 * class MyVeryLongClassNameWhichIsIndeedQuiteLongIsItNot { ... };
 *
 * // A function which takes a const reference to an instance of that class.
 * void f1 (const MyVeryLongClassNameWhichIsIndeedQuiteLongIsItNot& x)
 *
 * // A function which takes a nonconst reference to an instance of that class,
 * // but needs to call f1, which takes a const reference.
 * void f2 (MyVeryLongClassNameWhichIsIndeedQuiteLongIsItNot& x) {
 *   // ... do some stuff with x ...
 *   // Call f1.  Notice how few characters this takes,
 *   // compared with a const_cast.
 *   f1 (getConst (x));
 *   // ... do other stuff ...
 * }
 * \endcode
 *
 * \ingroup teuchos_language_support_grp
 */
template<class T>
inline const T& getConst( T& t ) {      return t; }

}       // end namespace Teuchos

#endif // TEUCHOS_GET_CONST_HPP
