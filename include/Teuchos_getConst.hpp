// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
