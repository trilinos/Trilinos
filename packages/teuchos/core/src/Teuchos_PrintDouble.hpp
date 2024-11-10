// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_PRINT_DOUBLE_HPP
#define TEUCHOS_PRINT_DOUBLE_HPP

#include <iosfwd>

/*! \file Teuchos_PrintDouble.hpp
    \brief Declares Teuchos::print_double
*/

namespace Teuchos {

/** \brief Prints a double-precision floating-point number exactly using minimal characters.

This function prints the value (v) to the stream (os) in such a way that,
when read back in, it will be bitwise exactly the same floating-point number
that was passed to Teuchos::print_double (regardless of rounding mode).
It will also choose the representation which minimizes the number characters
required to print the value.
*/

void print_double(std::ostream& os, double v);

}

#endif
