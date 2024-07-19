// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_FLOATING_POINT_TRAP_HPP
#define TEUCHOS_FLOATING_POINT_TRAP_HPP

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {


/** \defgroup Teuchos_FloatingPointTrap_grp Floating Point Trapping Support Code
 *
 * This is code that makes it easier to trap floating point errors in user
 * code.
 *
 * ToDo: Finish documentation!
 */
//@{


/** \brief Turn on or off a floating point trap.
 *
 * To use this from gdb, set:

 \verbatim
   handle SIGFPE stop nopass
 \endverbatim

 * before running the code in gdb.
 */
void doFloatingPointTrap(bool enableTrap);


//@}


} // namespace Teuchos


#endif // TEUCHOS_FLOATING_POINT_TRAP_HPP
