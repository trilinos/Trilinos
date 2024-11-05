// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_ENULL_HPP
#define TEUCHOS_ENULL_HPP

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos {

/** \brief Used to initialize a <tt>RCP</tt> object to NULL using an
 * implicit conversion!
 *
 * \relates RCP
 */
enum ENull { null };

} // end namespace Teuchos

#endif	// TEUCHOS_ENULL_HPP
