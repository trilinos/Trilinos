// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STATIC_SETUP_MACRO_HPP
#define TEUCHOS_STATIC_SETUP_MACRO_HPP


#include "Teuchos_ConfigDefs.hpp"


/** \brief Run setup code statically in a translation unit.
 *
 * NOTE: Make sure the call this in an anonymous namespace as:
 *
 \verbatim

 namespace {

 TEUCHOS_STATIC_SETUP()
 {
   // Some code you want to call before main() runs ...
   ...
 }

 } // namespace

 \endverbatim
 *
 */
#define TEUCHOS_STATIC_SETUP() \
  class StaticSetup { \
  public: \
    StaticSetup(); \
  } staticSetup; \
  \
  StaticSetup::StaticSetup()


#endif  // TEUCHOS_STATIC_SETUP_MACRO_HPP
