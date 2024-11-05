// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// ////////////////////////////////////////////////////////////
// Teuchos_StandardMemberCompositionMacros.hpp

#ifndef TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H
#define TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H

/*! \file Teuchos_StandardMemberCompositionMacros.hpp \brief Macro
        that adds <<std member comp>> members as attribute members for any
        class.
*/
#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos { class DummyDummyClass; } // see below

/** \brief Macro that adds <<std member comp>> attributes to any class
 *
 * For example, if you want to include a <<std member comp>> attribute
 * as a member object of type MyClass with the name my_attribute you
 * would include the macro in the public section of YourClass
 * declaration as follows:
 \verbatim

  class YourClass {
  public:
    STANDARD_MEMBER_COMPOSITION_MEMBERS( MyClass, my_attribute );
  };
 \endverbatim
 * This macro adds the following data member to the class declaration:
 \verbatim
        private:
                MyClass my_attribute_;
 \endverbatim
 * and the following methods to your class declaration:
 \verbatim
  public:
  void my_attribute( const My_Class & my_attribute_in )
    { my_attribute_ = my_attribute_in; }
  const My_Class& my_attribute() const
    { return my_attribute_; }
 \endverbatim
 * The advantage of using this type of declaration is that it saves
 * you a lot of typing and space.  Later if you need to override these
 * operations you can just implement the member functions by hand.
 */
#define STANDARD_MEMBER_COMPOSITION_MEMBERS( TYPE, NAME )\
  void NAME ( const TYPE & NAME ## _in ) { NAME ## _ = NAME ## _in ; }\
  const TYPE& NAME() const { return NAME ## _; }\
private:\
  TYPE NAME ## _;\
public: \
  typedef ::Teuchos::DummyDummyClass NAME ## DummyDummyClass_t

// Note: Above, the 'using Teuchos::DummyDummyClass' statement is just there
// to allow (and require) the macro use of the form:
//
//  STANDARD_MEMBER_COMPOSITION_MEMBERS( MyClass, my_attribute );
//
// which allows a semicolon at the end!
//

#endif  // TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H
