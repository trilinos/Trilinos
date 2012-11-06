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

// ////////////////////////////////////////////////////////////
// Teuchos_StandardMemberCompositionMacros.hpp

#ifndef TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H
#define TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H

/*! \file Teuchos_StandardMemberCompositionMacros.hpp \brief Macro
	that adds <<std member comp>> members as attribute members for any
	class.
*/
#include "Teuchos_ConfigDefs.hpp"

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

#endif	// TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H
