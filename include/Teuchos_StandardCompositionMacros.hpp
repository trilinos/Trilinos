// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STANDARD_COMPOSITION_MACROS_HPP
#define TEUCHOS_STANDARD_COMPOSITION_MACROS_HPP

/*! \file Teuchos_StandardCompositionMacros.hpp \brief Macro that adds
	<<std comp>> members as attribute members for any class.
*/

#include "Teuchos_RCP.hpp"

/** \brief Macro that adds <<std comp>> members for a composition association.
 *
 * This form is for when the object being held will have const attributes
 * the same as the <tt>this</tt> object.
 *
 * For example, if you want to include a <<std comp>> association
 * with an non-const object of type MyClass of the name my_object you
 * would include the macro in the public section of YourClass
 * declaration as follows:
 *
 \verbatim
 class YourClass {
 public:
   STANDARD_COMPOSITION_MEMBERS( MyClass, my_object );
 };
 \endverbatim
 *
 * Note that the macro adds the following data member
 * to the class declaration:<br>
 \verbatim
 private:
   Teuchos::RCP< TYPE > NAME_;		
 \endverbatim
 *
 * \ingroup StandardContainmentMacros_grp
 */
#define STANDARD_COMPOSITION_MEMBERS( TYPE, NAME ) \
	void set_ ## NAME (const Teuchos::RCP< TYPE >& NAME ## _in ) \
	{	NAME ## _ = NAME ## _in ; } \
	Teuchos::RCP< TYPE > get_ ## NAME() const \
	{	return NAME ## _; } \
	TYPE& NAME() \
	{	return *NAME ## _; } \
	const TYPE& NAME() const \
	{	return *NAME ## _; } \
private: \
	Teuchos::RCP< TYPE > NAME ## _; \
public: \
  typedef Teuchos::RCP< TYPE > NAME ## _ptr_t

/** \brief Macro that adds <<std comp>> members for a composition association.
 *
 * This form is for when the object being held will have non-const attributes
 * irrespective of the const of <tt>this</tt>.
 *
 * For example, if you want to include a <<std comp>> association
 * with an non-const object of type MyClass of the name my_object you
 * would include the macro in the public section of YourClass
 * declaration as follows:
 *
 \verbatim
 class YourClass {
 public:
   STANDARD_NONCONST_COMPOSITION_MEMBERS( MyClass, my_object );
 };
 \endverbatim
 *
 * Note that the macro adds the following data member
 * to the class declaration:<br>
 \verbatim
 private:
   Teuchos::RCP< TYPE > NAME_;		
 \endverbatim
 *
 * \ingroup StandardContainmentMacros_grp
 */
#define STANDARD_NONCONST_COMPOSITION_MEMBERS( TYPE, NAME ) \
	void set_ ## NAME ( const Teuchos::RCP< TYPE >& NAME ## _in ) \
	{	NAME ## _ = NAME ## _in ; } \
	Teuchos::RCP< TYPE > get_ ## NAME() const \
	{	return NAME ## _; } \
	TYPE& NAME() const \
	{	return *NAME ## _; } \
private: \
	Teuchos::RCP< TYPE > NAME ## _; \
public: \
  typedef Teuchos::RCP< TYPE > NAME ## _ptr_t

/** \brief Macro that adds <<std comp>> members for a composition association
 * where the contained object is always constant.
 *
 * This form is for when the object being held will have const attributes
 * irrespective of the const of <tt>this</tt>.
 *
 * For example, if you want to include a <<std comp>> association
 * with a const object of type MyClass of the name my_object you
 * would include the macro in the public section of YourClass
 * declaration as follows:
 *
 \verbatim
 class YourClass {
 public:
   STANDARD_CONST_COMPOSITION_MEMBERS( MyClass, my_object );
 };
 \endverbatim
 *
 * Note that the macro adds the following data member
 * to the class declaration:<br>
 \verbatim
 private:
   NAME_ptr_t NAME_;		
 \endverbatim
 *
 * \ingroup StandardContainmentMacros_grp
 */
#define STANDARD_CONST_COMPOSITION_MEMBERS( TYPE, NAME ) \
public: \
	void set_ ## NAME ( const Teuchos::RCP< const TYPE >& NAME ## _in ) \
	{	NAME ## _ = NAME ## _in ; } \
	Teuchos::RCP< const TYPE > get_ ## NAME() const \
	{	return NAME ## _; } \
	const TYPE& NAME() const \
	{	return *NAME ## _; } \
private: \
	Teuchos::RCP< const TYPE > NAME ## _; \
public: \
  typedef Teuchos::RCP< const TYPE > NAME ## _ptr_t

#endif	// TEUCHOS_STANDARD_COMPOSITION_MACROS_HPP
