// ////////////////////////////////////////////////////////////
// Teuchos_StandardMemberCompositionMacros.hpp

#ifndef TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H
#define TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H

#include "Teuchos_ConfigDefs.hpp"

///
/** Macro that addes << std member comp >> members for attribute member.
  * \ingroup StandardContainmentMacros_grp
  *
  * For example, if you want to include a <<std member comp>> attribute
  * as a member object of type MyClass with the name my_attribute you
  * would include the macro in the public section of YourClass
  * declaration as follows:
  \verbatim
	class YourClass {
	public:
		STANDARD_MEMBER_COMPOSITION_MEMBERS( MyClass, my_attribute )
	};
  \endverbatim
  * Note that the macro adds the following data member
  * to the class declaration:
  \verbatim
	private:
		TYPE NAME_;
  \endverbatim
  * The advantage of using this type of declaration is that it saves
  * you a lot of typing and space.  Later if you need to override these
  * operations you can just implement the member functions by hand.
  */
#define STANDARD_MEMBER_COMPOSITION_MEMBERS( TYPE, NAME )			\
public:																\
	void NAME (const TYPE & NAME )									\
	{	NAME ## _ = NAME ; }										\
	const TYPE& NAME() const										\
	{	return NAME ## _; }											\
private:															\
	TYPE NAME ## _;													\
public:

#endif	// TEUCHOS_STANDARD_MEMBER_COMPOSITION_MACROS_H
