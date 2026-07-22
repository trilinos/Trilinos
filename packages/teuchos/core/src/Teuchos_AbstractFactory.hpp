// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_ABSTRACT_FACTORY_HPP
#define TEUCHOS_ABSTRACT_FACTORY_HPP

#include "Teuchos_RCP.hpp"

namespace Teuchos {

/** \brief Simple, universal "Abstract Factory" interface for the
 * dynamic creation of objects.
 *
 * While <tt>RCP</tt> provides for specialized deallocation
 * policies it does not abstract, in any way, how an object is first
 * allocated.  The most general way to abstract how an object is
 * allocated is to use an "Abstract Factory".  This base class defines
 * the most basic "Abstract Factory" interface and defines only one
 * virtual function, <tt>create()</tt> that returns a
 * <tt>RCP</tt>-wrapped object.
 */
template<class T>
class AbstractFactory {
public:

#ifndef DOXYGEN_COMPILE
	/** \brief . */
	typedef Teuchos::RCP<T>   obj_ptr_t;
#endif

	/** \brief . */
	virtual ~AbstractFactory() {}

	/** \brief Create an object of type T returned as a smart reference
	 * counting pointer object.
	 */
	virtual obj_ptr_t create() const = 0;

}; // class AbstractFactory

} // end Teuchos

#endif // TEUCHOS_ABSTRACT_FACTORY_HPP
