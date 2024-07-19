// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#ifndef _TEUCHOS_DATAACCESS_HPP_
#define _TEUCHOS_DATAACCESS_HPP_

/*! \file Teuchos_DataAccess.hpp
    \brief Teuchos::DataAccess Mode enumerable type
*/

namespace Teuchos {
	
	/*! \enum DataAccess
      If set to Copy, user data will be copied at construction.
      If set to View, user data will be encapsulated and used throughout
      the life of the object.
	*/
	
	enum DataAccess {
		Copy, /*!< User data will be copied at construction. */
		View /*!< User data will be encapsulated and used throughout the life of the object. */
	};

} // namespace Teuchos

#endif /* _TEUCHOS_DATAACCESS_HPP_ */
