/*Paul
Modified 22-Jan-2003
*/

#ifndef _TPETRA_DATAACCESS_HPP_
#define _TPETRA_DATAACCESS_HPP_

/*! \file Tpetra_DataAccess.hpp 
    \brief Tpetra::DataAccess Mode enumerable type
*/

namespace Tpetra {
	
	/*! \enum DataAccess
      If set to Copy, user data will be copied at construction.
      If set to View, user data will be encapsulated and used throughout
      the life of the object.
	*/
	
	enum DataAccess {
		Copy, /*!< User data will be copied at construction. */
		View /*!< User data will be encapsulated and used throughout the life of the object. */
	};

} // namespace Tpetra

#endif /* _TPETRA_DATAACCESS_HPP_ */
