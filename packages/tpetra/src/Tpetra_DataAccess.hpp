/*Paul
27-May-2002 General cleanup. Checked for newNamingConvention (nothing changed).
06-August-2002 Changed to images (nothing changed).
*/

#ifndef _TPETRA_DATAACCESS_HPP_
#define _TPETRA_DATAACCESS_HPP_
/*! \file Tpetra_DataAccess.hpp 
    \brief Tpetra_DataAccess Mode enumerable type
*/

/*! \enum Tpetra_DataAccess
    If set to Copy, user data will be copied at construction.
    If set to View, user data will be encapsulated and used throughout
    the life of the object.
*/
enum Tpetra_DataAccess {Tpetra_Copy, /*!< User data will be copied at
                                  construction. */
                       Tpetra_View /*!< User data will be encapsulated and
                                 used throughout the life of the object. */
                       };

#endif /* _TPETRA_DATAACCESS_HPP_ */
