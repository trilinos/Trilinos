// 27-May-2002 General cleanup. Checked for newNamingConvention (nothing changed).

#ifndef _TPETRA_DATAACCESS_H_
#define _TPETRA_DATAACCESS_H_
/*! \file Tpetra_DataAccess.h 
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

#endif /* _TPETRA_DATAACCESS_H_ */
