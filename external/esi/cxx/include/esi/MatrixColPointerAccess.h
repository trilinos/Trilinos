#ifndef __ESI_MatrixColPointerAccess_h
#define __ESI_MatrixColPointerAccess_h

namespace esi {

/** Direct-access column-oriented interface for ESI matrix objects.

    The esi::MatrixColPointerAccess interface is not intended to be implemented 
    alone. It should be implemented along with esi::MatrixColWriteAccess and 
    esi::MatrixColReadAccess as a high-performance extension. This interface 
    lacks several necessary query functions which are assumed to be present in 
    the other two classes.

    The esi::MatrixColPointerAccess class is essentially a transpose of the 
    esi::MatrixRowPointerAccess class.

    All indices are zero-based global indices unless otherwise noted.

    Note that there are two sets of functions for obtaining pointers to column 
    data, namely functions for obtaining a 'read-lock', or a 'read-write-lock'.
    There is no practical distinction between these pointers, but the semantics 
    of the functions for obtaining them differ. This interface will grant 
    multiple simultaneous read-locks, but each read-lock must ultimately be 
    matched to a call to 'release*Lock'. On the other hand, only one 
    read-write-lock will be granted at a time. Furthermore, if a read-lock is 
    currently out, a read-write-lock will not be granted. If a lock is requested 
    under conditions where it is not granted, the function will return a 
    non-zero error-code, and the pointer argument will not be referenced.

    \verbatim
    Change log:

      11/25/2001 RLC Cloned this class from the esi::MatrixRowPointerAccess 
      class, and added 'col' into some of the method names to distinguish 
      row vs column orientation of the methods (in the unlikely event that 
      someone implements a class with both column and row oriented data 
      access.  Made global substitution from 'column number' to 'zero-base 
      column index' and cleaned up the documentation. 
    \endverbatim
*/
template<class Scalar, class Ordinal>
class MatrixColPointerAccess : public virtual MatrixData<Ordinal>
{
 public:

  /** Default destructor. */
  virtual ~MatrixColPointerAccess( void ) {};
    
  /** Get pointers to the coefficients and the global row indices for 
      a column.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param column       Input. Global column index.
      @param length       Output. Number of coefficients and row 
                          indices in this column.
      @param coefPointer  Output.
      @param indexPointer Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColPtrReadLock( Ordinal column, 
                                       Ordinal & length, 
                                       Scalar * & coefPointer, 
                                       Ordinal * & indexPointer ) = 0;

  /** Get a pointer to the coefficients for a column.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param column       Input. Global column index.
      @param length       Output. Number of coefficients in this column.
      @param coefPointer  Output. 
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColCoefPtrReadLock( Ordinal column, 
                                           Ordinal & length, 
                                           Scalar * & coefPointer ) = 0;

  /** Get a pointer to the global row indices for a column.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param column       Input. Global column index.
      @param length       Output. Number of coefficients in this column.
      @param indexPointer Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColIndicesPtrReadLock( Ordinal column, 
                                              Ordinal & length, 
                                              Ordinal * & indexPointer ) = 0;

  /** Get pointers to the coefficients and global row indices for 
      a column.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param column       Input. Global column index.
      @param length       Output. Number of coefficients and row-indices
                          in this column.
      @param coefPointer  Output. 
      @param indexPointer Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColPtrReadWriteLock( Ordinal column, 
                                            Ordinal & length, 
                                            Scalar * & coefPointer, 
                                            Ordinal * & indexPointer ) = 0;

  /** Get a pointer to the coefficients for a column.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param column       Input. Global column index.
      @param length       Output. Number of coefficients in this column.
      @param coefPointer  Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColCoefPtrReadWriteLock( Ordinal column, 
                                                Ordinal & length, 
                                                Scalar * & coefPointer ) = 0;

  /** Get pointer to just the global row indices for a column.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param column       Input. Global column index.
      @param length       Output. Number of coefficients in this column.
      @param indexPointer Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColIndicesPtrReadWriteLock( Ordinal column, 
                                                   Ordinal & length, 
                                                   Ordinal * & indexPointer ) = 0;

  /** Release previously-obtained coefficient and row-index pointers 
      for a column.

      @param column       Input. Global column index.
      @param coefPointer  Output. Set to 0/NULL.
      @param indexPointer Output. Set to 0/NULL.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode releaseColPtrLock( Ordinal column, 
                                       Ordinal & length, 
                                       Scalar * & coefPointer, 
                                       Ordinal * & indexPointer ) = 0;

  /** Release previously-obtained coefficient pointer for a column.

      @param column       Input. Global column index.
      @param coefPointer  Output. Set to 0/NULL.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode releaseColCoefPtrLock( Ordinal column, 
                                           Ordinal & length, 
                                           Scalar * & coefPointer ) = 0;

  /** Release previously-obtained row-index pointer for a column.

      @param column        Input. Global column index.
      @param indexPointer  Output. Set to 0/NULL.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode releaseColIndicesPtrLock( Ordinal column, 
                                              Ordinal & length, 
                                              Ordinal * & indexPointer ) = 0;

};     // esi::MatrixColPointerAccess class
};     // esi namespace
#endif // __ESI_MatrixColPointerAccess_h
