#ifndef __ESI_MatrixRowPointerAccess_h
#define __ESI_MatrixRowPointerAccess_h

namespace esi {

/** Direct-access row-oriented interface for ESI matrix objects.

    The esi::MatrixRowPointerAccess interface is not intended to be implemented 
    alone. It should be implemented along with esi::MatrixRowWriteAccess and 
    esi::MatrixRowReadAccess as a high-performance extension. This interface 
    lacks several necessary query functions which are assumed to be present in 
    the other two classes.

    The esi::MatrixRowPointerAccess class is essentially a transpose of the 
    esi::MatrixColPointerAccess class.

    All indices are zero-based global indices unless otherwise noted.

    Note that there are two sets of functions for obtaining pointers to row data,
    namely functions for obtaining a 'read-lock', or a 'read-write-lock'.
    There is no practical distinction between these pointers, but the semantics of
    the functions for obtaining them differ. This interface will grant multiple
    simultaneous read-locks, but each read-lock must ultimately be matched to a
    call to 'release*Lock'. On the other hand, only one read-write-lock will be
    granted at a time. Furthermore, if a read-lock is currently out, a 
    read-write-lock will not be granted. If a lock is requested under conditions
    where it is not granted, the function will return a non-zero error-code, and
    the pointer argument will not be referenced.

    \verbatim
    Change log:

      11/25/2001 RLC Added 'row' into some of the method names to 
      distinguish row vs column orientation of the methods (in the 
      unlikely event that someone implements a class with both column 
      and row oriented data access. Made global substitution from 
      'column number' to 'zero-base column index' and cleaned up the 
      documentation. 

      10/29/2001 RLC Added comments about the behavior of zero-length 
      arrays (pointers) which are returned as output args, per Ben Allan's
      suggestion. Cleaned up the documentation for consistency.
    \endverbatim
*/
template<class Scalar, class Ordinal>
class MatrixRowPointerAccess : public virtual MatrixData<Ordinal>
{
 public:

  /** Default destructor. */
  virtual ~MatrixRowPointerAccess( void ) {};
    
  /** Get pointers to the coefficients and column indices for a row.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param row          Input. Global row index.
      @param length       Output. Number of coefficients and column 
                          indices in this row.
      @param coefPointer  Output.
      @param indexPointer Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowPtrReadLock( Ordinal row, 
                                       Ordinal & length, 
                                       Scalar * & coefPointer, 
                                       Ordinal * & indexPointer ) = 0;

  /** Get a pointer to the coefficients for a row.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param row          Input. Global row index.
      @param length       Output. Number of coefficients in this row.
      @param coefPointer  Output. 
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowCoefPtrReadLock( Ordinal row, 
                                           Ordinal & length, 
                                           Scalar * & coefPointer ) = 0;

  /** Get a pointer to the column-indices for a row.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param row          Input. Global row index.
      @param length       Output. Number of coefficients in this row.
      @param indexPointer Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowIndicesPtrReadLock( Ordinal row, 
                                              Ordinal & length, 
                                              Ordinal * & indexPointer ) = 0;

  /** Get pointers to the coefficients and column indices for a row.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param row          Input. Global row index.
      @param length       Output. Number of coefficients and column
                          indices in this row.
      @param coefPointer  Output. 
      @param indexPointer Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowPtrReadWriteLock( Ordinal row, 
                                            Ordinal & length, 
                                            Scalar * & coefPointer, 
                                            Ordinal * & indexPointer ) = 0;

  /** Get a pointer to the coefficients for a row.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param row          Input. Global row index.
      @param length       Output. Number of coefficients in this row.
      @param coefPointer  Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowCoefPtrReadWriteLock( Ordinal row, 
                                                Ordinal & length, 
                                                Scalar * & coefPointer ) = 0;

  /** Get a pointer to the column indices for a row.

      As with ANSI malloc, each data pointer retrieved with this function
      shall not be 0/NULL and shall be properly aligned for the data type
      even if the data has 0 length.

      @param row          Input. Global row index.
      @param length       Output. Number of coefficients in this row.
      @param indexPointer Output.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowIndicesPtrReadWriteLock( Ordinal row, 
                                                   Ordinal & length, 
                                                   Ordinal * & indexPointer ) = 0;

  /** Release previously-obtained coefficient and column-index pointers 
      for a row.

      @param row          Input. Global row index.
      @param coefPointer  Output. Set to 0/NULL.
      @param indexPointer Output. Set to 0/NULL.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode releaseRowPtrLock( Ordinal row, 
                                       Ordinal & length, 
                                       Scalar * & coefPointer, 
                                       Ordinal * & indexPointer ) = 0;

  /** Release a previously-obtained coefficient pointer for a row.

      @param row          Input. Global row index.
      @param coefPointer  Output. Set to 0/NULL.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode releaseRowCoefPtrLock( Ordinal row, 
                                           Ordinal & length, 
                                           Scalar * & coefPointer ) = 0;
  
  /** Release a previously-obtained column-index pointer for a row.

      @param row          Input. Global row index.
      @param indexPointer Output. Set to 0/NULL.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode releaseRowIndicesPtrLock( Ordinal row, 
                                              Ordinal & length, 
                                              Ordinal * & indexPointer ) = 0;

};     // esi::MatrixRowPointerAccess class
};     // esi namespace
#endif // __ESI_MatrixRowPointerAccess_h
