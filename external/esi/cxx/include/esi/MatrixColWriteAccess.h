#ifndef __ESI_MatrixColwWriteAccess_h
#define __ESI_MatrixColwWriteAccess_h

namespace esi {

/** This is the base class for filling a column-oriented sparse matrix.

    The esi::MatrixColWriteAccess class is essentially a transpose of the 
    esi::MatrixRowWriteAccess class.

    All indices are zero-based global indices unless otherwise noted.

    \verbatim
    Change log:

      11/25/2001 RLC Cloned this class from esi::MatrixRowWriteAccess 
      class, and added 'col' into some of the method names to 
      distinguish row vs column orientation of the methods (in the 
      unlikely event that someone implements a class with both column 
      and row oriented data access. Made global substitution from 
      'column number' to 'zero-base column index' and cleaned up the 
      documentation. 

      10/29/2001 RLC Added comments about the behavior of zero-length 
      arrays (pointers) which are returned as output args, per Ben 
      Allan's suggestion. Cleaned up the documentation for consistency.

      10/25/2001 RLC Changed method 'copyInRow' to 'copyIntoRow' to 
      be consistent with 'sumIntoRow' method naming.  Cleaned up the 
      documentation for style consistency and added 'return 0 if 
      successful' comments where needed.
    \endverbatim
*/
template<class Scalar, class Ordinal>
class MatrixColWriteAccess : public virtual MatrixData<Ordinal>
{
 public:

  /** Default destructor. */
  virtual ~MatrixColWriteAccess( void ) {};
    
  /** Get the current allocated length for a column of the matrix.
      Note that this is not the number of currently-present valid 
      entries.

      @param column     INPUT: Global column index.
      @param length     OUTPUT: Allocated length.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColAllocatedLength( Ordinal column, 
                                           Ordinal & length ) = 0;

  /** Get the number of currently-stored non-zeros for a column 
      of the matrix.

      @param column     INPUT: Global column index.
      @param numEntries OUTPUT: Number of coefficient and row-index 
                        entries.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColNonzeros( Ordinal column, 
                                    Ordinal & numEntries ) = 0;

  /** Set all values of the matrix to be the specified scalar. 

      Typically, this method will be used to zero the matrix.

      @param s   INPUT: The value to be set throughout the matrix.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode setAllValues( Scalar s ) = 0;

  /** Copy coefficients for the specified rows into the specified column,
      replacing any values that may already be present for those positions.
      Other values are unchanged.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned for
      the data type.
      
      @param column      INPUT: Global column index.
      @param coefs       INPUT: List of coefficients.
      @param rowIndices  INPUT: List of global row indices.
      @param length      INPUT: Length of the above lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copyIntoCol( Ordinal column, 
                                 Scalar * coefs, 
                                 Ordinal * rowIndices, 
                                 Ordinal length ) = 0;

  /** Sum (accumulate) coefficients for the specified rows into the
      specified column, adding to any values that may already be present for
      those positions.  Other values are unchanged.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned for
      the data type.
      
      @param column     INPUT: Global column index.
      @param coefs      INPUT: List of coefficients.
      @param rowIndices INPUT: List of global row indices.
      @param length     INPUT: Length of the above lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode sumIntoCol( Ordinal column, 
                                Scalar * coefs, 
                                Ordinal * rowIndices, 
                                Ordinal length ) = 0;

  /** Query whether the 'loadComplete' function has been called on this
      matrix.

      @param state  OUTPUT: True if loadComplete has been called (since 
                    structure was last modified - i.e., if true, then
                    loadComplete need not be called again before using 
                    the matrix).
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode isLoaded( bool & state ) = 0;

  /** Query whether the internal arrays have been allocated for this
      matrix.

      @param state   OUTPUT: True if arrays have been allocated.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode isAllocated( bool & state ) = 0;

  /** Signal the matrix that all data-loading is complete. The object may
      now do any data consolidation and/or global synchronization necessary.

      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode loadComplete( void ) = 0;

  /** Request that the matrix object's internal arrays be allocated
      according to the given list of column lengths. 

      Some implementations may not allow this operation, allocation 
      sizes may be required at construction time. Other implementations 
      may allow columns to be allocated and re-allocated individually.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned for
      the data type.
      
      @param colLengths INPUT: List of length 'localSize', containing the
                        lengths to be allocated for each local column of the 
                        matrix.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode allocate( Ordinal * colLengths ) = 0;

  /** Request that the matrix object's internal arrays be allocated 
      according to the given column length. This requests that all 
      columns have the same length. 

      Some implementations may not allow this operation, allocation 
      sizes may be required at construction time. Other implementations 
      may allow columns to be allocated and re-allocated individually.

      @param colLengths INPUT: List of length 'localSize', containing the
                        lengths to be allocated for each local column of the 
                        matrix.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode allocateColsSameLength( Ordinal colLengths ) = 0;

  /** Set the length for a particular (locally-owned) column of the matrix.

      Some implementations may not allow this operation, if data is 
      stored in contiguous memory, etc.

      @param column     INPUT: Global column index.
      @param length     INPUT: Length for the specified column.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode setColLength( Ordinal column, 
                                  Ordinal length ) = 0;

};     // esi::MatrixColWriteAccess class
};     // esi namespace
#endif // __ESI_MatrixColWriteAccess_h
