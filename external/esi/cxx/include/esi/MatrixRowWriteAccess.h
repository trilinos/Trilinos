#ifndef __ESI_MatrixRowWriteAccess_h
#define __ESI_MatrixRowWriteAccess_h

namespace esi {

/** This is the base class for filling a row-oriented sparse matrix.

    The esi::MatrixRowWriteAccess class is essentially a transpose 
    of the esi::MatrixColWriteAccess class.

    All indices are zero-based global indices unless otherwise noted.

    \verbatim
    Change log:

      11/25/2001 RLC Made global substitution from 'row number' to 
      'zero-base row index' and cleaned up the documentation. 

      10/29/2001 RLC Added comments about the behavior of zero-length 
      arrays (pointers) which are returned as output args, per Ben 
      Allan's suggestion. Cleaned up the documentation for consistency.

      10/25/2001 RLC Changed method 'copyInRow' to 'copyIntoRow' to be
      consistent with 'sumIntoRow' method naming.  Cleaned up the 
      documentation for style consistency and added 'return 0 if 
      successful' comments where needed.
    \endverbatim
*/
template<class Scalar, class Ordinal>
class MatrixRowWriteAccess : public virtual MatrixData<Ordinal>
{
 public:

  /** Default destructor. */
  virtual ~MatrixRowWriteAccess( void ) {};
  
  /** Get the current allocated length for a row of the matrix.  Note 
      that this is not the number of currently-present valid entries.

      @param row     INPUT: Global row index.
      @param length  OUTPUT: Allocated length.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowAllocatedLength( Ordinal row, 
                                           Ordinal & length ) = 0;

  /** Get the number of currently-stored non-zeros for a row of the
      matrix.

      @param row        INPUT: Global row index.
      @param numEntries OUTPUT: Number of coefficient and column-index entries.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowNonzeros( Ordinal row, 
                                    Ordinal & numEntries ) = 0;

  /** Set all values of the matrix to be the specified scalar. 

      Typically, this method will be used to zero the matrix.

      @param s       INPUT: The value to be set throughout the matrix.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode setAllValues( Scalar s ) = 0;

  /** Copy coefficients for the specified columns into the specified row,
      replacing any values that may already be present for those positions.
      Other values are unchanged.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned for
      the data type.
      
      @param row         INPUT: Global row index.
      @param coefs       INPUT: List of the coefficients.
      @param colIndices  INPUT: List of the global column indices.
      @param length      INPUT: Length of the above lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copyIntoRow( Ordinal row, 
                                 Scalar * coefs, 
                                 Ordinal * colIndices, 
                                 Ordinal length ) = 0;

  /** Sum (accumulate) coefficients for the specified columns into the
      specified row, adding to any values that may already be present for
      those positions.  Other values are unchanged.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned for
      the data type.
      
      @param row        INPUT: Global row index.
      @param coefs      INPUT: List of the coefficients.
      @param colIndices INPUT: List of global column indices.
      @param length     INPUT: Length of the above lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode sumIntoRow( Ordinal row, 
                                Scalar * coefs, 
                                Ordinal * colIndices, 
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
      according to the given list of row lengths. 
      
      Some implementations may not allow this operation, allocation sizes 
      may be required at construction time. Other implementations may allow 
      rows to be allocated and re-allocated individually.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned for
      the data type.
      
      @param rowLengths INPUT: List of length 'localSize', containing the
                        lengths to be allocated for each (local) row of the 
                        matrix.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode allocate( Ordinal * rowLengths ) = 0;

  /** Request that the matrix object's internal arrays be allocated 
      according to the given row length. This requests that all rows have the
      same length. 

      Some implementations may not allow this operation, allocation sizes 
      may be required at construction time. Other implementations may allow 
      rows to be allocated and re-allocated individually.

      @param rowLengths INPUT: List of length 'localSize', containing the
                        lengths to be allocated for each local row of the 
                        matrix.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode allocateRowsSameLength( Ordinal rowLengths ) = 0;

  /** Set the length for a particular (locally-owned) row of the matrix.
      Some implementations may not allow this operation, if data is stored in
      contiguous memory, etc.

      @param row     INPUT: Global row index.
      @param length  INPUT: Length for the specified row.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode setRowLength( Ordinal row, 
                                  Ordinal length ) = 0;

};     // esi::MatrixRowWriteAccess class
};     // esi namespace
#endif // __ESI_MatrixRowWriteAccess_h

