#ifndef __ESI_MatrixRowReadAccess_h
#define __ESI_MatrixRowReadAccess_h

namespace esi {

/** This is the base class for reading a row-oriented sparse matrix.

    The esi::MatrixRowReadAccess class is essentially a transpose of the 
    esi::MatrixColReadAccess class.

    All indices are zero-based global indices unless otherwise noted.

    \verbatim
    Change log:

      11/25/2001 RLC Made global substitution from 'row number' 
      to 'zero-base row index' and cleaned up the documentation. 

      10/29/2001 RLC Added comments about the behavior of 
      zero-length arrays (pointers) which are returned as output 
      args, per Ben Allan's suggestion. Cleaned up the documentation 
      for consistency.
    \endverbatim
*/
template<class Scalar, class Ordinal>
class MatrixRowReadAccess : public virtual MatrixData<Ordinal>
{
  public:

  /** Default destructor. */
  virtual ~MatrixRowReadAccess( void ) {};
  
  /** Get the current allocated length for a row of the matrix.  Note 
      that this is not the number of currently-present valid entries.

      @param row    Input. Global row index.
      @param length Output. Allocated length.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowAllocatedLength( Ordinal row, 
                                           Ordinal & length ) = 0;

  /** Get the number of currently-stored non-zeros for a row of the
      matrix.

      @param row        Input. Global row index.
      @param numEntries Output. Number of coefficient and column-index entries.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getRowNonzeros( Ordinal row, 
                                    Ordinal & numEntries ) = 0;

  /** Get a copy of the coefficients and column indices for a row of
      the matrix. Note that this copies data into caller-allocated memory.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      @param row        Input. Global row index.
      @param coefs      Output. Coefficient data in caller-allocated list.
      @param colIndices Output. Global column indices in the 
                        caller-allocated list.
      @param length     Input. Length of the caller-allocated lists.
      @param rowLength  Output. Length of this row. If rowLength is greater 
                        than length, then only the first 'length' entries in 
                        the row will be copied into the coefs and colIndices 
                        lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copyOutRow( Ordinal row, 
                                Scalar * coefs, 
                                Ordinal * colIndices, 
                                Ordinal length, 
                                Ordinal & rowLength ) = 0;

  /** Get a copy of the coefficients for a row of the matrix. Note that 
      this copies data into caller-allocated memory.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      @param row       Input. Global row index.
      @param coefs     Output. Coefficient data in the caller-allocated list.
      @param length    Input. Length of the caller-allocated lists.
      @param rowLength Output. Length of this row. If rowLength is greater 
                       than length, then only the first 'length' entries in 
                       the row will be copied into the coefs and colIndices 
                       lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copyOutRowCoefficients( Ordinal row, 
                                            Scalar * coefs, 
                                            Ordinal length, 
                                            Ordinal & rowLength ) = 0;

  /** Get a copy of the column indices for a row of the matrix. Note that 
      this copies data into caller-allocated memory.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      @param row        Input. Global row index.
      @param colIndices Output. Global column-indices in the 
                        caller-allocated list.
      @param length     Input. Length of the caller-allocated lists.
      @param rowLength  Output. Length of this row. If rowLength is greater 
                        than length, then only the first 'length' entries in 
                        the row will be copied into the coefs and colIndices 
                        lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copyOutRowIndices( Ordinal row, 
                                       Ordinal * colIndices, 
                                       Ordinal length, 
                                       Ordinal & rowLength ) = 0;

  /** Query whether the 'loadComplete' function has been called on this 
      matrix.

      @param state Output. True if loadComplete has been called.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode isLoaded( bool & state ) = 0;

  /** Query whether the internal arrays have been allocated for this 
      matrix.

      @param state Output. True if arrays have been allocated.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode isAllocated( bool & state ) = 0;

};     // esi::MatrixRowReadAccess class
};     // esi namespace
#endif // __ESI_MatrixRowReadAccess_h

