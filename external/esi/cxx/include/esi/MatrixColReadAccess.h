#ifndef __ESI_MatrixColReadAccess_h
#define __ESI_MatrixColReadAccess_h

namespace esi {

/** This is the base class for reading a column-oriented sparse matrix.

    The esi::MatrixColReadAccess class is essentially a transpose of the 
    esi::MatrixRowReadAccess class.

    All indices are zero-based global indices unless otherwise noted.

    \verbatim
    Change log:

      11/25/2001 RLC Cloned this class from esi::MatrixRowReadAccess 
      class, and added 'col' into some of the method names to 
      distinguish row vs column orientation of the methods (in the 
      unlikely event that someone implements a class with both column 
      and row oriented data access.  Made global substitution from 
      'column number' to 'zero-base column index' and cleaned up the 
      documentation. 
    \endverbatim
*/
template<class Scalar, class Ordinal>
class MatrixColReadAccess : public virtual MatrixData<Ordinal>
{
  public:

  /** Default destructor. */
  virtual ~MatrixColReadAccess( void ) {};
  
  /** Get the current allocated length for a column of the matrix.
      Note that this is not the number of currently-present valid 
      entries.

      @param column    Input. Global column index.
      @param length    Output. Allocated column length.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColAllocatedLength( Ordinal column, 
                                           Ordinal & length ) = 0;

  /** Get the number of currently-stored non-zeros for a column of the
      matrix.

      @param column     Input. Global column index.
      @param numEntries Output. Number of coefficient and row-index 
                        entries.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode getColNonzeros( Ordinal column, 
                                    Ordinal & numEntries ) = 0;

  /** Get a copy of the coefficients and row indices for a column of
      the matrix.  Note that this copies data into caller-allocated 
      memory.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      @param column      Input. Global column index.
      @param coefs       Output. Coefficient data in caller-allocated list.
      @param rowIndices  Output. Global row indices in caller-allocated list.
      @param length      Input. Length of the caller-allocated lists.
      @param colLength   Output. Length of this column. If colLength is greater 
                         than length, then only the first 'length' entries in 
                         the column will be copied into the coefs and 
                         rowIndices lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copyOutCol( Ordinal column, 
                                Scalar * coefs, 
                                Ordinal * rowIndices, 
                                Ordinal length, 
                                Ordinal & colLength ) = 0;

  /** Get a copy of the coefficients for a column of the matrix.  Note 
      that this copies data into caller-allocated memory.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      @param column    Input. Global column index.
      @param coefs     Output. Coefficient data in caller-allocated list.
      @param length    Input. Length of the caller-allocated lists.
      @param colLength Output. Length of this column. If colLength is 
                       greater than length, then only the first 'length' 
                       entries in the column will be copied into the coefs 
                       and rowIndices lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copyOutColCoefficients( Ordinal column, 
                                            Scalar * coefs, 
                                            Ordinal length, 
                                            Ordinal & colLength ) = 0;

  /** Get a copy of the row indices for a column of the matrix. Note 
      that this copies data into caller-allocated memory.

      Even if the supplied data has 0 length, data pointer arguments to 
      this function shall not be 0/NULL and shall be properly aligned 
      for the data type.
      
      @param column     Input. Global column index.
      @param rowIndices Output. Global row indices in the caller-allocated 
                        list.
      @param length     Input. Length of the caller-allocated lists.
      @param colLength  Output. Length of this column. If colLength is greater 
                        than length, then only the first 'length' entries in 
                        the column will be copied into the coefs and 
                        rowIndices lists.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode copyOutColIndices( Ordinal column, 
                                       Ordinal * rowIndices, 
                                       Ordinal length, 
                                       Ordinal & colLength ) = 0;

  /** Query whether the 'loadComplete' function has been called on this 
      matrix.

      @param state   Output. True if loadComplete has been called.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode isLoaded( bool & state ) = 0;

  /** Query whether the internal arrays have been allocated for this 
      matrix.

      @param state   Output. True if arrays have been allocated.
      @return ErrorCode = 0 if successful.
  */
  virtual ErrorCode isAllocated( bool & state ) = 0;

};     // esi::MatrixColReadAccess class
};     // esi namespace
#endif // __ESI_MatrixColReadAccess_h

