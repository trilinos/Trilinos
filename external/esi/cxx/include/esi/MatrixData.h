#ifndef __ESI_MatrixData_h
#define __ESI_MatrixData_h

namespace esi {

/** The ESI MatrixData class.
    
    The esi::MatrixData base class has no data distribution assumptions.
    This is essentially a placeholder from which a particular
    matrix implementation can be derived.
*/
template<class Ordinal>
class MatrixData : public virtual Object
{
 public:

  /** Default destructor. */
  virtual ~MatrixData( void ) {};
    
  /** Query the global ('row' and 'column') sizes.
      @param rows        OUTPUT: The number of rows in the global matrix.
      @param columns     OUTPUT: The number of columns in the global matrix.
      @return error-code 0 if successful.
  */
  virtual ErrorCode getGlobalSizes( Ordinal & rows, 
                                    Ordinal & columns ) = 0;
  
  /** Query the local ('row' and 'column')  sizes.
      @param rows        OUTPUT: The number of local rows.
      @param columns     OUTPUT: The number of local columns.
      @return error-code 0 if successful.
  */
  virtual ErrorCode getLocalSizes( Ordinal & rows, 
                                   Ordinal & columns ) = 0;
  
  /** Get the esi::IndexSpace objects that describe the 
      row-space and the column-space of this matrix.

      @param rowIndexSpace  OUTPUT: The row 'index space' object.
      @param colIndexSpace  OUTPUT: The column 'index space' object.
      @return error-code 0 if successful.
  */
  virtual ErrorCode getIndexSpaces( IndexSpace<Ordinal> * & rowIndexSpace, 
                                    IndexSpace<Ordinal> * & colIndexSpace ) = 0;

};     // esi::MatrixData class
};     // esi namespace
#endif // __ESI_MatrixData_h
