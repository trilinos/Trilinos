#ifndef _Epetra_ESI_CrsMatrix_h_
#define _Epetra_ESI_CrsMatrix_h_

#include "Epetra_Object.h"
#include "Epetra_CrsMatrix.h"

//
//Epetra_ESI_Vector.h includes everything "above" it,
//including Epetra_ESI_Object.h, Epetra_Comm.h, Epetra_SerialComm.h,
//Epetra_MpiComm.h, Epetra_Map.h,
//Epetra_ESI_Vector.h, Epetra_Vector.h, and the ESI headers.
//
#include "Epetra_ESI_Vector.h"


namespace epetra_esi {

/** Petra's ESI RowMatrix implementation.
This class implements these ESI interfaces:
<ul>
<li> esi::Object
<li> esi::Operator
<li> esi::OperatorTranspose
<li> esi::MatrixData
<li> esi::RowMatrixReadAccess
<li> esi::RowMatrixWriteAccess Note that not all of the functions in this
interface are implemented. See individual function descriptions for details.
</ul>
Note that although this class is templated, it may only be instantiated on
the type-pair double,int.
*/

template<class Scalar, class Ordinal>
class CrsMatrix : public virtual esi::Object,
                  public virtual epetra_esi::Object,
                  public virtual esi::Operator<Scalar, Ordinal>,
                  public virtual esi::OperatorTranspose<Scalar, Ordinal>,
                  public virtual esi::MatrixData<Ordinal>,
                  public virtual esi::MatrixRowWriteAccess<Scalar, Ordinal>,
                  public virtual esi::MatrixRowReadAccess<Scalar, Ordinal>
{
 public:
  /** Constructor that takes a fully initialized Epetra_CrsGraph object.
     @param CV Specifies whether this object should take a copy or just a
       view of the incoming graph object. Valid values are simply 'Copy' or
       'View'.
     @param graph This should be a fully initialized object holding the
       complete structure of the matrix.
  */
  CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph& graph);

  /** Constructor that takes a IndexSpace object, and then allows dynamic
    insertion of indices and values into the rows of this matrix.
     @param CV Specifies whether this object should take a copy or just a
       view of the incoming graph object. Valid values are simply 'Copy' or
       'View'.
     @param indexspace Describes the distribution of the rows of this matrix
        across processors.
     @param estimatedNumEntriesPerRow An estimate of the size that each row
        will be. If 0 is provided here, the matrix rows will be dynamically
        allocated with each call to copyIntoRow.
  */
  CrsMatrix(Epetra_DataAccess CV,
            const epetra_esi::IndexSpace<Ordinal>& indexspace,
            Ordinal estimatedNumEntriesPerRow);

  /** Destructor. */
  virtual ~CrsMatrix();


  //
  //esi::Operator functions.
  //

  /** Function for signalling that the initialization and data-loading is
    complete. The epetra_esi::CrsMatrix object now has the internal
    Epetra_CrsMatrix perform the TransformToLocal() operation.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode setup( void )
    { int err = epetra_crsmatrix_->TransformToLocal();
      if (err) EPETRA_ESI_ERR_BEHAVIOR(-1);
      setupHasBeenCalled_ = true;    return(0); }


  /** Function for applying this operator to a vector. (This is a matvec.) If
   the setup() function has not yet been called, it will be called internally,
   automatically at this point.
   @param x Input. The run-time type of this vector is required to be 
      epetra_esi::Vector.
   @param y On exit, the contents of this will be the result of the matvec.
       The run-time of this vector must also be epetra_esi::Vector.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode apply(esi::Vector<Scalar, Ordinal>& x,
                               esi::Vector<Scalar, Ordinal>& y)
    { int err = 0;
      if (!setupHasBeenCalled_) err = setup();
      if (err) return(-1);
      const Epetra_Vector* px = dynamic_cast<const Epetra_Vector*>(&x);
      Epetra_Vector* py = dynamic_cast<Epetra_Vector*>(&y);
      if (px == NULL || py == NULL) EPETRA_ESI_ERR_BEHAVIOR(-2);
      err = epetra_crsmatrix_->Multiply(false, *px, *py);
      if (err != 0) EPETRA_ESI_ERR_BEHAVIOR(-3);
      return( 0 ); }


  //
  //esi::OperatorTranspose
  //

  /** Function for applying the transpose of this operator to a vector.
    (Transpose matvec.)
   @param x Input. The run-time type of this vector is required to be 
      epetra_esi::Vector.
   @param y On exit, the contents of this will be the result of the matvec.
       The run-time of this vector must also be epetra_esi::Vector.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode applyTranspose(esi::Vector<Scalar, Ordinal>& x,
                                        esi::Vector<Scalar, Ordinal>& y)
    { int err = 0;
      if (!setupHasBeenCalled_) err = setup();
      if (err) EPETRA_ESI_ERR_BEHAVIOR(-1);
      const Epetra_Vector* px = dynamic_cast<const Epetra_Vector*>(&x);
      Epetra_Vector* py = dynamic_cast<Epetra_Vector*>(&y);
      if (px == NULL || py == NULL) EPETRA_ESI_ERR_BEHAVIOR(-1);
      EPETRA_ESI_ERR_BEHAVIOR( epetra_crsmatrix_->Multiply(true, *px, *py) ); }


  //
  //esi::MatrixData
  //

  /** Query this operator's global dimensions.
    @param rows Output.
    @param columns Output. 
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getGlobalSizes(Ordinal& rows, Ordinal& columns)
    { rows = epetra_crsmatrix_->NumGlobalRows();
      columns = epetra_crsmatrix_->NumGlobalCols(); return(0); }

  /** Query this operator's local dimensions.
    @param rows Output.
    @param columns Output.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getLocalSizes(Ordinal& rows, Ordinal& columns)
    { rows = epetra_crsmatrix_->NumMyRows();
      columns = epetra_crsmatrix_->NumMyCols(); return(0); }


  /** Obtain the index-spaces that describe the row-space and column-space of
     this operator. !!!! Note that at the current time, the same index-space
     is returned
   for both the row-space and the column-space. This is clearly wrong. Our
   implementation always stores complete rows locally, so the column-space
   partitioning information is not correctly represented by the row-space. !!!!
    @param rowIndexSpace Output.
    @param colIndexSpace Output.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getIndexSpaces(esi::IndexSpace<Ordinal>*& rowIndexSpace, esi::IndexSpace<Ordinal>*& colIndexSpace)
    { Epetra_Map& rmap = (Epetra_Map&)(epetra_crsmatrix_->RowMap());
      Epetra_Map& cmap = (Epetra_Map&)(epetra_crsmatrix_->RowMap());
      rowIndexSpace = dynamic_cast<esi::IndexSpace<Ordinal>*>(&rmap);
      colIndexSpace = dynamic_cast<esi::IndexSpace<Ordinal>*>(&cmap);
      if (rowIndexSpace==NULL || colIndexSpace==NULL) {
        EPETRA_ESI_ERR_BEHAVIOR(-1);
      }
      return(0); 
    }


  //
  //esi::MatrixRowWriteAccess
  //

  /** Query the allocated length of a row.
    @param row Input, global row-number.
    @param length Output. Amount of space allocated for this row.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getRowAllocatedLength(Ordinal row, Ordinal& length)
    { length = epetra_crsmatrix_->NumAllocatedGlobalEntries(row);
      if (length<0) EPETRA_ESI_ERR_BEHAVIOR(-1);
      return(0); }


  /** Get the number of non-zeros currently stored in a row.
    @param row Input, global row-number.
    @param numEntries Output, number of non-zeros.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getRowNonzeros(Ordinal row, Ordinal& numEntries)
    { numEntries = epetra_crsmatrix_->NumGlobalEntries(row);
      if (numEntries<0) EPETRA_ESI_ERR_BEHAVIOR(-1); 
      return(0); }


  /** Set all values of the matrix to be the specified scalar.
    @param coef Input. The value to be set throughout the matrix.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode setAllValues(Scalar coef)
    { EPETRA_ESI_ERR_BEHAVIOR( epetra_crsmatrix_->PutScalar(coef) ); }

  /** Copy a set of coefficients and column-indices into a row.
    @param row Input, global row-number.
    @param coefs Input, coefficients.
    @param colIndices Input, column-indices.
    @param length Input, number of coefs and indices being copied in.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode copyIntoRow(Ordinal row, Scalar* coefs,
                                     Ordinal* colIndices, Ordinal length)
    {
      if (whichConstructor_ == 0) {
        EPETRA_ESI_ERR_BEHAVIOR( epetra_crsmatrix_->ReplaceGlobalValues(row, length, coefs, colIndices) );
      }
      else {
        EPETRA_ESI_ERR_BEHAVIOR( epetra_crsmatrix_->InsertGlobalValues(row, length, coefs, colIndices) );
      }
     }


  /** Sum (accumulate) a set of coefficients and column-indices into a row.
    Note that the specified column indices must already exist in the structure
    of the matrix. i.e., they must already have been inserted either in the
    Epetra_CrsGraph object provided at construction time, or if no graph was
    provided at construction time then they must already have been inserted
    via a call to the 'copyIntoRow' function.
    @param row Input, global row-number.
    @param coefs Input, coefficients.
    @param colIndices Input, column-indices.
    @param length Input, number of coefs and indices being copied in.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode sumIntoRow(Ordinal row, Scalar* coefs, Ordinal* colIndices, Ordinal length)
    { EPETRA_ESI_ERR_BEHAVIOR( epetra_crsmatrix_->SumIntoGlobalValues(row, length, coefs, colIndices) ); }


  /** Query whether the loadComplete function has been called.
    @param state Output.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode isLoaded(bool& state)
    { state = epetra_crsmatrix_->Filled(); return(0); }


  /** Query whether the internal data structures have been allocated.
    @param state Output. Always true, since epetra_esi::CrsMatrix uses the
       Epetra_CrsMatrix constructor that internally allocates at 
       construction time.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode isAllocated(bool& state)
    { state = true; return(0); }


  /** Signal the matrix that data-loading is complete, and that any necessary
    synchronization and/or data-consolidation can now be performed.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode loadComplete( void )
    { EPETRA_ESI_ERR_BEHAVIOR(setup()); }


  /** Allocate the matrix with the given list of row-lengths. This function is
     not available (always returns -1). epetra_esi::CrsMatrix is either
    allocated at construction time (if a graph is provided), or is allocated
    gradually by calls to the 'copyIntoRow' method.
    @param rowLengths Input, list of lengths of the local rows.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode allocate(Ordinal* rowLengths)
    { (void)rowLengths; EPETRA_ESI_ERR_BEHAVIOR(-1); }


  /** Allocate the matrix with all rows the same length. This function is not
    available (always returns -1). epetra_esi::CrsMatrix is either allocated at
    construction time (if a graph is provided), or is allocated gradually by
    calls to the 'copyIntoRow' method.
    @param rowLengths Input, length to be given to all rows.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode allocateRowsSameLength(Ordinal rowLengths)
    { (void)rowLengths; EPETRA_ESI_ERR_BEHAVIOR(-1); }


  /** Set the length of a particular row. This function is not available
    (always returns -1). epetra_esi::CrsMatrix is either allocated at
    construction time (if a graph is provided), or is allocated gradually by
    calls to the 'copyIntoRow' method.
    @param row Input, global row-number.
    @param length, Input, the length of the row.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode setRowLength(Ordinal row, Ordinal length)
    { (void)row; (void)length; EPETRA_ESI_ERR_BEHAVIOR(-1); }


  //
  //esi::MatrixRowReadAccess
  //

  /** Extract a copy of a row of the matrix (into caller-allocated memory).
    @param row Input, global row-number.
    @param coefs Input, coefficients.
    @param colIndices Input, column-indices.
    @param length Input, length of caller-allocated lists (coefs and colIndices)
    @param rowLength Output, Length of this row of the matrix.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode copyOutRow(Ordinal row, Scalar* coefs, Ordinal* colIndices, Ordinal length, Ordinal& rowLength)
    { EPETRA_ESI_ERR_BEHAVIOR(epetra_crsmatrix_->ExtractGlobalRowCopy(row, length, rowLength, coefs,colIndices));}


  /** Extract a copy of only the coefficients for a row of the matrix (into
     caller-allocated memory).
    @param row Input, global row-number.
    @param coefs Input, coefficients.
    @param length Input, length of caller-allocated list (coefs).
    @param rowLength Output, Length of this row of the matrix.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode copyOutRowCoefficients(Ordinal row, Scalar* coefs, Ordinal length, Ordinal& rowLength)
    { EPETRA_ESI_ERR_BEHAVIOR(epetra_crsmatrix_->ExtractGlobalRowCopy(row, length, rowLength, coefs)); }


  /** Extract a copy of only the column-indices for a row of the matrix (into
     caller-allocated memory).
    @param row Input, global row-number.
    @param colIndices Output, column-indices.
    @param length Input, length of caller-allocated list (coefs).
    @param rowLength Output, Length of this row of the matrix.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode copyOutRowIndices(Ordinal row, Ordinal* colIndices, Ordinal length, Ordinal& rowLength)
    { EPETRA_ESI_ERR_BEHAVIOR(epetra_crsmatrix_->Graph().ExtractGlobalRowCopy(row, length,rowLength,colIndices));}

  Epetra_CrsMatrix* getEpetra_CrsMatrix() { return(epetra_crsmatrix_); }

 private:
  bool setupHasBeenCalled_;
#ifdef EPETRA_MPI
  MPI_Comm* mpicommPtr_;
  MPI_Comm mpicomm_;
#endif
  int err1_, err2_, err3_;
  epetra_esi::IndexSpace<Ordinal>* ispace_;
  Ordinal globalSize_, localSize_;
  Epetra_Comm* comm_;
  Epetra_Map* rowMap_;
  Epetra_CrsMatrix* epetra_crsmatrix_;
  int whichConstructor_;
};

}; //namespace epetra_esi

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_CrsMatrix.cpp"
#endif

#endif

