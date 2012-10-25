/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_MatrixReducer_hpp_
#define _fei_MatrixReducer_hpp_

#include <fei_iosfwd.hpp>
#include <fei_mpi.h>
#include <fei_defs.h>

#include <fei_Matrix.hpp>
#include <fei_Reducer.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix_core.hpp>

#undef fei_file
#define fei_file "fei_MatrixReducer.hpp"
#include <fei_ErrMacros.hpp>

namespace fei {

  class MatrixReducer : public fei::Matrix {
  public:
    /** Constructor */
    MatrixReducer(fei::SharedPtr<fei::Reducer> reducer,
                  fei::SharedPtr<fei::Matrix> target);

    /** Destructor */
    virtual ~MatrixReducer();

    /** Query for the underlying target matrix. */
    fei::SharedPtr<fei::Matrix> getTargetMatrix()
      { return(target_); }

    /** Return a name describing the run-time type
	of this object.
    */
    const char* typeName() { return(target_->typeName()); }

    /** Parameters method
     */
    int parameters(const fei::ParameterSet& paramset);

    fei::SharedPtr<fei::MatrixGraph> getMatrixGraph() const
      {return( target_->getMatrixGraph() ); }

    /** Set the fei::MatrixGraph associated with this matrix */
    void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph);

    /** Get the global number of rows in the matrix.
     */
    int getGlobalNumRows() const;

    /** Get the local number of rows in the matrix.
     */
    int getLocalNumRows() const;

    /** Set a specified scalar throughout the matrix. */
    int putScalar(double scalar);

   /** Get the length of a row of the matrix.
       @param row Global 0-based equation number
       @param length Output. Length of the row.
       @return error-code non-zero if any error occurs.
   */
    int getRowLength(int row, int& length) const;

   /** Obtain a copy of the coefficients and indices for a row of the matrix.
       @param row Global 0-based equation number
       @param len Length of the caller-allocated coefs and indices arrays
       @param coefs Caller-allocated array, length 'len', to be filled with
       coefficients
       @param indices Caller-allocated array, length 'len', to be filled with
       indices. (These indices will be global 0-based equation numbers.)
       @return error-code non-zero if any error occurs.
   */
    int copyOutRow(int row, int len, double* coefs, int* indices) const;

    /** Sum coefficients into the matrix, adding them to any coefficients that
	may already exist at the specified row/column locations.

	@param numRows
	@param rows
	@param numCols
	@param cols
	@param values
	@param format For compatibility with old FEI elemFormat...
	0 means row-wise or row-major, 3 means column-major. Others not recognized
     */
    int sumIn(int numRows, const int* rows,
	      int numCols, const int* cols,
	      const double* const* values,
	      int format=0);

    /** Copy coefficients into the matrix, overwriting any coefficients that
	may already exist at the specified row/column locations.

	@param numRows
	@param rows
	@param numCols
	@param cols
	@param values
	@param format For compatibility with old FEI elemFormat...
	0 means row-wise or row-major, 3 means column-major. Others not recognized
    */
    int copyIn(int numRows, const int* rows,
	       int numCols, const int* cols,
	       const double* const* values,
	       int format=0);

    /** Sum coefficients into the matrix, specifying row/column locations by
	identifier/fieldID pairs.
	@param fieldID Input. field-identifier for which data is being input.
	@param idType Input. The identifier-type of the identifiers.
	@param rowID Input. Identifier in row-space, for which data is being
	input.
	@param colID Input. Identifier in column-space, for which data is being
	input.
	@param data Input. C-style table of data. num-rows is the field-size
	(i.e., number of scalar components that make up the field) of 'fieldID',
	as is num-columns.
	@param format For compatibility with old FEI elemFormat...
	0 means row-wise or row-major, 3 means column-major. Others not recognized
	@return error-code 0 if successful
    */
    int sumInFieldData(int fieldID,
		       int idType,
		       int rowID,
		       int colID,
		       const double* const* data,
		       int format=0);

    /** Sum coefficients into the matrix, specifying row/column locations by
	identifier/fieldID pairs.
	@param fieldID Input. field-identifier for which data is being input.
	@param idType Input. The identifier-type of the identifiers.
	@param rowID Input. Identifier in row-space, for which data is being
	input.
	@param colID Input. Identifier in column-space, for which data is being
	input.
	@param data Input. 1-D list representing a packed table of data. Data may
	be backed in row-major or column-major order and this may be specified with
	the 'format' argument. The "table" of data is of size num-rows X num-columns
	and num-rows is the field-size (i.e., number of scalar components that
	make up the field) of 'fieldID', as is num-columns.
	@param format For compatibility with old FEI elemFormat...
	0 means row-wise or row-major, 3 means column-major. Others not recognized
	@return error-code 0 if successful
    */
    int sumInFieldData(int fieldID,
		       int idType,
		       int rowID,
		       int colID,
		       const double* data,
		       int format=0);

    /** Sum coefficients, associated with a connectivity-block that was
	initialized on the MatrixGraph object, into this matrix.

	@param blockID
	@param connectivityID
	@param values
	@param format For compatibility with old FEI elemFormat...
	0 means row-wise or row-major, 3 means column-major. Others not recognized
     */
    int sumIn(int blockID, int connectivityID,
	      const double* const* values,
	      int format=0);

    /** Perform any necessary internal communications/synchronizations or other
	operations appropriate at end of data input. For some implementations this
	will be a no-op.
    */
    int globalAssemble();

    /** Form a matrix-vector product y = 'this' * x
     */
    int multiply(fei::Vector* x,
		 fei::Vector* y);

    void setCommSizes() { target_->setCommSizes(); }

    /** After local overlapping data has been input, (e.g., element-data for a
	finite-element application) call this method to have data that 
	corresponds to shared identifiers be communicated from sharing-but-not-
	owning processors, to owning processors.
    */
    int gatherFromOverlap(bool accumulate = true);

    /** Implementation of fei::Matrix::writeToFile */
    int writeToFile(const char* filename,
		    bool matrixMarketFormat=true);

    /** Implementation of fei::Matrix::writeToStream */

    int writeToStream(FEI_OSTREAM& ostrm,
		      bool matrixMarketFormat=true);

    bool usingBlockEntryStorage()
      { return(target_->usingBlockEntryStorage()); }

    /** for experts only */
    int giveToUnderlyingMatrix(int numRows, const int* rows,
			       int numCols, const int* cols,
			       const double* const* values,
			       bool sumInto,
			       int format);

    /** for experts only */
    int giveToUnderlyingBlockMatrix(int row,
				    int rowDim,
				    int numCols,
				    const int* cols,
				    const int* LDAs,
				    const int* colDims,
				    const double* const* values,
				    bool sumInto);

    void markState();

    bool changedSinceMark();

  private:
    int giveToMatrix(int numRows, const int* rows,
		     int numCols, const int* cols,
		     const double* const* values,
		     bool sumInto,
		     int format);
 
    int giveToBlockMatrix(int numRows, const int* rows,
			  int numCols, const int* cols,
			  const double* const* values,
			  bool sumInto);

    fei::SharedPtr<fei::Reducer> reducer_;
    fei::SharedPtr<fei::Matrix> target_;
    bool globalAssembleCalled_;
    bool changedSinceMark_;
  };//class MatrixReducer
}//namespace fei

#endif

