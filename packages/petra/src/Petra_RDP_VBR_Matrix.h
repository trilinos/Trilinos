#ifndef _PETRA_RDP_VBR_MATRIX_H_
#define _PETRA_RDP_VBR_MATRIX_H_
//! Petra_RDP_VBR_Matrix: A class for constructing and using real-valued double-precision sparse compressed row matrices.

/*! The Petra_RDP_VBR_Matrix enable the piecewise construction and use of real-valued double-precision sparse matrices
    where matrix entries are intended for row access.

    At this time, the primary function provided by Petra_RDP_VBR_Matrix is matrix time vector and matrix 
    times multi-vector multiplication.  It is also possible to extract matrix rows from a constructed matrix.

<b>Constructing Petra_RDP_VBR_Matrix objects</b>

Constructing Petra_RDP_VBR_Matrix objects is a multi-step process.  The basic steps are as follows:
<ol>
  <li> Create Petra_RDP_VBR_Matrix instance, including storage,  via constructor.
  <li> Enter values via one or more Put or SumInto functions.
  <li> Complete construction via FillComplete call.
</ol>

Note that, even after a matrix is constructed, it is possible to update existing matrix entries.  It is \e not possible to
create new entries.

<b> Counting Floating Point Operations </b>

Each Petra_RDP_VBR_Matrix object keep track of the number
of \e serial floating point operations performed using the specified object as the \e this argument
to the function.  The Flops() function returns this number as a double precision number.  Using this 
information, in conjunction with the Petra_Time class, one can get accurate parallel performance
numbers.  The ResetFlops() function resets the floating point counter.

\warning A Petra_BlockMap is required for the Petra_RDP_VBR_Matrix constructor.

*/    

#include "Petra_Petra.h" 
#include "Petra_Flops.h"
#include "Petra_BLAS.h"
#include "Petra_BlockMap.h"
#include "Petra_Import.h"
#include "Petra_Export.h"
#include "Petra_CRS_Graph.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_MultiVector.h"
#include "Petra_RDP_RowMatrix.h"


class Petra_RDP_VBR_Matrix: public Petra_Flops, public Petra_BLAS, public virtual Petra_RDP_RowMatrix{
      
  // Give ostream << function some access to private and protected data/functions.

  friend ostream& operator << (ostream& os, const Petra_RDP_VBR_Matrix& A);

 public:
  //! Petra_RDP_VBR_Matrix constuctor with variable number of indices per row.
  /*! Creates a Petra_RDP_VBR_Matrix object and allocates storage.  
    
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Petra_BlockMap.
    \param In
           NumBlockEntriesPerRow - An integer array of length NumRows
	   such that NumBlockEntriesPerRow[i] indicates the (approximate) number of Block entries in the ith row.
  */
  Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_BlockMap& RowMap, int *NumBlockEntriesPerRow);
  
  //! Petra_RDP_VBR_Matrix constuctor with fixed number of indices per row.
  /*! Creates a Petra_RDP_VBR_Matrix object and allocates storage.  
    
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Petra_BlockMap.
    \param In
           NumBlockEntriesPerRow - An integer that indicates the (approximate) number of Block entries in the each Block row.
	   Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	   
  */
  Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_BlockMap& RowMap, int NumBlockEntriesPerRow);

  //! Petra_RDP_VBR_Matrix constuctor with variable number of indices per row.
  /*! Creates a Petra_RDP_VBR_Matrix object and allocates storage.  
    
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Petra_BlockMap.
    \param In 
           ColMap - A Petra_BlockMap.
    \param In
           NumBlockEntriesPerRow - An integer array of length NumRows
	   such that NumBlockEntriesPerRow[i] indicates the (approximate) number of entries in the ith row.
  */
  Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_BlockMap& RowMap, const Petra_BlockMap& ColMap, int *NumBlockEntriesPerRow);
  
  //! Petra_RDP_VBR_Matrix constuctor with fixed number of indices per row.
  /*! Creates a Petra_RDP_VBR_Matrix object and allocates storage.  
    
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In 
           RowMap - A Petra_BlockMap.
    \param In 
           ColMap - A Petra_BlockMap.
    \param In
           NumBlockEntriesPerRow - An integer that indicates the (approximate) number of entries in the each row.
	   Note that it is possible to use 0 for this value and let fill occur during the insertion phase.
	   
  */
  Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_BlockMap& RowMap, const Petra_BlockMap& ColMap, int NumBlockEntriesPerRow);

  //! Construct a matrix using an existing Petra_CRS_Graph object.
  /*! Allows the nonzero structure from another matrix, or a structure that was
      constructed independently, to be used for this matrix.
    \param In
           CV - A Petra_DataAccess enumerated type set to Copy or View.
    \param In
           Graph - A Petra_CRS_Graph object, extracted from another Petra matrix object or constructed directly from
	   using the Petra_CRS_Graph constructors.
  */

  Petra_RDP_VBR_Matrix(Petra_DataAccess CV, const Petra_CRS_Graph & Graph);

  //! Copy constructor.
  Petra_RDP_VBR_Matrix(const Petra_RDP_VBR_Matrix & Matrix);

  //! Petra_RDP_VBR_Matrix Destructor
  virtual ~Petra_RDP_VBR_Matrix();

  //! Initialize all values in graph of the matrix with constant value.
  /*!
    \param In
           Scalar - Value to use.

    \return Integer error code, set to 0 if successful.
  */
    int PutScalar(double Scalar);

  //! Initiate insertion of a list of elements in a given global row of the matrix, values are inserted via SubmitEntry().
  /*!
    \param In
           BlockRow - Block Row number (in global coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginInsertGlobalValues(int BlockRow, int NumBlockEntries, int * BlockIndices);

  //! Initiate insertion of a list of elements in a given local row of the matrix, values are inserted via SubmitEntry().
  /*!
    \param In
           BlockRow - Block Row number (in local coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginInsertMyValues(int BlockRow, int NumBlockEntries, int * BlockIndices);

  //! Initiate replacement of current values with this list of entries for a given global row of the matrix, values are replaced via SubmitEntry()
  /*!
    \param In
           Row - Block Row number (in global coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginReplaceGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices);

  //! Initiate replacement of current values with this list of entries for a given local row of the matrix, values are replaced via SubmitEntry()
  /*!
    \param In
           Row - Block Row number (in local coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginReplaceMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices);

  //! Initiate summing into current values with this list of entries for a given global row of the matrix, values are replaced via SubmitEntry()
  /*!
    \param In
           Row - Block Row number (in global coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Global column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginSumIntoGlobalValues(int BlockRow, int NumBlockEntries, int *BlockIndices);

  //! Initiate summing into current values with this list of entries for a given local row of the matrix, values are replaced via SubmitEntry()
  /*!
    \param In
           Row - Block Row number (in local coordinates) to put elements.
    \param In
           NumBlockEntries - Number of entries.
    \param In
           Indices - Local column indices corresponding to values.

    \return Integer error code, set to 0 if successful.
  */
    int BeginSumIntoMyValues(int BlockRow, int NumBlockEntries, int *BlockIndices);

    //! Submit a block entry to the indicated block row and column specified in the Begin routine.
    /* Submit a block entry that will recorded in the block row that was initiated by one of the
       Begin routines listed above.  Once a one of the following routines: BeginInsertGlobalValues(),
       BeginInsertMyValues(), BeginReplaceGlobalValues(), BeginReplaceMyValues(), BeginSumIntoGlobalValues(),
       BeginSumIntoMyValues(), you \e must call SubmitBlockEntry() NumBlockEntries times to register the values 
       corresponding to the block indices passed in to the Begin routine.  If the Petra_RDP_VBR_Matrix constuctor
       was called in Copy mode, the values will be copied.  However, no copying will be done until the EndSubmitEntries()
       function is call to complete submission of the current block row.  If the constructor was called in View mode, all
       block entries passed via SubmitBlockEntry() will not be copied, but a pointer will be set to point to the argument Values
       that was passed in by the user.

       For performance reasons, SubmitBlockEntry() does minimal processing of data.  Any processing that can be
       delayed is performed in EndSubmitEntries().

    \param In
           Values - The starting address of the values.
    \param In
           LDA - The stride between successive columns of Values.
    \param In
           NumRows - The number of rows passed in.
    \param In
           NumCols - The number of columns passed in.

    \return Integer error code, set to 0 if successful.
    */
    int SubmitBlockEntry(double *Values, int LDA, int NumRows, int NumCols);

    //! Completes processing of all data passed in for the current block row.
    /*! This function completes the processing of all block entries submitted via SubmitBlockEntry().  
        It also checks to make sure that SubmitBlockEntry was called the correct number of times as
	specified by the Begin routine that initiated the entry process.
    */

    int EndSubmitEntries();

    //! Signal that data entry is complete, perform transformations to local index space.
    /* This version of TransformToLocal assumes that the domain and range distributions are
       identical to the matrix row distributions.
    */
    int TransformToLocal();

    //! Signal that data entry is complete, perform transformations to local index space.
    /* This version of TransformToLocal requires the explicit specification of the domain
       and range distribution maps.  These maps are used for importing and exporting vector
       and multi-vector elements that are needed for distributed matrix computations.  For
       example, to compute y = Ax in parallel, we would specify the DomainMap as the distribution
       of the vector x and the RangeMap as the distribution of the vector y.
    \param In
           DomainMap - Map that describes the distribution of vector and multi-vectors in the
	               matrix domain.
    \param In
           RangeMap - Map that describes the distribution of vector and multi-vectors in the
	               matrix range.
    */
    int TransformToLocal(Petra_BlockMap *DomainMap, Petra_BlockMap *RangeMap);

    //! If FillComplete() has been called, this query returns true, otherwise it returns false.
    bool Filled() const {return(Graph_->Filled());};

    // Matrix data extraction routines

    //! Copy the block indices into user-provided array, set pointers for rest of data for specified global block row.
    /*! 
      This function provides the lightest weight approach to accessing a global block row when the matrix may be
      be stored in local or global index space.  In other words, this function will always work because the block
      indices are returned in user-provided space.  All other array arguments are independent of whether or not
      indices are local or global.  Other than the BlockIndices array, all other array argument are returned as 
      pointers to internal data.

    \param In
           BlockRow - Global block row to extract.
    \param In
	   MaxNumBlockEntries - Length of user-provided BlockIndices array.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries actually extracted.
    \param Out
	   BlockIndices - Extracted global column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that pointed to by Values.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that is pointed to by Values.
    \param Out
	   Values - Pointer to list of pointers to block entries. Note that the actual values are not copied.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractGlobalBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
				      int & RowDim,  int & NumBlockEntries, 
				      int * BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const;

    //! Copy the block indices into user-provided array, set pointers for rest of data for specified local block row.
    /*! 
      This function provides the lightest weight approach to accessing a local block row when the matrix may be
      be stored in local or global index space.  In other words, this function will always work because the block
      indices are returned in user-provided space.  All other array arguments are independent of whether or not
      indices are local or global.  Other than the BlockIndices array, all other array argument are returned as 
      pointers to internal data.

    \param In
           BlockRow - Local block row to extract.
    \param In
	   MaxNumBlockEntries - Length of user-provided BlockIndices array.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries actually extracted.
    \param Out
	   BlockIndices - Extracted local column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that pointed to by Values.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that is pointed to by Values.
    \param Out
	   Values - Pointer to list of pointers to block entries. Note that the actual values are not copied.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
				       int & RowDim, int & NumBlockEntries, 
				       int * BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const;

    //! Initiates a copy of the specified global row in user-provided arrays.
    /*! 
    \param In
           BlockRow - Global block row to extract.
    \param In
	   MaxNumBlockEntries - Length of user-provided BlockIndices, ColDims, and LDAs arrays.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries actually extracted.
    \param Out
	   BlockIndices - Extracted global column indices for the corresponding block entries.
    \param Out
	   ColDim - List of column dimensions for each corresponding block entry that will be extracted.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractGlobalBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
				       int & RowDim,  int & NumBlockEntries, 
				       int * BlockIndices, int * ColDims) const;

    //! Initiates a copy of the specified local row in user-provided arrays.
    /*! 
    \param In
           BlockRow - Local block row to extract.
    \param In
	   MaxNumBlockEntries - Length of user-provided BlockIndices, ColDims, and LDAs arrays.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries actually extracted.
    \param Out
	   BlockIndices - Extracted local column indices for the corresponding block entries.
    \param Out
	   ColDim - List of column dimensions for each corresponding block entry that will be extracted.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractMyBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
				       int & RowDim, int & NumBlockEntries, 
				       int * BlockIndices, int * ColDims) const;

    //! Extract a copy of an entry from the block row specified by one of the BeginExtract routines.
    /*! Once BeginExtractGlobalBlockRowCopy() or BeginExtractMyBlockRowCopy() is called, you can extract
        the block entries of specified block row one-entry-at-a-time.  The entries will be extracted
	in an order corresponding to the BlockIndices list that was returned by the BeginExtract routine.

    \param In
           SizeOfValues - Amount of memory associated with Values.  This must be at least as big as
	                  LDA*NumCol, where NumCol is the column dimension of the block entry being copied
    \param InOut
           Values - Starting location where the block entry will be copied.  
    \param In
           LDA - Specifies the stride that will be used when copying columns into Values.
    \param In
           SumInto - If set to true, the block entry values will be summed into existing values.
    */

    int ExtractEntryCopy(int SizeOfValues, double * Values, int LDA, bool SumInto) const;

    //! Initiates a view of the specified global row, only works if matrix indices are in global mode.
    /*! 
    \param In
           BlockRow - Global block row to view.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries to be viewed.
    \param Out
	   BlockIndices - Pointer to global column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that will be viewed.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that will be viewed.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
				       int * & BlockIndices, int * & ColDims, int * & LDAs) const;

    //! Initiates a view of the specified local row, only works if matrix indices are in local mode.
    /*! 
    \param In
           BlockRow - Local block row to view.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries to be viewed.
    \param Out
	   BlockIndices - Pointer to local column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that will be viewed.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that will be viewed.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
				       int * & BlockIndices, int * & ColDims, int * & LDAs) const;


    //! Returns a pointer to, and leading dimension o, the current block entry.
    /*! After a call to BeginExtractGlobal() or BlockRowViewBeginExtractMyBlockRowView(),
        ExtractEntryView() can be called up to NumBlockEntries times to get pointer and stride
	information for each block entry in the specified block row.
    \param InOut
           Values - A pointer that will be set to point to the starting address of the current block entry.
    */
    
    int ExtractEntryView(double * & Values) const;

    //! Initiates a view of the specified global row, only works if matrix indices are in global mode.
    /*! 
    \param In
           BlockRow - Global block row to view.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries to be viewed.
    \param Out
	   BlockIndices - Pointer to global column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that will be viewed.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that will be viewed.
    \param Out
	   Values - Pointer to an array of pointers to the block entries in the specified block row.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractGlobalBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
				  int * & BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const;

    //! Initiates a view of the specified local row, only works if matrix indices are in local mode.
    /*! 
    \param In
           BlockRow - Local block row to view.
    \param Out
	   RowDim - Number of equations in the requested block row.
    \param Out
	   NumBlockEntries - Number of nonzero entries to be viewed.
    \param Out
	   BlockIndices - Pointer to local column indices for the corresponding block entries.
    \param Out
	   ColDim - Pointer to list of column dimensions for each corresponding block entry that will be viewed.
    \param Out
	   LDAs - Pointer to list of leading dimensions for each corresponding block entry that will be viewed.
    \param Out
	   Values - Pointer to an array of pointers to the block entries in the specified block row.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
			      int * & BlockIndices, int * & ColDims, int * & LDAs, double ** & Values) const;


    //! Returns a copy of the main diagonal in a user-provided vector.
    /*! 
    \param Out
	   Diagonal - Extracted main diagonal.

    \return Integer error code, set to 0 if successful.
  */
    int ExtractDiagonalCopy(Petra_RDP_Vector & Diagonal) const;

    //! Initiates a copy of the block diagonal entries to user-provided arrays.
    /*! 
    \param In
	   MaxNumBlockDiagonalEntries - Length of user-provided RowColDims array.
    \param Out
	   NumBlockDiagonalEntries - Number of block diagonal entries that can actually be extracted.
    \param Out
	   RowColDim - List  of row and column dimension for corresponding block diagonal entries.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractBlockDiagonalCopy(int MaxNumBlockDiagonalEntries, 
				      int & NumBlockDiagonalEntries, int * RowColDims ) const;
    //! Extract a copy of a block diagonal entry from the matrix.
    /*! Once BeginExtractBlockDiagonalCopy() is called, you can extract
        the block diagonal entries one-entry-at-a-time.  The entries will be extracted
	in ascending order.

    \param In
           SizeOfValues - Amount of memory associated with Values.  This must be at least as big as
	                  LDA*NumCol, where NumCol is the column dimension of the block entry being copied
    \param InOut
           Values - Starting location where the block entry will be copied.  
    \param In
           LDA - Specifies the stride that will be used when copying columns into Values.
    \param In
           SumInto - If set to true, the block entry values will be summed into existing values.
    */

    int ExtractBlockDiagonalEntryCopy(int SizeOfValues, double * Values, int LDA, bool SumInto) const;

    //! Initiates a view of the block diagonal entries.
    /*! 
    \param Out
	   NumBlockDiagonalEntries - Number of block diagonal entries that can be viewed.
    \param Out
	   RowColDim - Pointer to list  of row and column dimension for corresponding block diagonal entries.
	  
    \return Integer error code, set to 0 if successful.
  */
    int BeginExtractBlockDiagonalView(int & NumBlockDiagonalEntries, int * & RowColDims ) const;

    //! Extract a view of a block diagonal entry from the matrix.
    /*! Once BeginExtractBlockDiagonalView() is called, you can extract a view of
        the block diagonal entries one-entry-at-a-time.  The views will be extracted
	in ascending order.

    \param Out
           Values - Pointer to internal copy of block entry.  
    \param Out
           LDA - Column stride of Values.
    */

    int ExtractBlockDiagonalEntryView(double * & Values, int & LDA) const;

    // Mathematical functions.


    //! Returns the result of a Petra_RDP_VBR_Matrix multiplied by a Petra_RDP_Vector x in y.
    /*! 
    \param In
	   TransA - If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   x - A Petra_RDP_Vector to multiply by.
    \param Out
	   y - A Petra_RDP_Vector containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply1(bool TransA, const Petra_RDP_Vector& x, Petra_RDP_Vector& y) const;
    //! Returns a copy of the specified global row in user-provided arrays.
    /*! 
    \param In
           GlobalRow - Global row to extract.
    \param In
	   Length - Length of Values and Indices.
    \param Out
	   NumEntries - Number of nonzero entries extracted.
    \param Out
	   Values - Extracted values for this row.
    \param Out
	   Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractGlobalRowCopy(int GlobalRow, int Length, int & NumEntries, double *Values, int * Indices) const;

    //! Returns a copy of the specified local row in user-provided arrays.
    /*! 
    \param In
           MyRow - Local row to extract.
    \param In
	   Length - Length of Values and Indices.
    \param Out
	   NumEntries - Number of nonzero entries extracted.
    \param Out
	   Values - Extracted values for this row.
    \param Out
	   Indices - Extracted global column indices for the corresponding values.
	  
    \return Integer error code, set to 0 if successful.
  */
    int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;


    //! Returns the result of a Petra_RDP_VBR_Matrix multiplied by a Petra_RDP_MultiVector X in Y.
    /*! 
    \param In
	   TransA -If true, multiply by the transpose of matrix, otherwise just use matrix.
    \param In
	   X - A Petra_RDP_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Petra_RDP_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
    int Multiply(bool TransA, const Petra_RDP_MultiVector& X, Petra_RDP_MultiVector& Y) const;

    //! Returns the result of a solve using the Petra_RDP_VBR_Matrix on a Petra_RDP_Vector x in y.
    /*! 
    \param In
	   Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
	   Trans -If true, solve transpose problem.
    \param In
	   UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
	   x -A Petra_RDP_Vector to solve for.
    \param Out
	   y -A Petra_RDP_Vector containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Petra_RDP_Vector& x, Petra_RDP_Vector& y) const;

    //! Returns the result of a Petra_RDP_VBR_Matrix multiplied by a Petra_RDP_MultiVector X in Y.
    /*! 
    \param In
	   Upper -If true, solve Ux = y, otherwise solve Lx = y.
    \param In
	   Trans -If true, solve transpose problem.
    \param In
	   UnitDiagonal -If true, assume diagonal is unit (whether it's stored or not).
    \param In
	   X - A Petra_RDP_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Petra_RDP_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
    int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Petra_RDP_MultiVector& X, Petra_RDP_MultiVector& Y) const;


    //! Computes the sum of absolute values of the rows of the Petra_RDP_VBR_Matrix, results returned in x.
    /*! The vector x will return such that x[i] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to LeftScale()
	will make the infinity norm of the resulting matrix exactly 1.
    \param Out
	   x -A Petra_RDP_Vector containing the row sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    int InvRowSums(Petra_RDP_Vector& x) const;

    //! Scales the Petra_RDP_VBR_Matrix on the left with a Petra_RDP_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(i)*A(i,j) where i denotes the row number of A
        and j denotes the column number of A.
    \param In
	   x -A Petra_RDP_Vector to solve for.

    \return Integer error code, set to 0 if successful.
  */
    int LeftScale(const Petra_RDP_Vector& x);

    //! Computes the sum of absolute values of the columns of the Petra_RDP_VBR_Matrix, results returned in x.
    /*! The vector x will return such that x[j] will contain the inverse of sum of the absolute values of the 
        \e this matrix will be sca such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.  Using the resulting vector from this function as input to 
	RighttScale() will make the one norm of the resulting matrix exactly 1.
    \param Out
	   x -A Petra_RDP_Vector containing the column sums of the \e this matrix. 
	   \warning It is assumed that the distribution of x is the same as the rows of \e this.

    \return Integer error code, set to 0 if successful.
  */
    int InvColSums(Petra_RDP_Vector& x) const ;

    //! Scales the Petra_RDP_VBR_Matrix on the right with a Petra_RDP_Vector x.
    /*! The \e this matrix will be scaled such that A(i,j) = x(j)*A(i,j) where i denotes the global row number of A
        and j denotes the global column number of A.
    \param In
	   x -The Petra_RDP_Vector used for scaling \e this.

    \return Integer error code, set to 0 if successful.
  */
    int RightScale(const Petra_RDP_Vector& x);
    // Atribute access functions

    //! Returns the infinity norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
    */ 
    double NormInf() const;

    //! Returns the one norm of the global matrix.
    /* Returns the quantity \f$ \| A \|_1\f$ such that
       \f[\| A \|_1 = \max_{1\lej\len} \sum_{i=1}^m |a_{ij}| \f].
    */ 
    double NormOne() const;


    //! Returns the maximum row dimension of all block entries on this processor.
    int MaxRowDim() const {return(Graph_->MaxRowDim());};

    //! Returns the maximum column dimension of all block entries on this processor.
    int MaxColDim() const {return(Graph_->MaxColDim());};

    //! Returns the maximum row dimension of all block entries across all processors.
    int GlobalMaxRowDim() const {return(Graph_->GlobalMaxRowDim());};

    //! Returns the maximum column dimension of all block entries across all processors.
    int GlobalMaxColDim() const {return(Graph_->GlobalMaxColDim());};

    //! Returns the number of matrix rows owned by the calling processor.
    int NumMyRows() const {return(Graph_->NumMyRows());};
    //! Returns the number of matrix columns owned by the calling processor.
    int NumMyCols() const {return(Graph_->NumMyCols());};

    //! Returns the number of nonzero entriesowned by the calling processor .
    int NumMyNonzeros() const {return(Graph_->NumMyNonzeros());};

    //! Returns the number of global matrix rows.
    int NumGlobalRows() const {return(Graph_->NumGlobalRows());};

    //! Returns the number of global matrix columns.
    int NumGlobalCols() const {return(Graph_->NumGlobalCols());};

    //! Returns the number of nonzero entries in the global matrix.
    int NumGlobalNonzeros() const {return(Graph_->NumGlobalNonzeros());};

    //! Returns the number of Block matrix rows owned by the calling processor.
    int NumMyBlockRows() const {return(Graph_->NumMyBlockRows());};

    //! Returns the number of Block matrix columns owned by the calling processor.
    int NumMyBlockCols() const {return(Graph_->NumMyBlockCols());};
    
    //! Returns the number of nonzero block entries in the calling processor's portion of the matrix.
    int NumMyBlockEntries() const {return(Graph_->NumMyEntries());};
    
    //! Returns the number of local nonzero block diagonal entries.
    int NumMyBlockDiagonals() const {return(Graph_->NumMyBlockDiagonals());};
    
    //! Returns the number of local nonzero diagonal entries.
    int NumMyDiagonals() const {return(Graph_->NumMyDiagonals());};
    
    //! Returns the number of global Block matrix rows.
    int NumGlobalBlockRows() const {return(Graph_->NumGlobalBlockRows());};
    
    //! Returns the number of global Block matrix columns.
    int NumGlobalBlockCols() const {return(Graph_->NumGlobalBlockCols());};
    
    //! Returns the number of nonzero block entries in the global matrix.
    int NumGlobalBlockEntries() const {return(Graph_->NumGlobalNonzeros());};
    
    //! Returns the number of global nonzero block diagonal entries.
    int NumGlobalBlockDiagonals() const {return(Graph_->NumGlobalBlockDiagonals());};
    
    //! Returns the number of global nonzero diagonal entries.
    int NumGlobalDiagonals() const {return(Graph_->NumGlobalDiagonals());};

    //! Returns the current number of nonzero Block entries in specified global row on this processor.
    int NumGlobalBlockEntries(int Row) const {return(Graph_->NumGlobalIndices(Row));};

    //! Returns the allocated number of nonzero Block entries in specified global row on this processor.
int NumAllocatedGlobalBlockEntries(int Row) const{return(Graph_->NumAllocatedGlobalIndices(Row));};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int MaxNumBlockEntries() const {return(Graph_->MaxNumIndices());};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
int GlobalMaxNumBlockEntries() const {return(Graph_->GlobalMaxNumIndices());};

    //! Returns the current number of nonzero Block entries in specified local row on this processor.
    int NumMyBlockEntries(int Row) const {return(Graph_->NumMyIndices(Row));};

    //! Returns the allocated number of nonzero Block entries in specified local row on this processor.
    int NumAllocatedMyBlockEntries(int Row) const {return(Graph_->NumAllocatedMyIndices(Row));};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int MaxNumNonzeros() const {return(Graph_->MaxNumNonzeros());};

    //! Returns the maximum number of nonzero entries across all rows on this processor.
    int GlobalMaxNumNonzeros() const {return(Graph_->GlobalMaxNumNonzeros());};

    //! Returns the index base for row and column indices for this graph.
    int IndexBase() const {return(Graph_->IndexBase());};

    //! Sort column entries, row-by-row, in ascending order.
    int SortEntries();

    //! If SortEntries() has been called, this query returns true, otherwise it returns false.
    bool Sorted() const {return(Graph_->Sorted());};

    //! Add entries that have the same column index. Remove redundant entries from list.
    int MergeRedundantEntries();

    //! If MergeRedundantEntries() has been called, this query returns true, otherwise it returns false.
    bool NoRedundancies() const {return(Graph_->NoRedundancies());};
    //! Eliminates memory that is used for construction.  Make consecutive row index sections contiguous.
    int OptimizeStorage();

    //! If OptimizeStorage() has been called, this query returns true, otherwise it returns false.
    bool StorageOptimized() const {return(Graph_->StorageOptimized());};

    //! If matrix indices has not been transformed to local, this query returns true, otherwise it returns false.
    bool IndicesAreGlobal() const {return(Graph_->IndicesAreGlobal());};

    //! If matrix indices has been transformed to local, this query returns true, otherwise it returns false.
    bool IndicesAreLocal() const {return(Graph_->IndicesAreLocal());};

    //! If matrix indices are packed into single array (done in OptimizeStorage()) return true, otherwise false.
    bool IndicesAreContiguous() const {return(Graph_->IndicesAreContiguous());};

    //! If matrix is lower triangular, this query returns true, otherwise it returns false.
    bool LowerTriangular() const {return(Graph_->LowerTriangular());};

    //! If matrix is upper triangular, this query returns true, otherwise it returns false.
    bool UpperTriangular() const {return(Graph_->UpperTriangular());};

    //! If matrix is lower triangular, this query returns true, otherwise it returns false.
    bool NoDiagonal() const {return(Graph_->NoDiagonal());};

    //! Returns the local row index for given global row index, returns -1 if no local row for this global row.
    int LRID( int GRID) const {return(Graph_->LRID(GRID));};

    //! Returns the global row index for give local row index, returns IndexBase-1 if we don't have this local row.
    int GRID( int LRID) const {return(Graph_->GRID(LRID));};

    //! Returns the local column index for given global column index, returns -1 if no local column for this global column.
    int LCID( int GCID) const {return(Graph_->LCID(GCID));};

    //! Returns the global column index for give local column index, returns IndexBase-1 if we don't have this local column.
    int GCID( int LCID) const {return(Graph_->GCID(LCID));};
 
    //! Returns true if the GRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyGRID(int GRID) const {return(Graph_->MyGRID(GRID));};
   
    //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyLRID(int LRID) const {return(Graph_->MyLRID(LRID));};

    //! Returns true if the GCID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyGCID(int GCID) const {return(Graph_->MyGCID(GCID));};
   
    //! Returns true if the LRID passed in belongs to the calling processor in this map, otherwise returns false.
    bool  MyLCID(int LCID) const {return(Graph_->MyLCID(LCID));};

    //! Returns true of GID is owned by the calling processor, otherwise it returns false.
    bool MyGlobalBlockRow(int GID) const {return(Graph_->MyGlobalRow(GID));};

    //! Returns a pointer to the Petra_CRS_Graph object associated with this matrix.
    const Petra_CRS_Graph & Graph() const {return(*Graph_);};

    //! Returns the Petra_BlockMap object associated with the rows of this matrix.
    const Petra_BlockMap & RowMap() const {return((Petra_BlockMap &)Graph_->RowMap());};


    //! Returns the Petra_BlockMap object associated with the rows of this matrix for RowMatrix interface.
    const Petra_BlockMap & BlockRowMap() const {return((Petra_BlockMap &)Graph_->RowMap());};

    //! Returns the Petra_BlockMap object associated with columns of this matrix.
    const Petra_BlockMap & ColMap() const {return((Petra_BlockMap &)Graph_->RowMap());};

    //! Returns the Petra_BlockMap object that describes the import vector for distributed operations.
    const Petra_BlockMap & ImportMap() const {return((Petra_BlockMap &) Graph_->ImportMap());};

    //! Returns the import vector map for distributed operations for the RowMatrix interface.
    const Petra_BlockMap & BlockImportMap() const {return((Petra_BlockMap &) Graph_->ImportMap());};

    //! Returns the Petra_Import object that contains the import operations for distributed operations.
    const Petra_Import * Importer() const {return(Graph_->Importer());};

    //! Returns the Petra_BlockMap object that describes the export vector for distributed operations.
    const Petra_BlockMap & ExportMap() const {return((Petra_BlockMap &) Graph_->ExportMap());};

    //! Returns the Petra_Export object that contains the export operations for distributed operations.
    const Petra_Export * Exporter() const {return(Graph_->Exporter());};

    //! Fills a matrix with rows from a source matrix based on the specified importer.
  /*!
    \param In
           SourceMatrix - Matrix from which values are imported into the "\e this" matrix.
    \param In
           Importer - A Petra_Import object specifying the communication required.

    \param In
           CombineMode - A Petra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
    int Import(const Petra_RDP_VBR_Matrix& SourceMatrix, const Petra_Import & Importer, 
	       Petra_CombineMode CombineMode);

    //! Fills a matrix with rows from a source matrix based on the specified exporter.
  /*!
    \param In
           SourceMatrix - Matrix from which values are imported into the "\e this" matrix.
    \param In
           Exporter - A Petra_Export object specifying the communication required.

    \param In
           CombineMode - A Petra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
    int Export(const Petra_RDP_VBR_Matrix& SourceMatrix, 
	       const Petra_Export & Exporter, Petra_CombineMode CombineMode);

    //! Fills a matrix with rows from a source matrix based on the specified exporter.
  /*!
    \param In
           SourceMatrix - Matrix from which values are imported into the "\e this" matrix.
    \param In
           Exporter - A Petra_Export object specifying the communication required.  Communication is done
	   in reverse of an export.

    \param In
           CombineMode - A Petra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
    int Import(const Petra_RDP_VBR_Matrix& SourceMatrix, 
	       const Petra_Export & Exporter, Petra_CombineMode CombineMode);

    //! Fills a matrix with rows from a source matrix based on the specified importer.
  /*!
    \param In
           SourceMatrix - Matrix from which values are imported into the "\e this" matrix.
    \param In
           Importer - A Petra_Import object specifying the communication required.Communication is done
	   in reverse of an import.

    \param In
           CombineMode - A Petra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
    int Export(const Petra_RDP_VBR_Matrix& SourceMatrix, 
	       const Petra_Import & Importer, Petra_CombineMode CombineMode);

    //! Returns a pointer to the Petra_Comm communicator associated with this matrix.
    const Petra_Comm & Comm() const {return(Graph_->Comm());};

 protected:
    bool Allocated() const {return(Allocated_);};
    int SetAllocated(bool Flag) {Allocated_ = Flag; return(0);};
    double *** Values() const {return(Values_);};

  // Internal utilities

  void InitializeDefaults();
  int Allocate();
  int BeginInsertValues(int BlockRow, int NumBlockEntries, 
			int * BlockIndices, bool IndicesAreLocal);
  int BeginReplaceValues(int BlockRow, int NumBlockEntries, 
			 int *BlockIndices, bool IndicesAreLocal);
  int BeginSumIntoValues(int BlockRow, int NumBlockEntries, 
			 int *BlockIndices, bool IndicesAreLocal);
  int SetupForSubmits(int BlockRow, int NumBlockEntries, int * BlockIndices, 
		      bool IndicesAreLocal, Petra_CombineMode SubmitMode);
  int EndReplaceSumIntoValues();
  int EndInsertValues();

  int CopyMat(double * A, int LDA, int NumRows, int NumCols, 
	      double * B, int LDB, bool SumInto) const;
  int BeginExtractBlockRowCopy(int BlockRow, int MaxNumBlockEntries, 
			       int & RowDim, int & NumBlockEntries, 
			       int * BlockIndices, int * ColDims, 
			       bool IndicesAreLocal) const;
  int SetupForExtracts(int BlockRow, int & RowDim, int NumBlockEntries,
		       bool ExtractView, bool IndicesAreLocal) const;
  int ExtractBlockDimsCopy(int NumBlockEntries, int * ColDims) const;
  int ExtractBlockRowPointers(int BlockRow, int MaxNumBlockEntries, 
				  int & RowDim, int & NumBlockEntries, 
				  int * BlockIndices, int * & ColDims, 
				  int * & LDAs, double ** & Values, bool IndicesAreLocal) const;
  int BeginExtractBlockRowView(int BlockRow, int & RowDim, int & NumBlockEntries, 
			       int * & BlockIndices, int * & ColDims, 
			       int * & LDAs, bool IndicesAreLocal) const;
  int ExtractBlockDimsView(int NumBlockEntries, int * & ColDims, int * & LDAs) const;
  int CopyMatDiag(double * A, int LDA, int NumRows, int NumCols, 
		  double * Diagonal) const;

  void BlockRowMultiply(bool TransA, int RowDim, int NumEntries, 
			int * BlockIndices, int RowOff,
			int * FirstElementEntryList, int * ElementSizeList,
			double Alpha, double ** As, int * LDAs, 
			double ** X, double Beta, double ** Y, int NumVectors) const;
  int InverseSums(bool DoRows, Petra_RDP_Vector& x) const;
  int Scale(bool DoRows, const Petra_RDP_Vector& x);
  void BlockRowNormInf(int RowDim, int NumEntries, 
		       int * ColDims, int * LDAs, double ** As, 
		       double * Y) const;
  void BlockRowNormOne(int RowDim, int NumEntries, int * BlockRowIndices,
		       int * ColDims, int * LDAs, double ** As, 
		       int * ColFirstElementEntryList, double * x) const;
  void SetStaticGraph(bool Flag) {StaticGraph_ = Flag;};
  int DoTransfer(const Petra_RDP_VBR_Matrix& SourceMatrix, 
		 Petra_CombineMode CombineMode,
		 int NumSameIDs, int NumPermuteIDs, int NumRemoteIDs, 
		 int NumExportIDs, 
		 int *PermuteToLIDs, int *PermuteFromLIDs, int *RemoteLIDs, 
		 int * ExportLIDs,
		 int Nsend, int Nrecv, int SizeOfPacket, int IntNsend, int IntNrecv, int IntSizeOfPacket,
		 int & LenExports, double * & Exports, int & LenIntExports, int * & IntExports,
		 int & LenImports, double * & Imports, int & LenIntImports, int * & IntImports,
#ifdef PETRA_MPI
		 GSComm_Plan & Plan, 
#endif
		 bool DoReverse);

  int CopyAndPermute(Petra_RDP_VBR_Matrix & Target, 
					   const Petra_RDP_VBR_Matrix & Source,
					   int NumSameIDs, 
					   int NumPermuteIDs, int * PermuteToLIDs,
					   int *PermuteFromLIDs);

  int Pack(const Petra_RDP_VBR_Matrix & Source, int SizeOfPacket, int IntSizeOfPacket,
				 int NumSendIDs, int * SendLIDs, double * Sends, 
				 int * IntSends);

  int UnpackAndCombine(Petra_RDP_VBR_Matrix & Target, int SizeOfPacket, int IntSizeOfPacket,
					     int NumRecvIDs, int * RecvLIDs, 
					     double * Recvs, int * IntRecvs, Petra_CombineMode CombineMode);

  bool StaticGraph() const {return(StaticGraph_);};
  Petra_CRS_Graph * Graph_;
  bool Allocated_;
  bool StaticGraph_;
  
  int NumMyBlockRows_;

  int * NumEntriesPerRow_;
  int * NumAllocatedEntriesPerRow_;
  Petra_DataAccess CV_;


  int * NumBlockEntriesPerRow_;
  int * NumAllocatedBlockEntriesPerRow_;
  int ** Indices_;
  int * ElementSizeList_;
  int * FirstElementEntryList_;

  double ***Values_;
  int ** ColDims_;
  int ** LDAs_;
  double **All_Values_;
  mutable double NormInf_;
  mutable double NormOne_;

  mutable Petra_RDP_MultiVector * ImportVector_;
  mutable Petra_RDP_MultiVector * ExportVector_;
  mutable int LenIntImports_;
  mutable int LenIntExports_;
  mutable int LenImports_;
  mutable int LenExports_;
  mutable double * Imports_;
  mutable double * Exports_;
  mutable int * IntImports_;
  mutable int * IntExports_;

  // State variables needed for constructing matrix entry-by-entry
  mutable int *TempRowDims_;
  mutable int *TempColDims_;
  mutable int *TempLDAs_;
  mutable double **TempValues_;
  mutable int LenTemps_;
  mutable int CurBlockRow_;
  mutable int CurNumBlockEntries_;
  mutable int * CurBlockIndices_;
  mutable int CurEntry_;
  mutable bool CurIndicesAreLocal_;
  mutable Petra_CombineMode CurSubmitMode_;
  
  // State variables needed for extracting entries
  mutable int CurExtractBlockRow_;
  mutable int CurExtractEntry_; 
  mutable int CurExtractNumBlockEntries_;
  mutable bool CurExtractIndicesAreLocal_;
  mutable bool CurExtractView_;
  mutable int CurRowDim_;

  // State variable for extracting block diagonal entries
  mutable int CurBlockDiag_;

};

//! << operator will work for Petra_RDP_VBR_Matrix.
ostream& operator << (ostream& os, const Petra_RDP_VBR_Matrix& A);

#endif /* _PETRA_RDP_VBR_MATRIX_H_ */
