#ifndef _SSGraph_hpp_
#define _SSGraph_hpp_
/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/


#include <fei_macros.hpp>
#include <fei_fwd.hpp>
#include <feiArray.hpp>

enum { SS_Constr_Default, SS_Constr_EqnBuf,
       SS_Constr_RawArrays, SS_Constr_RawArrays2, SS_Constr_RawArraysSymm };


/** SSGraph stands for Super-Sparse Graph. It is a data structure that can
hold the graph or structure of matrices that are arbitrarily sparse, i.e.,
whose rows/columns don't necessarily form a contiguous set, don't necessarily
start at 0 or 1, and aren't even necessarily sorted. E.g., an SSGraph instance
may be a matrix having rows 94, 38, and 1123, with each of those rows 
containing an equally arbitrary set of column-indices. An SSGraph may contain an
element-contribution (a finite-element stiffness array), whose row/column
numbers are the scatter indices into the global system matrix being assembled.

What is the point of this? I need to do matrix operations on partial (e.g.,
element-wise) contributions to a global system matrix during the assembly 
process.

SSGraph provides a couple of constructors for "wrapping" an SSGraph object
around existing data. These constructors are intended to be as light-weight
as possible, so the SSGraph object keeps pointers to the existing data. This
means there is an inherent danger that the data may be destroyed before the
SSGraph, leaving the SSGraph holding bad pointers. USER BEWARE.
*/

class SSGraph {
 public:
  /** Default constructor. Creates an empty SSGraph object.
   */
  SSGraph();

  /** Constructor to create an SSGraph from raw arrays of indices.
  */
  SSGraph(int numRows, const int* rowNumbers,
	int numCols, const int* colIndices);

  /** Constructor to create an SSGraph from raw arrays in a different
      format...
      @param numRows
      @param rowNumbers
      @param numColsPerRow Number of column-indices for each row.
      @param rowColOffsets i-th entry provides the offset into the colIndices
      list at which the column-indices for the i-th row may be found. This is
      redundant, strictly speaking, since each row must be of length 
      numColsPerRow.
      @param colIndices Packed list of column-indices, length is 
      numRows*numColsPerRow
  */
  SSGraph(int numRows, const int* rowNumbers,
	int numColsPerRow, const int* rowColOffsets,
	const int* colIndices);

  /** Standard destructor. */
  ~SSGraph();

  /** Function to clear the structure but not delete any memory. Only works if
      this SSGraph object was created with the default constructor.
  */
  void logicalClear();

  /** Add position (row,col) for a coefficient. This function only works
      if this SSGraph object was constructed using the default constructor.
      @param row Global row number.
      @param col Global column number.
  */
  void createPosition(int row, int col);

  /** Return the row-numbers contained in this SSGraph. */
  feiArray<int>& getRows() { return( *rows_ ); }

  /** Return the row-lengths associated with this SSGraph. */
  feiArray<int>& getRowLengths() { return( *rowLengths_ ); }

  /** Return the column-indices for this SSGraph. */
  feiArray<feiArray<int>*>& getIndices() { return( *indices_ ); }

 private:
  /** Create position (row,col) and return the indices into the coefs_ table
      at which the position is located. If the position already exists, simply
      return those indices. If this SSGraph instance wasn't created using the
      default constructor, then rowIndex and colIndex aren't referenced.
  */
  void createPosition(int row, int col, int& rowIndex, int& colIndex);

  void appendRow(int row);
  void insertRow(int row, int index);
  int whichConstructor_;

  feiArray<int>* rows_;
  feiArray<int>* rowLengths_;
  feiArray<feiArray<int>*>* indices_;
};

#endif // _SSGraph_hpp_
