#ifndef _SSMat_hpp_
#define _SSMat_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iosfwd.hpp>
#include <feiArray.hpp>
#include <fei_SSVec.hpp>

/** SSMat stands for Super-Sparse Matrix. It is a data structure that can
hold and perform operations on matrices that are arbitrarily sparse, i.e.,
whose rows/columns don't necessarily form a contiguous set, don't necessarily
start at 0 or 1, and aren't even necessarily sorted. E.g., an SSMat instance
may be a matrix having rows 94, 38, and 1123, with each of those rows 
containing an equally arbitrary set of column-indices. An SSMat may contain an
element-contribution (a finite-element stiffness array), whose row/column
numbers are the scatter indices into the global system matrix being assembled.

What is the point of this? I need to do matrix operations on partial (e.g.,
element-wise) contributions to a global system matrix during the assembly 
process.

SSMat provides methods for forming the matrix-matrix product of two SSMat
objects, the matrix-vector product of a SSMat and a SSVec object, etc.

Since the contents of these matrices can be so arbitrary, dimension-checking
prior to forming a matrix-matrix product is meaningless. Mathematically
speaking, each matrix is assumed to extend as far in either direction (row-
space or column-space) as would be necessary for the product to be performed,
and coefficients that don't explicitly appear, are assumed to be zero.

An unfortunate by-product of the complete generality allowed for above, is that
the matrix operations are going to be very slow. (But the size of an element-
stiffness will generally be small, comparatively speaking...)

SSMat provides a couple of constructors for "wrapping" an SSMat object
around existing data. These constructors are intended to be as light-weight
as possible, so the SSMat object keeps pointers to the existing data. This
means there is an inherent danger that the data may be destroyed before the
SSMat, leaving the SSMat holding bad pointers. USER BEWARE.
*/

class SSMat {
 public:
  /** Default constructor. Creates an empty SSMat object.
   */
  SSMat(int alloc_increment=64, int row_alloc_increment=64);

  /** Constructor to create an SSMat from the data in an existing EqnBuffer
      object. 
  */
  SSMat(EqnBuffer& eqnBuf);

  /** Copy constructor */
  SSMat(const SSMat& src);

  /** Standard destructor. */
  virtual ~SSMat();

  /** Function to clear the values but not delete any memory. Only works if
      this SSMat object was created with the default constructor.
      @return error-code 0 if successful, -1 if this SSMat object was created
      using one of the non-default constructors.
  */
  void logicalClear();

  /** Matrix-matrix product between 'this' and 'inMat' with result stored
      in 'result'.
  */
  int matMat(SSMat& inMat, SSMat& result, bool storeResultZeros=true);

  /** Matrix-transpose-matrix product between transpose('this') and 'inMat'
      with result stored in 'result'.
  */
  int matTransMat(SSMat& inMat, SSMat& result, bool storeResultZeros=true);

  /** Matrix-vector product between 'this' and 'inVec' with result stored in
      'result'.
  */
  int matVec(SSVec& inVec, SSVec& result);

  /** Matrix-transpose-vector product between transpose('this') and 'inVec'
      with result stored in 'result'.
  */
  int matTransVec(SSVec& inVec, SSVec& result);

  /** An in-efficient but convenient read-write coefficient access function.
      @param row Global row number.
      @param col Global column number.
      @param coefPtr Output. A pointer to the coefficient.
      @return error-code 0 if position (row,col) exists, -1 if it doesn't.
   */
  int coefficientPointer(int row, int col, double*& coefPtr);

  /** Accumulate a coefficient into the matrix. (Only if this matrix was
      constructed using the default constructor.) Creates a new position if it
      didn't already exist.
      @return error-code 0 if successful
  */
  int sumInCoef(int row, int col, double coef);

  /**Put a coefficient into the matrix. (Only if this matrix was
      constructed using the default constructor.) Creates a new position if it
      didn't already exist.
      @return error-code 0 if successful
  */
  int putCoef(int row, int col, double coef);

  /** Accumulate a row of coefs into the matrix. (Only if this matrix was
      constructed using the default constructor.) Creates any new positions that
      don't already exist.
  */
  int sumInRow(int row, const int* cols, const double* coefs, int len);

  /** Put a row of coefs into the matrix, replacing any that were already
      present. (Only if this matrix was
      constructed using the default constructor.) Creates any new positions that
      don't already exist.
  */
  int putRow(int row, const int* cols, const double* coefs, int len);

  /** Add position (row,col) for a coefficient. This function only works
      if this SSMat object was constructed using the default constructor.
      @param row Global row number.
      @param col Global column number.
  */
  void createPosition(int row, int col);

  /** Return the row-numbers contained in this SSMat. */
  feiArray<int>& getRowNumbers() { return( *rowNumbers_ ); }

  /** Return the row-numbers contained in this SSMat. */
  const feiArray<int>& getRowNumbers() const { return( *rowNumbers_ ); }

  /** Return the rows for this SSMat. */
  feiArray<SSVec*>& getRows() { return( *rows_ ); }

  /** Return the rows for this SSMat. */
  const feiArray<SSVec*>& getRows() const { return( *rows_ ); }

  /** Return a specified row, or NULL if the specified row doesn't exist. If
      the default argument create_if_necessary is true, return a newly created
      row if it doesn't already exist. create_if_necessary defaults to false. */
  SSVec* getRow(int row, bool create_if_necessary=false);

  SSVec* rowContainingCol(int col, int& offsetInRow);

  bool structurallySymmetric;

  int numNonzeros();

  bool operator==(const SSMat& lhs);

  bool operator!=(const SSMat& lhs);

  int getMinCol();

  int getMaxCol();

  void writeToStream(FEI_OSTREAM& os);

  SSMat& operator=(const SSMat& src);
  SSMat& operator+=(const SSMat& src);

 private:

  /** Create position (row,col) and return the indices into the coefs_ table
      at which the position is located. If the position already exists, simply
      return those indices.
  */
  void createPosition(int row, int col, int& rowIndex, int& colIndex);

  void appendRow(int row);
  SSVec* insertRow(int row, int index);
  int whichConstructor_;

  feiArray<int>* rowNumbers_;
  feiArray<SSVec*>* rows_;

  int highWaterMark_;
  int row_alloc_incr_;
};

#ifndef _fei_ostream_ops_hpp_
#include <fei_ostream_ops.hpp>
#endif

#endif
