/*--------------------------------------------------------------------*/
/*    Copyright 2001 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_LinSysCore_flexible_hpp_
#define _fei_LinSysCore_flexible_hpp_

#include "fei_LinearSystemCore.hpp"

/** Abstract interface that derives from LinearSystemCore and adds new
functions related to changing constraint relations.
*/

class LinSysCore_flexible : public virtual LinearSystemCore {
 public:
  virtual ~LinSysCore_flexible() {}

  /** Reset any previously-loaded lagrange multiplier constraint-relations.
   */
  virtual int resetConstraints(double s) = 0;

  /** Signal that we're done calling the setMultCREqns function.
   */
  virtual int setMultCRComplete() = 0;

  /** Supply LinSysCore_flexible with information defining the structure of
      the constraint section of the global matrix. This function is similar to the
      LinearSystemCore::setMatrixStructure function, except that only the
      constraint section of the matrix is supplied, *AND* only the "row" portion.
      i.e., the structure of C is supplied, but not C^T. Note also, that only
      the *local* rows of the C matrix are supplied.
      @param numLocalRows Number of local rows in C.
      @param globalRowNumbers Specifies whichrows in the global system matrix are
                 occupied by the constraint matrix C, but only those rows which
		 are locally stored.
      @param rowLengths List of length numLocalRows. Specifies how many
                        column-entries are in each row of the constraint matrix C.
      @param globalColIndices Table containing the column-indices in C. This
                      "table" has number-of-rows = numRows, and row i is of
                      length rowLengths[i].
      @return error-code 0 if successful
  */
  //virtual int setConstraintMatrixStructure(int numRows,
  //					   int* globalRowNumbers,
  //					   int* rowLengths,
  //					   int** globalColIndices) = 0;

  /** Signal the underlying solver library that the loading of constraints is
      now complete. (Constraint-loading functions are in LinearSystemCore.)
      Important Note: This method, LinSysCore_flexible::constraintsLoadComplete()
      will be called *INSTEAD OF* LinearSystemCore::matrixLoadComplete().
  */
  virtual int constraintsLoadComplete() = 0;
};

#endif // _fei_LinSysCore_flexible_hpp_
