/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_LinearProblemManager_hpp_
#define _fei_LinearProblemManager_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_mpi.h>

#include <vector>

namespace fei {
  class ParameterSet;
  class SparseRowGraph;

  /** Linear system assembly and solution manager.
  */
  class LinearProblemManager {
   public:
    //@{ \name Destructor

    /** Destructor. */
    virtual ~LinearProblemManager(){}

    //@}
    //@{ \name Initialization of size/structure

    /** Set the linear-system's global row distribution.

      @param ownedGlobalRows List of row-numbers to be owned by local processor.
    */
    virtual void setRowDistribution(const std::vector<int>& ownedGlobalRows)=0;

    /** Set the matrix-graph structure. This is the nonzero structure for
        locally-owned matrix rows.
    */
    virtual void setMatrixGraph(fei::SharedPtr<fei::SparseRowGraph> matrixGraph)=0;

    //@}
    //@{ \name Matrix access

    /** Set a specified scalar value throughout the matrix.
     */
    virtual void setMatrixValues(double scalar)=0;

    /** Query the number of local rows. This is expected to be the number of
        point-entry rows on the local processor.
    */
    virtual int getLocalNumRows()=0;

    /** Given a locally-owned global row number, query the length (number of
        nonzeros) of that row.
     */
    virtual int getRowLength(int row)=0;

    /** Given a locally-owned global row number, pass out a copy of the
        contents of that row.

        @param row Global row number

        @param len Length of user-allocated 'coefs' and 'indices' arrays.
                   if 'len' != 'getRowLength(row)', then the number of
                   coefs/indices returned will be max(len,getRowLength(row)).

        @param coefs Output array of matrix coefficients.
        @param indices Output array of column-indices.

        @return error-code 0 if successful. Non-zero return-value may indicate
        that the specified row is not locally owned.
    */
    virtual int copyOutMatrixRow(int row,
                                 int len,
                                 double* coefs,
                                 int* indices)=0;

    /** Put a C-style table (array of pointers) of coefficient data into the
        matrix.  This is a rectangular array of coefficients for
        rows/columns defined by the 'rows' and 'cols' lists.
        If the sum_into argument is true, values should be added to any that
        already exist at the specified locations. Otherwise (if sum_into is
        false) incoming values should overwrite already-existing values.
     */
    virtual int insertMatrixValues(int numRows, const int* rows,
                                   int numCols, const int* cols,
                                   const double* const* values,
                                   bool sum_into)=0;

    //@}
    //@{ \name Vector access (both soln and rhs)

    /** Set a specified scalar value throughout the vector.

        @param scalar Value to be used.

        @param soln_vector If true, scalar should be set in solution vector,
                           otherwise set rhs vector.
     */
    virtual void setVectorValues(double scalar, bool soln_vector)=0;

    /** Put coefficient data into a vector at the specified global indices.
      If any specified indices are out of range (negative or too large) the
      corresponding positions in the values array will not be referenced,
      and a positive warning code will be returned.

      @param numValues Length of caller-allocated 'globalIndices' and
                       'values' arrays.

      @param globalIndices List of global-indices specifying the locations in
                the vector for incoming values to be placed.

      @param values List of incoming values.

      @param sum_into If true, incoming values should be added to values that
                 may already be in the specified locations. If sum_into is
                 false, then incoming values should overwrite existing values.
 
      @param soln_vector If true, incoming values should be placed in the
                   solution vector. Otherwise, they should be placed in the
                   rhs vector.

      @param vectorIndex If the linear system has multiple rhs/soln vectors,
                       then this parameter specifies which vector the incoming
                       values should be put into.
    */
    virtual int insertVectorValues(int numValues,
                                   const int* globalIndices,
                                   const double* values,
                                   bool sum_into,
                                   bool soln_vector,
                                   int vectorIndex=0)=0;

    /** Copy values for the specified vector indices into the caller-allocated
      'values' array.
    */
    virtual int copyOutVectorValues(int numValues,
                                    const int* globalIndices,
                                    double* values,
                                    bool soln_vector,
                                    int vectorIndex=0) = 0;

    /** Dangerous, high-performance vector access. Return a pointer to
      local vector values. Implementations that can't support this may
      return NULL, in which case the caller will revert to using the
      copyOutVectorValues method.
    */
    virtual double* getLocalVectorValuesPtr(bool soln_vector,
                                            int vectorIndex=0) = 0;
    //@}
    //@{ \name Problem finalization

    /** Perform any necessary internal communications/synchronizations or other
        operations appropriate at the end of data input. For some
        implementations this may be a no-op.  (Trilinos/Epetra implementations
        would call 'FillComplete' on the matrix at this point.)
    */
    virtual int globalAssemble() = 0;
    //@}
  };//class LinearProblemManager

}//namespace fei

#endif // _fei_LinearProblemManager_hpp_

