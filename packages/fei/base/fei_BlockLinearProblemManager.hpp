/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_BlockLinearProblemManager_hpp_
#define _fei_BlockLinearProblemManager_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_mpi.h>

namespace fei {
  class ParameterSet;
  class SparseRowGraph;

  /** Linear system assembly and solution manager.

      This interface assumes that the parallel distribution of the
      linear system across processors is 1-D in the row-dimension,
      meaning that each processor owns a collection of complete rows of
      the matrix, and corresponding entries in the rhs and soln vectors.
  */
  class BlockLinearProblemManager {
   public:
    //@{ \name Destructor

    /** Destructor. */
    virtual ~BlockLinearProblemManager(){}

    //@}
    //@{ \name Initialization of size/structure

    /** Set the linear-system's global row distribution.

      @param ownedIDs List of IDs to be owned by local processor.

      @param dofPerOwnedID List giving num-degrees-of-freedom for each ID in
                      the 'ownedGlobalIDs' list.

      @param ghostIDs List of IDs that the local processor will need to know
        about, but which are not in the 'ownedIDs' list.  ghostIDs correspond
        to matrix columns for which the corresponding row is owned by a remote
        processor.

      @param dofPerGhostID List giving num-degrees-of-freedom for each ID in
                      the 'ownedGlobalIDs' list.
    */
    virtual void setRowDistribution(const std::vector<int>& ownedIDs,
                                    const std::vector<int>& dofPerOwnedID,
                                    const std::vector<int>& ghostIDs,
                                    const std::vector<int>& dofPerGhostID)=0;

    /** Set the matrix-graph structure. This is the nonzero structure for
        locally-owned matrix rows.
    */
    virtual void setMatrixGraph(fei::SharedPtr<fei::SparseRowGraph> matrixGraph)=0;

    //@}
    //@{ \name Matrix access

    /** Set a specified scalar value throughout the matrix.
     */
    virtual void setMatrixValues(double scalar)=0;

    /** Query the number of locally-owned matrix block-rows.
    */
    virtual int getNumOwnedIDs()=0;

    /** Given a locally-owned matrix row, query the point-entry-length
       (number of scalar nonzeros) of that row.  Note that although the
       block-row 'ownedID' corresponds to multiple point-rows, each of
       those point-rows has the same length.
     */
    virtual int getRowPointLength(int ownedID)=0;

    /** Given a locally-owned matrix row, query the block-entry-length
       (number of block columns) of that row.
     */
    virtual int getRowBlockLength(int ownedID)=0;

    /** Given a locally-owned point-entry row, pass out a copy of the
        contents of that row.

        @param ownedID Global row id.

        @param dofOffset Offset into block-row 'ownedID'.

        @param numColIDs Length of caller-allocated colIDs and dofPerColID
                         arrays.

        @param numCoefs Length of caller-allocated coefs array.

        @param colIDs Output array of block-column-indices.

        @param dofPerColID Output array which will hold the number of
                           scalar dof per column-id. i.e., the number of
                           scalar coefficients that correspond to colIDs[i]
                           is given by dofPerID[i].

        @param coefs Output array of matrix coefficients.

        @return error-code 0 if successful. Non-zero return-value may indicate
                that the specified row is not locally owned, or that the
                lengths of the user-allocated arrays aren't consistent.
    */
    virtual int copyOutMatrixRow(int ownedID, int dofOffset,
                                 int numColIDs, int numCoefs,
                                 int* colIDs,
                                 int* dofPerColID,
                                 double* coefs);

    /** Put a C-style table (array of pointers) of coefficient data into the
        matrix.  This is a rectangular array of coefficients for a single
        block-entry specified by rowID and colID.
        If the sum_into argument is true, values should be added to any that
        already exist at the specified locations. Otherwise (if sum_into is
        false) incoming values should overwrite already-existing values.
     */
    virtual int insertMatrixValues(int rowID, int numRowDof,
                                   int colID, int numColDof,
                                   const double* const* values,
                                   bool sum_into)=0;

    /** Put a single coefficient into the matrix.
        If the sum_into argument is true, value should be added to any that
        already exists at the specified location. Otherwise (if sum_into is
        false) incoming value should overwrite the already-existing value.
     */
    virtual int insertMatrixValues(int rowID, int rowDofOffset,
                                   int colID, int colDofOffset,
                                   double value,
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
      If any specified indices are out of range (negative or too large), a
      positive warning code will be returned and the corresponding positions
      in the values array will not be referenced.

      @param ID Global id of block-entry position at which incoming values
                are to be stored.

      @param numDof Number of scalar coefficients that correspond to ID.

      @param values Input list of coefficients.

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
    virtual int insertVectorValues(int ID,
                                   int numDof,
                                   const double* values,
                                   bool sum_into,
                                   bool soln_vector,
                                   int vectorIndex=0)=0;

    /** Copy values for the specified vector id into the caller's
      'values' array.
    */
    virtual int copyOutVectorValues(int ID,
                                    int numDof,
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
    //@{ \name Solution

    /** Solve the linear-system. (Most implementations will first internally
      perform the globalAssemble() operation if the caller hasn't called it
      explicitly.)

      @return Error/status code. At the implementation's discretion, the
        return-code should be a value indicating the success/error status
        of the solve. Generally a return-value of 0 will indicate a successful
        solve.
    */
    virtual int solve(const fei::ParameterSet& parameters) = 0;

    //@}
  };//class BlockLinearProblemManager

}//namespace fei

#endif // _fei_BlockLinearProblemManager_hpp_

