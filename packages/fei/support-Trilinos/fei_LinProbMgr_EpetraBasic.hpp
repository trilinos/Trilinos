/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_LinProbMgr_EpetraBasic_hpp_
#define _fei_LinProbMgr_EpetraBasic_hpp_

#include <fei_LinearProblemManager.hpp>

#include <fei_Include_Trilinos.hpp>
#include <fei_SharedPtr.hpp>

#include <vector>

class LinProbMgr_EpetraBasic : public fei::LinearProblemManager {
 public:
  LinProbMgr_EpetraBasic(MPI_Comm comm);
  virtual ~LinProbMgr_EpetraBasic();

  /** Set the linear-system's global row distribution.

    @param ownedGlobalRows List of row-numbers to be owned by local processor.
  */
  void setRowDistribution(const std::vector<int>& ownedGlobalRows);

  /** Set the matrix-graph structure. This is the nonzero structure for
      locally-owned matrix rows.
  */
  void setMatrixGraph(fei::SharedPtr<fei::SparseRowGraph> matrixGraph);

  /** Set a specified scalar value throughout the matrix.
   */
  void setMatrixValues(double scalar);

  /** Set a specified scalar value throughout the vector.

      @param scalar Value to be used.

      @param soln_vector If true, scalar should be set in solution vector,
                         otherwise set rhs vector.
  */
  void setVectorValues(double scalar, bool soln_vector);

  /** Query the number of local rows. This is expected to be the number of
      point-entry rows on the local processor.
  */
  int getLocalNumRows();

  /** Given a locally-owned global row number, query the length (number of
      nonzeros) of that row.
   */
  int getRowLength(int row);

  /** Given a locally-owned global row number, pass out a copy of the
      contents of that row.
      @param row Global row number
      @param len Length of the user-allocated arrays coefs and indices.
      @param coefs User-allocated array which will hold matrix coefficients
      on output.
      @param indices User-allocated array which will hold column-indices on
      output.

      @return error-code 0 if successful. Non-zero return-value may indicate
      that the specified row is not locally owned.
  */
  int copyOutMatrixRow(int row, int len,
                       double* coefs, int* indices);

  /** Put a C-style table (array of pointers) of coefficient data into the
      matrix.  This is a rectangular array of coefficients for
      rows/columns defined by the 'rows' and 'cols' lists.
      If the sum_into argument is true, values should be added to any that
      already exist at the specified locations. Otherwise (if sum_into is
      false) incoming values should overwrite already-existing values.
   */
  int insertMatrixValues(int numRows, const int* rows,
                         int numCols, const int* cols,
                         const double* const* values,
                         bool sum_into);

  /** Put coefficient data into a vector at the specified global indices.
    If any specified indices are out of range (negative or too large), the
     corresponding positions in the values array will not be referenced and
     a positive warning code will be returned.

    @param numValues Number of coefficient values being input.

    @param globalIndices List of global-indices specifying the locations in
              the vector for incoming values to be placed.

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
  int insertVectorValues(int numValues,
                         const int* globalIndices,
                         const double* values,
                         bool sum_into,
                         bool soln_vector,
                         int vectorIndex=0);

  /** Copy values for the specified vector indices into the caller-allocated
    'values' array.
  */
  int copyOutVectorValues(int numValues,
                           const int* globalIndices,
                           double* values,
                           bool soln_vector,
                           int vectorIndex=0);

  /** Dangerous, high-performance vector access. Return a pointer to
    local vector values. Implementations that can't support this may
    return NULL, in which case the caller will revert to using the
    copyOutVectorValues method.
  */
  double* getLocalVectorValuesPtr(bool soln_vector,
                                  int vectorIndex=0);

  /** Perform any necessary internal communications/synchronizations or other
      operations appropriate at the end of data input. For some
      implementations this may be a no-op.  (Trilinos/Epetra implementations
      would call 'FillComplete' on the matrix at this point.)
  */
  int globalAssemble();

  /** Return the Epetra matrix.
  */
  fei::SharedPtr<Epetra_CrsMatrix> get_A_matrix();

  /** Return the rhs Epetra vector.
  */
  fei::SharedPtr<Epetra_MultiVector> get_rhs_vector();

  /** Return the solution Epetra vector.
  */
  fei::SharedPtr<Epetra_MultiVector> get_solution_vector();

 private:
  MPI_Comm comm_;
  std::vector<int> ownedRows_;
  fei::SharedPtr<Epetra_Comm> epetra_comm_;
  fei::SharedPtr<Epetra_Map> epetra_rowmap_;
  fei::SharedPtr<fei::SparseRowGraph> fei_srgraph_;
  fei::SharedPtr<Epetra_CrsGraph> crsgraph_;
  fei::SharedPtr<Epetra_CrsMatrix> A_;
  int numVectors_;
  fei::SharedPtr<Epetra_MultiVector> x_;
  fei::SharedPtr<Epetra_MultiVector> b_;
};

#endif // _LinProbMgr_EpetraBasic_hpp_

