#ifndef _fei_Matrix_core_hpp_
#define _fei_Matrix_core_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <fei_EqnComm.hpp>
#include <fei_fwd.hpp>
#include <fei_Vector.hpp>
#include <fei_CommUtils.hpp>
#include <fei_SSMat.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Logger.hpp>

#include <vector>

namespace fei {
class Matrix_core : protected fei::Logger {
 public:

  Matrix_core(fei::SharedPtr<fei::MatrixGraph> matrixGraph, int numLocalEqns);

  virtual ~Matrix_core();

  virtual int giveToMatrix(int numRows, const int* rows,
			   int numCols, const int* cols,
			   const double* const* values,
			   bool sumInto,
			   int format) = 0;

  /** Given a 2-D table (C-style list-of-pointers) of size numRows X numCols,
      copy the transpose of the data into the given 1-D work array and set the
      given 2-D array of pointers to point into the 1-D data. In other words,
      columns in the input will be rows in the output, and the entries in
      work_2D will point to the beginnings of those rows (former columns).

      @param numRows First dimension of 'values'. i.e., number of pointers in the
      list of pointers.

      @param numCols Second dimension of 'values'. i.e., the length of each of
      the rows (all rows have the same length).

      @param values 2-D C-style table of data. List of pointers of length numRows,
      with each pointer pointing to 'numCols' items.

      @param work_1D on exit, will have length numRows X numCols, and will contain
      the data from values. The data from 'values' will be transposed such that the
      entries from each column of 'values' will now lie contiguously.

      @param work_2D on exit, will have length numCols, and will represent a C-style
      list-of-pointers where work_2D[i] = &(work_1D[i*numRows]). In other words,
      each entry of work_2D points to a row which was a column of 'values'.
      work_2D represents a table with 'numCols' rows, each of length 'numRows'.
  */
  static void copyTransposeToWorkArrays(int numRows, int numCols,
				       const double*const* values,
				       std::vector<double>& work_1D,
				       std::vector<const double*>& work_2D);

  /** If slave-constraints have been registered with the matrix-graph, and if
      the constraints have a non-zero right-hand-side coefficient, then this
      matrix needs access to the linear-system's right-hand-side vector for
      assembling data into. For more details, see the SAND report that documents
      the FEI's constraint reduction.
  */
  void setRHS(fei::SharedPtr<fei::Vector> rhsvector);

  /** Instruct the Matrix_core to set its slave-constraint information (such
      as dependency matrix and constraint-right-hand-side vector) from the
      specified matrix-graph object.
  */
  void setSlaveInfo(fei::SharedPtr<fei::MatrixGraph> matrixGraph);

  fei::SharedPtr<fei::MatrixGraph> getMatrixGraph() { return( matrixGraph_ ); }

  std::vector<SSMat*>& getRemotelyOwnedMatrix() { return( remotelyOwned_ ); }

  void setEqnComm(fei::SharedPtr<fei::EqnComm> eqnComm);

 protected:
  void parameters(const fei::ParameterSet& paramset);

  virtual int parameters(int numParams,
			 const char* const* paramStrings);

  virtual int giveToUnderlyingMatrix(int numRows, const int* rows,
				     int numCols, const int* cols,
				     const double* const* values,
				     bool sumInto,
				     int format) = 0;

  virtual int giveToUnderlyingBlockMatrix(int row,
					  int rowDim,
					  int numCols,
					  const int* cols,
					  const int* LDAs,
					  const int* colDims,
					  const double* const* values,
					  bool sumInto) = 0;

  void setName(const char* name);

  int gatherFromOverlap(bool accumulate);

  void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph);

  /** Copy a rectangular (non-ragged) table of coefficients representing a set of
      point-entry matrix rows into a block-entry matrix row which is laid out as
      follows:
      blkValues[i], i in [0 .. numBlkCols-1], is an array containing the
      values for a single block-entry. The dimensions of the block-entry are
      numPtRows X blkColDims[i], and the block-entry values will be arranged in
      column-major order, meaning that each column of the block-entry occupies
      contiguous memory.
  */
  int copyPointRowsToBlockRow(int numPtRows,
			      int numPtCols,
			      const double*const* ptValues,
			      int numBlkCols,
			      const int* blkColDims,
			      double** blkValues);

  int convertPtToBlk(int numRows,
		       const int* rows,
		       int numCols,
		       const int* cols,
		       int* blkRows,
		       int* blkRowOffsets,
		       int* blkCols,
		       int* blkColOffsets);

  virtual int sumIntoMatrix(SSMat& mat) = 0;

  MPI_Comm getCommunicator() const { return( comm_ ); }

  fei::SharedPtr<fei::VectorSpace> vecSpace() { return( vecSpace_ ); }

  std::vector<int>& globalOffsets() { return( globalOffsets_ ); }
  const std::vector<int>& globalOffsets() const { return( globalOffsets_ ); }

  int firstLocalOffset() const { return( firstLocalOffset_ ); }
  int lastLocalOffset() const  { return( lastLocalOffset_ ); }

  int numProcs() const { return( numProcs_ ); }
  int localProc() const { return( localProc_ ); }

  bool haveBlockMatrix() const { return( haveBlockMatrix_ ); }
  void setBlockMatrix(bool flag) {haveBlockMatrix_ = flag; }

  bool haveFEMatrix() const { return( haveFEMatrix_ ); }
  void setFEMatrix(bool flag) {haveFEMatrix_ = flag; }

  int getOwnerProc(int globalEqn);

  std::string name_;

  std::vector<int> work_indices_;
  std::vector<int> work_indices2_;

  std::vector<int> work_ints_;

  std::vector<double> work_data1D_;
  std::vector<const double*> work_data2D_;

 protected:
  fei::SharedPtr<fei::EqnComm> eqnComm_;

 private:
  fei::SharedPtr<fei::Vector> rhsVector_;

  MPI_Comm comm_;

  int localProc_, numProcs_;

  fei::SharedPtr<fei::VectorSpace> vecSpace_;
  fei::SharedPtr<fei::MatrixGraph> matrixGraph_;

  std::vector<SSMat*> remotelyOwned_;

  bool haveBlockMatrix_;
  bool haveFEMatrix_;

  std::vector<int> globalOffsets_;
  int firstLocalOffset_, lastLocalOffset_;
};//class Matrix_core
}//namespace snl_fei

#endif

