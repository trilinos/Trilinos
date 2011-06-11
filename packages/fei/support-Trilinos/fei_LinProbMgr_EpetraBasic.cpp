/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_trilinos_macros.hpp>
#include <fei_SparseRowGraph.hpp>

#ifdef HAVE_FEI_EPETRA

#include <fei_LinProbMgr_EpetraBasic.hpp>

LinProbMgr_EpetraBasic::LinProbMgr_EpetraBasic(MPI_Comm comm)
 : comm_(comm),
   ownedRows_(),
   epetra_comm_(),
   epetra_rowmap_(),
   fei_srgraph_(),
   crsgraph_(),
   A_(),
   numVectors_(1),
   x_(),
   b_()
{
#ifdef FEI_SER
  fei::SharedPtr<Epetra_Comm> ecomm(new Epetra_SerialComm);
#else
  fei::SharedPtr<Epetra_Comm> ecomm(new Epetra_MpiComm(comm));
#endif

  epetra_comm_ = ecomm;
}

LinProbMgr_EpetraBasic::~LinProbMgr_EpetraBasic()
{
}


void LinProbMgr_EpetraBasic
::setRowDistribution(const std::vector<int>& ownedGlobalRows)
{
  if (ownedGlobalRows == ownedRows_) {
    return;
  }

  if (!ownedRows_.empty()) {
    throw std::runtime_error("setRowDistribution called multiple times with different distributions. not allowed.");
  }

  int* rows = const_cast<int*>(&ownedGlobalRows[0]);
  epetra_rowmap_.reset(new Epetra_Map(-1, ownedGlobalRows.size(),
                                      rows, 0, //indexBase
                                      *epetra_comm_));

  x_.reset(new Epetra_MultiVector(*epetra_rowmap_, numVectors_));

  b_.reset(new Epetra_MultiVector(*epetra_rowmap_, numVectors_));
}

void LinProbMgr_EpetraBasic
::setMatrixGraph(fei::SharedPtr<fei::SparseRowGraph> matrixGraph)
{
  if (fei_srgraph_.get() != NULL) {
    if (*fei_srgraph_ != *matrixGraph) {
      throw std::runtime_error("setMatrixGraph called multiple times with different graphs. not allowed.");
    }
    return;
  }
  else {
    fei_srgraph_ = matrixGraph;
  }

  if (epetra_rowmap_.get() == NULL) {
    setRowDistribution(matrixGraph->rowNumbers);
  }

  if ((int)fei_srgraph_->rowNumbers.size() != epetra_rowmap_->NumMyElements()) {
    throw std::runtime_error("setMatrixGraph: num-rows not consistent with value from setRowDistribution");
  }

  //We'll create and populate a Epetra_CrsGraph object.

  std::vector<int>& rowNumbers = fei_srgraph_->rowNumbers;
  std::vector<int>& rowOffsets = fei_srgraph_->rowOffsets;

  //First create an array of num-indices-per-row.
  std::vector<int> numIndicesPerRow; numIndicesPerRow.reserve(rowNumbers.size());

  int err;
  unsigned i;
  for(i=0; i<numIndicesPerRow.size(); ++i) {
    numIndicesPerRow.push_back(rowOffsets[i+1] - rowOffsets[i]);
  }

  bool static_profile = true;

  crsgraph_.reset(new Epetra_CrsGraph(Copy, *epetra_rowmap_,
                  &numIndicesPerRow[0], static_profile));

  //Now put in all the column-indices

  std::vector<int>& colIndices = fei_srgraph_->packedColumnIndices;
  for(i=0; i<rowNumbers.size(); ++i) {
    int offset = rowOffsets[i];
    err = crsgraph_->InsertGlobalIndices(rowNumbers[i], numIndicesPerRow[i],
                                        &colIndices[offset]);
    if (err != 0) {
      throw std::runtime_error("setMatrixGraph: err from Epetra_CrsGraph::InsertGlobalIndices.");
    }
  }

  err = crsgraph_->FillComplete();
  if (err != 0) {
    throw std::runtime_error("setMatrixGraph: err from Epetra_CrsGraph::FillComplete.");
  }
 
  //and finally, create a matrix.
  A_.reset(new Epetra_CrsMatrix(Copy, *crsgraph_));
}

void LinProbMgr_EpetraBasic::setMatrixValues(double scalar)
{
  int err = A_->PutScalar(scalar);
  if (err != 0) {
    throw std::runtime_error("error in Epetra_CrsMatrix->PutScalar");
  }
}

void LinProbMgr_EpetraBasic::setVectorValues(double scalar,
                                             bool soln_vector)
{
  int err = soln_vector ?
    x_->PutScalar(scalar) : b_->PutScalar(scalar);
  if (err != 0) {
    throw std::runtime_error("error in Epetra_MultiVector->PutScalar");
  }
}

int LinProbMgr_EpetraBasic::getLocalNumRows()
{
  if (epetra_rowmap_.get() == NULL) return(-1);

  return(epetra_rowmap_->NumMyElements());
}

int LinProbMgr_EpetraBasic::getRowLength(int row)
{
  if (A_.get() == NULL) return(-1);

  return( A_->NumGlobalEntries(row) );
}

int LinProbMgr_EpetraBasic::copyOutMatrixRow(int row, int len,
                                             double* coefs, int* indices)
{
  int dummy;
  return( A_->ExtractGlobalRowCopy(row, len, dummy, coefs, indices) );
}

int LinProbMgr_EpetraBasic::insertMatrixValues(int numRows, const int* rows,
                                               int numCols, const int* cols,
                                               const double* const* values,
                                               bool sum_into)
{
  int* nc_cols = const_cast<int*>(cols);
  double** nc_values = const_cast<double**>(values);
  int err = 0;
  if (sum_into) {
    for(int i=0; i<numRows; ++i) {
      err = A_->SumIntoGlobalValues(rows[i], numCols, nc_values[i], nc_cols);
      if (err < 0) {
        return(err);
      }
    }
  }
  else {
    for(int i=0; i<numRows; ++i) {
      err = A_->ReplaceGlobalValues(rows[i], numCols, nc_values[i], nc_cols);
      if (err < 0) {
        return(err);
      }
    }
  }
  return(err);
}

int LinProbMgr_EpetraBasic::insertVectorValues(int numValues,
                                               const int* globalIndices,
                                               const double* values,
                                               bool sum_into,
                                               bool soln_vector,
                                               int vectorIndex)
{
  double* localvaluesptr = soln_vector ?
    x_->Pointers()[vectorIndex] : b_->Pointers()[vectorIndex];

  int min_my_gid = epetra_rowmap_->MinMyGID();
  int returnValue = 0;

  if (sum_into) {
    for(int i=0; i<numValues; ++i) {
      int offset = globalIndices[i] - min_my_gid;
      if (offset < 0) {
        returnValue = 1;
        continue;
      }
      localvaluesptr[offset] += values[i];
    }
  }
  else {
    for(int i=0; i<numValues; ++i) {
      int offset = globalIndices[i] - min_my_gid;
      if (offset < 0) {
        returnValue = 1;
        continue;
      }
      localvaluesptr[offset] = values[i];
    }
  }

  return(returnValue);
}

int LinProbMgr_EpetraBasic::copyOutVectorValues(int numValues,
                                                const int* globalIndices,
                                                double* values,
                                                bool soln_vector,
                                                int vectorIndex)
{
  double* localvaluesptr = soln_vector ?
    x_->Pointers()[vectorIndex] : b_->Pointers()[vectorIndex];

  int min_my_gid = epetra_rowmap_->MinMyGID();

  for(int i=0; i<numValues; ++i) {
    int offset = globalIndices[i] - min_my_gid;
    values[i] = localvaluesptr[offset];
  }
  return(0);
}

double* LinProbMgr_EpetraBasic::getLocalVectorValuesPtr(bool soln_vector,
                                                        int vectorIndex)
{
  double* localvaluesptr = soln_vector ?
    x_->Pointers()[vectorIndex] : b_->Pointers()[vectorIndex];

  return(localvaluesptr);
}

int LinProbMgr_EpetraBasic::globalAssemble()
{
  if (!A_->Filled()) {
    //assumes the matrix is square!
    int err = A_->FillComplete();
    if (err != 0) {
      return(err);
    }
  }

  if (!A_->StorageOptimized()) {
    A_->OptimizeStorage();
  }

  return(0);
}

fei::SharedPtr<Epetra_CrsMatrix>
LinProbMgr_EpetraBasic::get_A_matrix()
{
  return( A_ );
}

fei::SharedPtr<Epetra_MultiVector>
LinProbMgr_EpetraBasic::get_rhs_vector()
{
  return( b_ );
}

fei::SharedPtr<Epetra_MultiVector>
LinProbMgr_EpetraBasic::get_solution_vector()
{
  return( x_ );
}

//HAVE_FEI_EPETRA
#endif

