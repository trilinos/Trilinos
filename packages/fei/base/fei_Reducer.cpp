
/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_Reducer.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix.hpp>
#include <fei_Matrix_core.hpp>
#include <fei_Vector.hpp>
#include <fei_Graph_Impl.hpp>
#include <snl_fei_ArrayUtils.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_SparseRowGraph.hpp>
#include <fei_Vector.hpp>

namespace fei {
Reducer::Reducer(fei::SharedPtr<SSMat> globalSlaveDependencyMatrix,
                 fei::SharedPtr<SSVec> g_vector,
                 MPI_Comm comm)
 : D_(globalSlaveDependencyMatrix),
   slavesPtr_(NULL),
   Kii_(),
   Kid_(),
   Kdi_(),
   Kdd_(),
   tmpMat1_(),
   tmpMat2_(),
   fi_(),
   fd_(),
   xi_(),
   xd_(),
   tmpVec1_(),
   tmpVec2_(),
   g_(g_vector),
   g_nonzero_(false),
   localUnreducedEqns_(),
   localReducedEqns_(),
   nonslaves_(),
   reverse_(),
   isSlaveEqn_(NULL),
   numGlobalSlaves_(0),
   numLocalSlaves_(0),
   firstLocalReducedEqn_(0),
   lastLocalReducedEqn_(0),
   lowestGlobalSlaveEqn_(0),
   highestGlobalSlaveEqn_(0),
   localProc_(0),
   numProcs_(1),
   comm_(comm),
   dbgprefix_("Reducer: "),
   mat_counter_(0),
   soln_vec_counter_(0),
   rhs_vec_counter_(0),
   bool_array_(0),
   int_array_(0),
   double_array_(0),
   array_len_(0),
   work_1D_(),
   work_2D_()
{
  initialize();
}

void
Reducer::initialize()
{
  numGlobalSlaves_ = D_->getRowNumbers().length();
  slavesPtr_ = D_->getRowNumbers().dataPtr();
  lowestGlobalSlaveEqn_ = slavesPtr_[0];
  highestGlobalSlaveEqn_ = slavesPtr_[numGlobalSlaves_-1];

  if (g_.get() != NULL) {
    if (g_->length() > 0) {
      double* gptr = g_->coefs().dataPtr();
      for(int i=0; i<g_->length(); ++i) {
        if (gptr[i] != 0.0) {
          g_nonzero_ = true;
          break;
        }
      }
    }
  }

  //The following code sets up the mappings (nonslaves_ and reverse_)
  //necessary to implement the method 'translateFromReducedEqn'.
  //This code is ugly and inelegant, and it took me quite a while to get
  //it right. I won't even bother trying to comment it, but note that its
  //correctness is tested in test_utils/test_Reducer.cpp.

  nonslaves_.resize(numGlobalSlaves_*2);

  int nonslave_counter = 0;
  int slaveOffset = 0;
  int num_consecutive = 0;
  while(slaveOffset < numGlobalSlaves_-1) {
    int gap = slavesPtr_[slaveOffset+1]-slavesPtr_[slaveOffset]-1;
    if (gap > 0) {
      nonslaves_[nonslave_counter] = slavesPtr_[slaveOffset]+1;
      nonslaves_[numGlobalSlaves_+nonslave_counter++] = num_consecutive;
      num_consecutive = 0;
    }
    else {
      ++num_consecutive;
    }
    ++slaveOffset;
  }

  nonslaves_[nonslave_counter] = highestGlobalSlaveEqn_+1;
  nonslaves_[numGlobalSlaves_+nonslave_counter++] = num_consecutive;

  reverse_.resize(nonslave_counter);
  int first = lowestGlobalSlaveEqn_;
  reverse_[0] = first;
  for(int i=1; i<nonslave_counter; ++i) {
    reverse_[i] = reverse_[i-1] +
         (nonslaves_[i]-nonslaves_[i-1] - nonslaves_[numGlobalSlaves_+i] - 1);
  }

#ifndef FEI_SER
    MPI_Comm_rank(comm_, &localProc_);
    MPI_Comm_size(comm_, &numProcs_);
#endif

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<< "ctor, numGlobalSlaves="<<numGlobalSlaves_
      << FEI_ENDL;
  }

  if (numGlobalSlaves_ < 1) {
    throw std::runtime_error("ERROR: don't use fei::Reducer when numGlobalSlaves==0. Report to Alan Williams.");
  }
}

Reducer::Reducer(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
 : D_(matrixGraph->getSlaveDependencyMatrix()),
   slavesPtr_(NULL),
   Kii_(),
   Kid_(),
   Kdi_(),
   Kdd_(),
   tmpMat1_(),
   tmpMat2_(),
   fi_(),
   fd_(),
   xi_(),
   xd_(),
   tmpVec1_(),
   tmpVec2_(),
   g_(),
   g_nonzero_(false),
   localUnreducedEqns_(),
   localReducedEqns_(),
   nonslaves_(),
   reverse_(),
   isSlaveEqn_(NULL),
   numGlobalSlaves_(0),
   numLocalSlaves_(0),
   firstLocalReducedEqn_(0),
   lastLocalReducedEqn_(0),
   lowestGlobalSlaveEqn_(0),
   highestGlobalSlaveEqn_(0),
   localProc_(0),
   numProcs_(1),
   comm_(),
   dbgprefix_("Reducer: "),
   mat_counter_(0),
   soln_vec_counter_(0),
   rhs_vec_counter_(0),
   bool_array_(0),
   int_array_(0),
   double_array_(0),
   array_len_(0)
{
  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph->getRowSpace();
  comm_ = vecSpace->getCommunicator();
  initialize();

  int num = vecSpace->getNumIndices_Owned();
  std::vector<int> indices(num);
  vecSpace->getIndices_Owned(num, &indices[0], num);
  setLocalUnreducedEqns(indices);
}

Reducer::~Reducer()
{
  delete [] isSlaveEqn_;   isSlaveEqn_ = 0;
  delete [] bool_array_;   bool_array_ = 0;
  delete [] int_array_;    int_array_ = 0;
  delete [] double_array_; double_array_ = 0;
  array_len_ = 0;
}

void
Reducer::setLocalUnreducedEqns(const std::vector<int>& localUnreducedEqns)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<< "setLocalUnreducedEqns, numLocalEqns="
      <<localUnreducedEqns.size() << FEI_ENDL;
  }

  if (localUnreducedEqns_ == localUnreducedEqns) {
    return;
  }

  localUnreducedEqns_ = localUnreducedEqns;

  int num = localUnreducedEqns_.size();

  isSlaveEqn_ = num > 0 ? new bool[localUnreducedEqns_.size()] : NULL;

  numLocalSlaves_ = 0;

  for(unsigned i=0; i<localUnreducedEqns_.size(); ++i) {
    int idx = snl_fei::binarySearch(localUnreducedEqns_[i],
                                    slavesPtr_, numGlobalSlaves_);
    if (idx < 0) {
      isSlaveEqn_[i] = false;
    }
    else {
      isSlaveEqn_[i] = true;
      ++numLocalSlaves_;

      if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
        FEI_OSTREAM& os = *output_stream_;
        os << dbgprefix_<<" slave " << localUnreducedEqns_[i] << " depends on ";
        feiArray<int>& indices = D_->getRows()[idx]->indices();
        for(int j=0; j<indices.length(); ++j) {
          os << indices[j] << " ";
        }
        os << FEI_ENDL;
      }
    }

  }

  int num_slaves_on_lower_procs = 0;


#ifndef FEI_SER

  if (numProcs_ > 1) {
    std::vector<int> procNumLocalSlaves(numProcs_);

    MPI_Allgather(&numLocalSlaves_, 1, MPI_INT, &procNumLocalSlaves[0],
                  1, MPI_INT, comm_);

    for(int p=0; p<localProc_; ++p) {
      num_slaves_on_lower_procs += procNumLocalSlaves[p];
    }
  }
#endif

  unsigned first_non_slave_offset = 0;
  while(first_non_slave_offset < localUnreducedEqns_.size() &&
        isSlaveEqn_[first_non_slave_offset] == true) {
    ++first_non_slave_offset;
  }

  firstLocalReducedEqn_ = localUnreducedEqns_[first_non_slave_offset]
      - num_slaves_on_lower_procs - first_non_slave_offset;

  int num_local_eqns = localUnreducedEqns_.size() - numLocalSlaves_;

  lastLocalReducedEqn_ = firstLocalReducedEqn_ + num_local_eqns - 1;

  localReducedEqns_.resize(num_local_eqns);

  unsigned offset = 0;
  int eqn = firstLocalReducedEqn_;
  for(unsigned i=0; i<localUnreducedEqns_.size(); ++i) {
    if (isSlaveEqn_[i] == false) {
      localReducedEqns_[offset++] = eqn++;
    }
  }

  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    if (localUnreducedEqns_.size() > 0) {
      os << dbgprefix_<< "first local eqn="
         <<localUnreducedEqns_[0]<<", last local eqn="
      <<localUnreducedEqns_[localUnreducedEqns_.size()-1] << FEI_ENDL;
    }
    os << dbgprefix_<<"numLocalSlaves_="<<numLocalSlaves_
       <<", firstLocalReducedEqn_="<<firstLocalReducedEqn_
       <<", lastLocalReducedEqn_="<<lastLocalReducedEqn_<<FEI_ENDL;
  }
}

void
Reducer::addGraphEntries(fei::SharedPtr<fei::SparseRowGraph> matrixGraph)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addGraphEntries"<<FEI_ENDL;
  }

  //iterate through the incoming matrixGraph, putting its contents into
  //Kdd, Kdi, Kid and Kii as appropriate.

  bool create_if_necessary = true;

  std::vector<int>& rowNumbers = matrixGraph->rowNumbers;
  std::vector<int>& rowOffsets = matrixGraph->rowOffsets;
  std::vector<int>& packedCols = matrixGraph->packedColumnIndices;

  for(unsigned i=0; i<rowNumbers.size(); ++i) {
    int row = rowNumbers[i];

    bool slave_row = isSlaveEqn(row);

    int rowLength = rowOffsets[i+1]-rowOffsets[i];
    int* cols = &packedCols[rowOffsets[i]];

    if (slave_row) {
      SSVec* Kdd_row = Kdd_.getRow(row, create_if_necessary);
      SSVec* Kdi_row = Kdi_.getRow(row, create_if_necessary);

      for(int j=0; j<rowLength; ++j) {
        int col = cols[j];
        bool slave_col = isSlaveEqn(col);

        if (slave_col) {
          Kdd_row->addEntry(col, 0.0);
        }
        else {
          Kdi_row->addEntry(col, 0.0);
        }
      }
    }
    else {
      //not a slave row, so add slave columns to Kid, and non-slave
      //columns to graph.
      SSVec* Kid_row = Kid_.getRow(row, create_if_necessary);
      SSVec* Kii_row = Kii_.getRow(row, create_if_necessary);

      for(int j=0; j<rowLength; ++j) {
        int col = cols[j];
        bool slave_col = isSlaveEqn(col);

        if (slave_col) {
          Kid_row->addEntry(col, 0.0);
        }
        else {
          Kii_row->addEntry(col, 0.0);
        }
      }
    }
  }
}

void
Reducer::expand_work_arrays(int size)
{
  if (size <= array_len_) return;

  array_len_ = size;
  delete [] bool_array_;
  delete [] int_array_;
  delete [] double_array_;
  bool_array_ = new bool[array_len_];
  int_array_ = new int[array_len_];
  double_array_ = new double[array_len_];
}

void
Reducer::addGraphIndices(int numRows, const int* rows,
                         int numCols, const int* cols,
                         fei::Graph& graph)
{
  bool create_if_necessary = true;

  expand_work_arrays(numCols);

  bool no_slave_cols = true;
  for(int i=0; i<numCols; ++i) {
    bool_array_[i] = isSlaveEqn(cols[i]);
    if (bool_array_[i]) no_slave_cols = false;
  }

  for(int i=0; i<numRows; ++i) {
    bool slave_row = isSlaveEqn(rows[i]);

    if (slave_row) {
      SSVec* Kdd_row = Kdd_.getRow(rows[i], create_if_necessary);
      SSVec* Kdi_row = Kdi_.getRow(rows[i], create_if_necessary);

      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          Kdd_row->addEntry(cols[j], 0.0);
        }
        else {
          Kdi_row->addEntry(cols[j], 0.0);
        }
      }
      ++mat_counter_;
    }
    else {
      //not a slave row, so add slave columns to Kid, and non-slave
      //columns to graph.
      SSVec* Kid_row = no_slave_cols ?
        NULL : Kid_.getRow(rows[i], create_if_necessary);
  
      unsigned num_non_slave_cols = 0;

      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          Kid_row->addEntry(cols[j], 0.0);
          ++mat_counter_;
        }
        else {
          int_array_[num_non_slave_cols++] = translateToReducedEqn(cols[j]);
        }
      }

      if (num_non_slave_cols > 0) {
        int reduced_row = translateToReducedEqn(rows[i]);
        graph.addIndices(reduced_row, num_non_slave_cols, int_array_);
      }
    }
  }

  if (mat_counter_ > 600) {
    assembleReducedGraph(&graph, false);
  }
}

void
Reducer::addSymmetricGraphIndices(int numIndices, const int* indices,
                                  bool diagonal,
                                  fei::Graph& graph)
{
  if (diagonal) {
    for(int i=0; i<numIndices; ++i) {
      addGraphIndices(1, &indices[i], 1, &indices[i], graph);
    }
  }
  else {
    addGraphIndices(numIndices, indices, numIndices, indices, graph);
  }
}

void
Reducer::assembleReducedGraph(fei::Graph* graph,
                              bool global_gather)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"assembleReducedGraph(fei::Graph)"<<FEI_ENDL;
    if (output_level_ >= fei::FULL_LOGS) {
      os << dbgprefix_<<"Kdi:"<<FEI_ENDL<<Kdi_;
    }
  }

  //This function performs the appropriate manipulations (matrix-matrix
  //products, etc., on the Kid,Kdi,Kdd matrices and then assembles the
  //results into the input graph structure.
  //

  //form tmpMat1_ = Kid*D
  if ( Kid_.matMat(*D_, tmpMat1_) != 0) {
    throw std::runtime_error("fei::Reducer::assembleReducedGraph ERROR 1.");
  }

  //form tmpMat2_ = D^T*Kdi
  if ( D_->matTransMat(Kdi_, tmpMat2_) != 0) {
    throw std::runtime_error("fei::Reducer::assembleReducedGraph ERROR 2.");
  }

  translateSSMatToReducedEqns(tmpMat1_);
  translateSSMatToReducedEqns(tmpMat2_);

  addSSMatToGraph(tmpMat1_, graph);
  addSSMatToGraph(tmpMat2_, graph);

  //form tmpMat1_ = D^T*Kdd
  if ( D_->matTransMat(Kdd_, tmpMat1_) != 0) {
    throw std::runtime_error("fei::Reducer::assembleReducedGraph ERROR 3.");
  }

  //form tmpMat2_ = tpmMat1_*D = D^T*Kdd*D
  if ( tmpMat1_.matMat(*D_, tmpMat2_) != 0) {
    throw std::runtime_error("fei::Reducer::assembleReducedGraph ERROR 4.");
  }

  translateSSMatToReducedEqns(tmpMat2_);

  addSSMatToGraph(tmpMat2_, graph);

  //lastly, translate Kii and add it to the graph.
  translateSSMatToReducedEqns(Kii_);
  addSSMatToGraph(Kii_, graph);

  Kii_.logicalClear();
  Kdi_.logicalClear();
  Kid_.logicalClear();
  Kdd_.logicalClear();

  mat_counter_ = 0;

  if (global_gather) {
    graph->gatherFromOverlap();
  }
}

void
Reducer::assembleReducedGraph(fei::SparseRowGraph* srgraph)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"assembleReducedGraph(fei::SparseRowGraph)"<<FEI_ENDL;
  }

  fei::Graph_Impl graph(comm_, firstLocalReducedEqn_, lastLocalReducedEqn_);
  assembleReducedGraph(&graph);
  fei::copyToSparseRowGraph(*(graph.getLocalGraph()), *srgraph);
}

void
Reducer::assembleReducedMatrix(fei::Matrix& matrix)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"assembleReducedMatrix(fei::Matrix)"<<FEI_ENDL;
    if (output_level_ >= fei::FULL_LOGS) {
      os << dbgprefix_<<"Kid:"<<FEI_ENDL<<Kid_;
      os << dbgprefix_<<"Kdi:"<<FEI_ENDL<<Kdi_;
    }
  }

  //form tmpMat1_ = Kid_*D
  if ( Kid_.matMat(*D_, tmpMat1_, false) != 0) {
    throw std::runtime_error("fei::Reducer::assembleReducedMatrix ERROR 1.");
  }

  //form tmpMat2_ = D^T*Kdi_
  if ( D_->matTransMat(Kdi_, tmpMat2_, false) != 0) {
    throw std::runtime_error("fei::Reducer::assembleReducedMatrix ERROR 2.");
  }

  //accumulate the above two results into the global system matrix.
  translateSSMatToReducedEqns(tmpMat1_);
  addSSMatToMatrix(tmpMat1_, true, matrix);
  translateSSMatToReducedEqns(tmpMat2_);
  addSSMatToMatrix(tmpMat2_, true, matrix);

  //form tmpMat1_ = D^T*Kdd_
  if ( D_->matTransMat(Kdd_, tmpMat1_, false) != 0) {
    throw std::runtime_error("fei::Reducer::assembleReducedMatrix ERROR 3.");
  }

  //form tmpMat2_ = tmpMat1_*D = D^T*Kdd_*D
  if ( tmpMat1_.matMat(*D_, tmpMat2_, false) != 0) {
    throw std::runtime_error("fei::Reducer::assembleReducedMatrix ERROR 4.");
  }

  if (g_nonzero_) {
    //form tmpVec1_ = Kid_*g_
    if (Kid_.matVec(*g_, tmpVec1_) != 0) {
      throw std::runtime_error("fei::Reducer::assembleReducedMatrix ERROR 4g1.");
    }

    //add tmpVec1_ to fi_
    fi_.addEntries(tmpVec1_.length(), tmpVec1_.coefs().dataPtr(),
                   tmpVec1_.indices().dataPtr());

    //we already have tmpMat1_ = D^T*Kdd which was computed above, and we need
    //to form tmpVec1_ = D^T*Kdd*g_.
    //So we can simply form tmpVec1_ = tmpMat1_*g_.
    if (tmpMat1_.matVec(*g_, tmpVec1_) != 0) {
      throw std::runtime_error("fei::Reducer::assembleReducedMatrix ERROR 4g2.");
    }

    //now add tmpVec1_ to the right-hand-side fi_
    fi_.addEntries(tmpVec1_.length(), tmpVec1_.coefs().dataPtr(),
                   tmpVec1_.indices().dataPtr());
  }

  //accumulate tmpMat2_ = D^T*Kdd_*D into the global system matrix.
  translateSSMatToReducedEqns(tmpMat2_);
  addSSMatToMatrix(tmpMat2_, true, matrix);

  //lastly, translate Kii and add it to the graph.
  translateSSMatToReducedEqns(Kii_);
  addSSMatToMatrix(Kii_, true, matrix);

  Kii_.logicalClear();
  Kdi_.logicalClear();
  Kid_.logicalClear();
  Kdd_.logicalClear();

  mat_counter_ = 0;
}

bool
Reducer::isSlaveEqn(int unreducedEqn)
{
  int num = localUnreducedEqns_.size();

  int offset = num>0 ? unreducedEqn - localUnreducedEqns_[0] : -1;
  if (offset < 0 || offset >= (int)localUnreducedEqns_.size()) {
    return(isSlaveCol(unreducedEqn));
  }

  return(isSlaveEqn_[offset]);
}

std::vector<int>
Reducer::getSlaveMasterEqns(int slaveEqn)
{
  std::vector<int> masters;
  int idx = snl_fei::binarySearch(slaveEqn, slavesPtr_, numGlobalSlaves_);
  if (idx<0) {
    return(masters);
  }

  SSVec* row = D_->getRows()[idx];
  feiArray<int>& indices = row->indices();
  for(int i=0; i<indices.length(); ++i) {
    masters.push_back(indices[i]);
  }

  return(masters);
}

bool
Reducer::isSlaveCol(int unreducedEqn)
{
  int idx = snl_fei::binarySearch(unreducedEqn,
                                  slavesPtr_, numGlobalSlaves_);
  
  return(idx>=0);
}

int
Reducer::translateToReducedEqn(int eqn)
{
  if (eqn < lowestGlobalSlaveEqn_) {
    return(eqn);
  }

  if (eqn > highestGlobalSlaveEqn_) {
    return(eqn - numGlobalSlaves_);
  }

  int index = 0;
  int foundOffset = snl_fei::binarySearch(eqn, slavesPtr_, numGlobalSlaves_,
                                          index);

  if (foundOffset >= 0) {
    throw std::runtime_error("Reducer::translateToReducedEqn ERROR, input is slave eqn.");
  }

  return(eqn - index);
}

int
Reducer::translateFromReducedEqn(int reduced_eqn)
{
  int index = -1;
  int offset = snl_fei::binarySearch(reduced_eqn, &reverse_[0],
                                     reverse_.size(), index);
  if (offset >= 0) {
    return(nonslaves_[offset]);
  }

  if (index == 0) {
    return(reduced_eqn);
  }

  int adjustment = reduced_eqn - reverse_[index-1];

  return(nonslaves_[index-1] + adjustment);
}

void
Reducer::translateSSMatToReducedEqns(SSMat& mat)
{
  //Given a matrix in global numbering, convert all of its contents to the
  //"reduced" equation space. If any of the row-numbers or column-indices in
  //this matrix object are slave equations, an exception will be thrown.

  feiArray<int>& rowNumbers = mat.getRowNumbers();
  SSVec** rows = mat.getRows().dataPtr();

  for(int i=0; i<rowNumbers.length(); i++) {
    rowNumbers[i] = translateToReducedEqn(rowNumbers[i]);

    feiArray<int>& row_indices = rows[i]->indices();
    int* inds = row_indices.dataPtr();
    int numInds = row_indices.length();
    for(int j=0; j<numInds; j++) {
      inds[j] = translateToReducedEqn(inds[j]);
    }
  }
}

void
Reducer::translateSSVecToReducedEqns(SSVec& vec)
{
  feiArray<int>& indices = vec.indices();
  int* indPtr = indices.dataPtr();
  for(int i=0; i<indices.length(); ++i) {
    indPtr[i] = translateToReducedEqn(indPtr[i]);
  }
}

void
Reducer::addSSMatToGraph(SSMat& mat, fei::Graph* graph)
{
   //This function must be called with a SSMat object that already has its
  //contents numbered in "reduced" equation-numbers.
  //
  //This function has one simple task -- for each row,col pair stored in 'mat',
  //add the index pair to the graph_ object.
  //

  int numRows = mat.getRowNumbers().length();
  if (numRows == 0) return;

  int* rowNumbers = mat.getRowNumbers().dataPtr();
  SSVec** rows = mat.getRows().dataPtr();

  for(int i=0; i<numRows; i++) {
    int rowLen = rows[i]->length();
    int* indicesRow = rows[i]->indices().dataPtr();

    if ( graph->addIndices(rowNumbers[i], rowLen, indicesRow) != 0) {
      throw std::runtime_error("fei::Reducer::addSSMatToGraph ERROR in graph->addIndices.");
    }
  }
}

void
Reducer::addSSMatToMatrix(SSMat& mat, bool sum_into,
                             fei::Matrix& matrix)
{
  int numRows = mat.getRowNumbers().length();
  if (numRows == 0) return;

  int* rowNumbers = mat.getRowNumbers().dataPtr();
  SSVec** rows = mat.getRows().dataPtr();

  for(int i=0; i<numRows; i++) {
    int rowLen = rows[i]->length();
    int* indicesRow = rows[i]->indices().dataPtr();

    double* coefPtr = rows[i]->coefs().dataPtr();
    if (sum_into) {
      if ( matrix.sumIn(1, &rowNumbers[i],
                         rowLen, indicesRow,
                         &coefPtr, FEI_DENSE_ROW) != 0) {
        throw std::runtime_error("fei::Reducer::addSSMatToMatrix ERROR in matrix->sumIn.");
      }
    }
    else {
      if ( matrix.copyIn(1, &rowNumbers[i],
                          rowLen, indicesRow,
                          &coefPtr, FEI_DENSE_ROW) != 0) {
        throw std::runtime_error("fei::Reducer::addSSMatToMatrix ERROR in matrix->copyIn.");
      }
    }
  }
}

int
Reducer::addMatrixValues(int numRows, const int* rows,
                         int numCols, const int* cols,
                         const double* const* values,
                         bool sum_into)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addMatrixValues()"<<FEI_ENDL;
  }

  expand_work_arrays(numCols);

  for(int j=0; j<numCols; ++j) {
    bool_array_[j] = isSlaveEqn(cols[j]);
  }

  for(int i=0; i<numRows; ++i) {
    if (isSlaveEqn(rows[i])) {
      //slave row: slave columns go into Kdd, non-slave columns go
      //into Kdi.
      SSVec* Kdd_row = Kdd_.getRow(rows[i], true);
      SSVec* Kdi_row = Kdi_.getRow(rows[i], true);

      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          Kdd_row->addEntry(cols[j], values[i][j]);
        }
        else {
          Kdi_row->addEntry(cols[j], values[i][j]);
        }
      }
    }
    else {//not slave row
      //put non-slave columns into Kii,
      //and slave columns into Kid.
      SSVec* Kid_row = Kid_.getRow(rows[i], true);
      SSVec* Kii_row = Kii_.getRow(rows[i], true);

      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          Kid_row->addEntry(cols[j], values[i][j]);
        }
        else {
          Kii_row->addEntry(cols[j], values[i][j]);
        }
      }
    }
  }

  return(0);
}

int
Reducer::addMatrixValues(int numRows, const int* rows,
                         int numCols, const int* cols,
                         const double* const* values,
                         bool sum_into,
                         fei::Matrix& feimat,
                         int format)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addMatrixValues(fei::Matrix)"<<FEI_ENDL;
  }

  expand_work_arrays(numCols+numRows);

  const double** myvalues = const_cast<const double**>(values);
  if (format != FEI_DENSE_ROW) {
    if (format != FEI_DENSE_COL) {
      throw std::runtime_error("fei::Reducer::addMatrixValues ERROR, submatrix format must be either FEI_DENSE_ROW or FEI_DENSE_COL. Other formats not supported with slave constraints.");
    }

    fei::Matrix_core::copyTransposeToWorkArrays(numRows, numCols, values,
                                                work_1D_, work_2D_);
    myvalues = &work_2D_[0];
  }

  bool no_slave_cols = true;
  unsigned num_non_slave_cols = 0;
  for(int j=0; j<numCols; ++j) {
    bool_array_[j] = isSlaveEqn(cols[j]);
    if (bool_array_[j]) no_slave_cols = false;
    else int_array_[num_non_slave_cols++] = translateToReducedEqn(cols[j]);
  }

  bool no_slave_rows = true;
  for(int i=0; i<numRows; ++i) {
    bool_array_[numCols+i] = isSlaveEqn(rows[i]);
    if (bool_array_[numCols+i]) no_slave_rows = false;
    else int_array_[numCols+i] = translateToReducedEqn(rows[i]);
  }

  if (no_slave_rows && no_slave_cols) {
    if (sum_into) {
      feimat.sumIn(numRows, int_array_+numCols,
                   numCols, int_array_, myvalues, FEI_DENSE_ROW);
    }
    else {
      feimat.copyIn(numRows, int_array_+numCols,
                    numCols, int_array_, myvalues, FEI_DENSE_ROW);
    }

    return(0);
  }

  for(int i=0; i<numRows; ++i) {
    if (bool_array_[numCols+i]) {
      //slave row: slave columns go into Kdd, non-slave columns go
      //into Kdi.
      SSVec* Kdd_row = Kdd_.getRow(rows[i], true);
      SSVec* Kdi_row = Kdi_.getRow(rows[i], true);

      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          Kdd_row->addEntry(cols[j], myvalues[i][j]);
        }
        else {
          Kdi_row->addEntry(cols[j], myvalues[i][j]);
        }
      }
      ++mat_counter_;
    }
    else {//not slave row
      if (no_slave_cols) {
        int reduced_row = int_array_[numCols+i];
        const double* rowvals = myvalues[i];
        if (sum_into) {
          feimat.sumIn(1, &reduced_row, numCols, int_array_,
                       &rowvals, format);
        }
        else {
          feimat.copyIn(1, &reduced_row, num_non_slave_cols, int_array_,
                        &rowvals, format);
        }
        continue;
      }

      //put non-slave columns into Kii,
      //and slave columns into Kid.
      SSVec* Kid_row = Kid_.getRow(rows[i], true);

      unsigned offset = 0;
      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          Kid_row->addEntry(cols[j], myvalues[i][j]);
          ++mat_counter_;
        }
        else {
          double_array_[offset++] = myvalues[i][j];
        }
      }

      if (num_non_slave_cols > 0) {
        int reduced_row = int_array_[numCols+i];
        if (sum_into) {
          feimat.sumIn(1, &reduced_row, num_non_slave_cols, int_array_,
                       &double_array_, format);
        }
        else {
          feimat.copyIn(1, &reduced_row, num_non_slave_cols, int_array_,
                        &double_array_, format);
        }
      }
    }
  }

  if (mat_counter_ > 600) {
    assembleReducedMatrix(feimat);
  }

  return(0);
}

int
Reducer::addVectorValues(int numValues,
                         const int* globalIndices,
                         const double* values,
                         bool sum_into,
                         bool soln_vector,
                         int vectorIndex)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addVectorValues()"<<FEI_ENDL;
  }

  for(int i=0; i<numValues; ++i) {
    if (isSlaveEqn(globalIndices[i])) {
      if (sum_into) {
        if (!soln_vector) fd_.addEntry(globalIndices[i], values[i]);
        //if soln_vector==true we discard values for slave eqns...
      }
      else {
        if (!soln_vector) fd_.putEntry(globalIndices[i], values[i]);
        //if soln_vector==true we discard values for slave eqns...
      }
      if (!soln_vector) ++rhs_vec_counter_;
    }
    else {
      if (sum_into) {
        if (soln_vector) xi_.addEntry(globalIndices[i], values[i]);
        else fi_.addEntry(globalIndices[i], values[i]);
      }
      else {
        if (soln_vector) xi_.putEntry(globalIndices[i], values[i]);
        else fi_.putEntry(globalIndices[i], values[i]);
      }
    }
  }

  return(0);
}

int
Reducer::addVectorValues(int numValues,
                         const int* globalIndices,
                         const double* values,
                         bool sum_into,
                         bool soln_vector,
                         int vectorIndex,
                         fei::Vector& feivec)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"addVectorValues(fei::Vector)"<<FEI_ENDL;
  }

  for(int i=0; i<numValues; ++i) {
    if (isSlaveEqn(globalIndices[i])) {
      if (sum_into) {
        if (!soln_vector) fd_.addEntry(globalIndices[i], values[i]);
        //if soln_vector==true we discard values for slave eqns...
      }
      else {
        if (!soln_vector) fd_.putEntry(globalIndices[i], values[i]);
        //if soln_vector==true we discard values for slave eqns...
      }
      if (!soln_vector) ++rhs_vec_counter_;
    }
    else {
      int reduced_index = translateToReducedEqn(globalIndices[i]);

      if (sum_into) {
        feivec.sumIn(1, &reduced_index, &values[i], vectorIndex);
      }
      else {
        feivec.copyIn(1, &reduced_index, &values[i], vectorIndex);
      }
    }
  }

  if (soln_vector && soln_vec_counter_ > 600) {
    assembleReducedVector(soln_vector, feivec, sum_into);
  }
  else if (!soln_vector && rhs_vec_counter_ > 600) {
    assembleReducedVector(soln_vector, feivec, sum_into);
  }

  return(0);
}

void
Reducer::assembleReducedVector(bool soln_vector,
                               fei::Vector& feivec,
                               bool sum_into)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"assembleReducedVector(fei::Vector)"<<FEI_ENDL;
  }

  SSVec& vec = soln_vector ? xd_ : fd_;

  if (!soln_vector && vec.length() > 0) {
    //form tmpVec1 = D^T*vec.
    if ( D_->matTransVec(vec, tmpVec1_) != 0) {
      throw std::runtime_error("fei::Reducer::assembleReducedVec ERROR.");
    }

    vec.logicalClear();

    translateSSVecToReducedEqns(tmpVec1_);

    if (sum_into) {
      feivec.sumIn(tmpVec1_.length(), tmpVec1_.indices().dataPtr(),
                   tmpVec1_.coefs().dataPtr(), 0);
    }
    else {
      feivec.copyIn(tmpVec1_.length(), tmpVec1_.indices().dataPtr(),
                   tmpVec1_.coefs().dataPtr(), 0);
    }
  }

  SSVec& vec_i = soln_vector ? xi_ : fi_;

  if (vec_i.length() > 0) {
    translateSSVecToReducedEqns(vec_i);

    if (sum_into) {
      feivec.sumIn(vec_i.length(), vec_i.indices().dataPtr(),
                   vec_i.coefs().dataPtr(), 0);
    }
    else {
      feivec.copyIn(vec_i.length(), vec_i.indices().dataPtr(),
                    vec_i.coefs().dataPtr(), 0);
    }

    vec_i.logicalClear();
  }

  if (soln_vector) {
    soln_vec_counter_ = 0;
  }
  else {
    rhs_vec_counter_ = 0;
  }
}

int
Reducer::copyOutVectorValues(int numValues,
                             const int* globalIndices,
                             double* values,
                             bool soln_vector,
                             int vectorIndex,
                             fei::Vector& feivec)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"copyOutVectorValues"<<FEI_ENDL;
  }

  tmpVec1_.logicalClear();
  std::vector<int> reduced_indices;
  std::vector<int> offsets;

  for(int i=0; i<numValues; ++i) {
    if (isSlaveEqn(globalIndices[i])) {
      tmpVec1_.putEntry(globalIndices[i], 0.0);
      offsets.push_back(i);
    }
    else {
      int reduced_idx = translateToReducedEqn(globalIndices[i]);
      feivec.copyOut(1, &reduced_idx, &values[i], vectorIndex);
    }
  }

  if (tmpVec1_.length() > 0) {
    D_->matTransVec(tmpVec1_, tmpVec2_);
    int* tmpVec2Indices = tmpVec2_.indices().dataPtr();
    for(int i=0; i<tmpVec2_.length(); ++i) {
      reduced_indices.push_back(translateToReducedEqn(tmpVec2Indices[i]));
    }

    feivec.copyOut(tmpVec2_.length(), &reduced_indices[0],
                   tmpVec2_.coefs().dataPtr(), vectorIndex);

    D_->matVec(tmpVec2_, tmpVec1_);

    if (g_nonzero_) {
      int* ginds = g_->indices().dataPtr();
      double* gcoefs = g_->coefs().dataPtr();
      for(int ii=0; ii<g_->length(); ++ii) {
        tmpVec1_.addEntry(ginds[ii], -gcoefs[ii]);
      }
    }

    int len = tmpVec1_.length();
    int* indices = tmpVec1_.indices().dataPtr();
    double* coefs = tmpVec1_.coefs().dataPtr();

    for(unsigned ii=0; ii<offsets.size(); ++ii) {
      int index = globalIndices[offsets[ii]];
      int idx = snl_fei::binarySearch(index, indices, len);
      if (idx < 0) {
        continue;
      }

      values[offsets[ii]] = coefs[idx];
    }
  }

  return(0);
}

std::vector<int>&
Reducer::getLocalReducedEqns()
{
  return(localReducedEqns_);
}

}//namespace fei

