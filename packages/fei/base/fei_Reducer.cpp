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



#include <fei_Reducer.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix.hpp>
#include <fei_Matrix_core.hpp>
#include <fei_Vector.hpp>
#include <fei_Graph_Impl.hpp>
#include <fei_ArrayUtils.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_SparseRowGraph.hpp>
#include <fei_Vector.hpp>
#include <fei_impl_utils.hpp>

namespace fei {

Reducer::Reducer(fei::SharedPtr<FillableMat> globalSlaveDependencyMatrix,
                 fei::SharedPtr<CSVec> g_vector,
                 MPI_Comm comm)
 : csrD_(),
   slavesPtr_(NULL),
   Kii_(),
   Kid_(),
   Kdi_(),
   Kdd_(),
   csrKii(),
   csrKid(),
   csrKdi(),
   csrKdd(),
   fi_(),
   fd_(),
   csfi(),
   csvec(),
   csvec_i(),
   tmpMat1_(),
   tmpMat2_(),
   tmpVec1_(),
   tmpVec2_(),
   csg_(),
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
   rhs_vec_counter_(0),
   bool_array_(0),
   int_array_(0),
   double_array_(0),
   array_len_(0),
   work_1D_(),
   work_2D_()
{
  csrD_ = *globalSlaveDependencyMatrix;
  if (g_vector.get() != NULL) {
    csg_ = *g_vector;
  }

  initialize();
}

void
Reducer::initialize()
{
  numGlobalSlaves_ = csrD_.getNumRows();
  slavesPtr_ = &((csrD_.getGraph().rowNumbers)[0]);
  lowestGlobalSlaveEqn_ = slavesPtr_[0];
  highestGlobalSlaveEqn_ = slavesPtr_[numGlobalSlaves_-1];

  if (csg_.size() > 0) {
    double* gptr = &(csg_.coefs()[0]);
    for(size_t i=0; i<csg_.size(); ++i) {
      if (gptr[i] != 0.0) {
        g_nonzero_ = true;
        break;
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
 : csrD_(),
   slavesPtr_(NULL),
   Kii_(),
   Kid_(),
   Kdi_(),
   Kdd_(),
   fi_(),
   fd_(),
   tmpMat1_(),
   tmpMat2_(),
   tmpVec1_(),
   tmpVec2_(),
   csg_(),
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
   rhs_vec_counter_(0),
   bool_array_(0),
   int_array_(0),
   double_array_(0),
   array_len_(0)
{
  fei::SharedPtr<fei::VectorSpace> vecSpace = matrixGraph->getRowSpace();
  comm_ = vecSpace->getCommunicator();
  initialize();

  std::vector<int> indices;
  vecSpace->getIndices_Owned(indices);
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

  if (isSlaveEqn_ != NULL) delete [] isSlaveEqn_;

  isSlaveEqn_ = num > 0 ? new bool[localUnreducedEqns_.size()] : NULL;

  numLocalSlaves_ = 0;

  for(size_t i=0; i<localUnreducedEqns_.size(); ++i) {
    int idx = fei::binarySearch(localUnreducedEqns_[i],
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
        int offset = csrD_.getGraph().rowOffsets[idx];
        int rowlen = csrD_.getGraph().rowOffsets[idx+1]-offset;
        int* indices = &(csrD_.getGraph().packedColumnIndices[offset]);
        for(int j=0; j<rowlen; ++j) {
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

  std::vector<int>& rowNumbers = matrixGraph->rowNumbers;
  std::vector<int>& rowOffsets = matrixGraph->rowOffsets;
  std::vector<int>& packedCols = matrixGraph->packedColumnIndices;

  for(unsigned i=0; i<rowNumbers.size(); ++i) {
    int row = rowNumbers[i];

    bool slave_row = isSlaveEqn(row);

    int rowLength = rowOffsets[i+1]-rowOffsets[i];
    int* cols = &packedCols[rowOffsets[i]];

    if (slave_row) {
      fei::CSVec* Kdd_row = Kdd_.create_or_getRow(row);
      fei::CSVec* Kdi_row = Kdi_.create_or_getRow(row);

      for(int j=0; j<rowLength; ++j) {
        int col = cols[j];
        bool slave_col = isSlaveEqn(col);

        if (slave_col) {
          add_entry(*Kdd_row, col, 0.0);
        }
        else {
          add_entry(*Kdi_row, col, 0.0);
        }
      }
    }
    else {
      //not a slave row, so add slave columns to Kid, and non-slave
      //columns to graph.
      fei::CSVec* Kid_row = Kid_.create_or_getRow(row);
      fei::CSVec* Kii_row = Kii_.create_or_getRow(row);

      for(int j=0; j<rowLength; ++j) {
        int col = cols[j];
        bool slave_col = isSlaveEqn(col);

        if (slave_col) {
          add_entry(*Kid_row, col, 0.0);
        }
        else {
          add_entry(*Kii_row, col, 0.0);
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
  expand_work_arrays(numCols);

  bool no_slave_cols = true;
  for(int i=0; i<numCols; ++i) {
    bool_array_[i] = isSlaveEqn(cols[i]);
    if (bool_array_[i]) no_slave_cols = false;
  }

  for(int i=0; i<numRows; ++i) {
    bool slave_row = isSlaveEqn(rows[i]);

    if (slave_row) {
      fei::CSVec* Kdd_row = Kdd_.create_or_getRow(rows[i]);
      fei::CSVec* Kdi_row = Kdi_.create_or_getRow(rows[i]);

      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          add_entry(*Kdd_row, cols[j], 0.0);
        }
        else {
          add_entry(*Kdi_row, cols[j], 0.0);
        }
      }
      ++mat_counter_;
    }
    else {
      //not a slave row, so add slave columns to Kid, and non-slave
      //columns to graph.
      fei::CSVec* Kid_row = no_slave_cols ?
        NULL : Kid_.create_or_getRow(rows[i]);
  
      unsigned num_non_slave_cols = 0;

      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          add_entry(*Kid_row, cols[j], 0.0);
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
  csrKid = Kid_;
  fei::multiply_CSRMat_CSRMat(csrKid, csrD_, tmpMat1_, true);

  csrKdi = Kdi_;
  fei::multiply_trans_CSRMat_CSRMat(csrD_, csrKdi, tmpMat2_, true);

  fei::impl_utils::translate_to_reduced_eqns(*this, tmpMat1_);
  fei::impl_utils::translate_to_reduced_eqns(*this, tmpMat2_);

  fei::impl_utils::add_to_graph(tmpMat1_, *graph);
  fei::impl_utils::add_to_graph(tmpMat2_, *graph);

  //form tmpMat1_ = D^T*Kdd
  csrKdd = Kdd_;
  fei::multiply_trans_CSRMat_CSRMat(csrD_, csrKdd, tmpMat1_, true);

  //form tmpMat2_ = tmpMat1_*D = D^T*Kdd*D
  fei::multiply_CSRMat_CSRMat(tmpMat1_, csrD_, tmpMat2_, true);

  fei::impl_utils::translate_to_reduced_eqns(*this, tmpMat2_);

  fei::impl_utils::add_to_graph(tmpMat2_, *graph);

  //lastly, translate Kii and add it to the graph.
  csrKii = Kii_;
  fei::impl_utils::translate_to_reduced_eqns(*this, csrKii);
  fei::impl_utils::add_to_graph(csrKii, *graph);

  Kii_.clear();
  Kdi_.clear();
  Kid_.clear();
  Kdd_.clear();

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
  csrKid = Kid_;

  fei::multiply_CSRMat_CSRMat(csrKid, csrD_, tmpMat1_);

  //form tmpMat2_ = D^T*Kdi_
  csrKdi = Kdi_;
  fei::multiply_trans_CSRMat_CSRMat(csrD_, csrKdi, tmpMat2_);

  //accumulate the above two results into the global system matrix.
  fei::impl_utils::translate_to_reduced_eqns(*this, tmpMat1_);

  fei::impl_utils::add_to_matrix(tmpMat1_, true, matrix);

  fei::impl_utils::translate_to_reduced_eqns(*this, tmpMat2_);
  fei::impl_utils::add_to_matrix(tmpMat2_, true, matrix);

  //form tmpMat1_ = D^T*Kdd_
  csrKdd = Kdd_;
  fei::multiply_trans_CSRMat_CSRMat(csrD_, csrKdd, tmpMat1_);

  //form tmpMat2_ = tmpMat1_*D = D^T*Kdd_*D
  fei::multiply_CSRMat_CSRMat(tmpMat1_, csrD_, tmpMat2_);

  if (g_nonzero_) {
    //form tmpVec1_ = Kid_*g_
    fei::multiply_CSRMat_CSVec(csrKid, csg_, tmpVec1_);

    //add tmpVec1_ to fi_
    csfi = fi_;
    fei::add_CSVec_CSVec(tmpVec1_, csfi);

    //we already have tmpMat1_ = D^T*Kdd which was computed above, and we need
    //to form tmpVec1_ = D^T*Kdd*g_.
    //So we can simply form tmpVec1_ = tmpMat1_*g_.
    fei::multiply_CSRMat_CSVec(tmpMat1_, csg_, tmpVec1_);

    //now add tmpVec1_ to the right-hand-side fi_
    fei::add_CSVec_CSVec(tmpVec1_, csfi);
  }

  //accumulate tmpMat2_ = D^T*Kdd_*D into the global system matrix.
  fei::impl_utils::translate_to_reduced_eqns(*this, tmpMat2_);
  fei::impl_utils::add_to_matrix(tmpMat2_, true, matrix);

  //lastly, translate Kii and add it to the graph.
  csrKii = Kii_;
  fei::impl_utils::translate_to_reduced_eqns(*this, csrKii);
  fei::impl_utils::add_to_matrix(csrKii, true, matrix);

  Kii_.clear();
  Kdi_.clear();
  Kid_.clear();
  Kdd_.clear();

  mat_counter_ = 0;
}

bool
Reducer::isSlaveEqn(int unreducedEqn) const
{
  int num = localUnreducedEqns_.size();

  int offset = num>0 ? unreducedEqn - localUnreducedEqns_[0] : -1;
  if (offset < 0 || offset >= (int)localUnreducedEqns_.size()) {
    return(isSlaveCol(unreducedEqn));
  }

  return(isSlaveEqn_[offset]);
}

void
Reducer::getSlaveMasterEqns(int slaveEqn, std::vector<int>& masterEqns)
{
  masterEqns.clear();

  std::vector<int>& rows = csrD_.getGraph().rowNumbers;

  std::vector<int>::iterator iter =
    std::lower_bound(rows.begin(), rows.end(), slaveEqn);

  if (iter == rows.end() || *iter != slaveEqn) {
    return;
  }

  size_t offset = iter - rows.begin();

  int rowBegin = csrD_.getGraph().rowOffsets[offset];
  int rowEnd = csrD_.getGraph().rowOffsets[offset+1];
  std::vector<int>& cols = csrD_.getGraph().packedColumnIndices;

  for(int j=rowBegin; j<rowEnd; ++j) {
    masterEqns.push_back(cols[j]);
  }
}

bool
Reducer::isSlaveCol(int unreducedEqn) const
{
  int idx = fei::binarySearch(unreducedEqn,
                                  slavesPtr_, numGlobalSlaves_);
  
  return(idx>=0);
}

int
Reducer::translateToReducedEqn(int eqn) const
{
  if (eqn < lowestGlobalSlaveEqn_) {
    return(eqn);
  }

  if (eqn > highestGlobalSlaveEqn_) {
    return(eqn - numGlobalSlaves_);
  }

  int index = 0;
  int foundOffset = fei::binarySearch(eqn, slavesPtr_, numGlobalSlaves_,
                                          index);

  if (foundOffset >= 0) {
    throw std::runtime_error("Reducer::translateToReducedEqn ERROR, input is slave eqn.");
  }

  return(eqn - index);
}

int
Reducer::translateFromReducedEqn(int reduced_eqn) const
{
  int index = -1;
  int offset = fei::binarySearch(reduced_eqn, &reverse_[0],
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
      fei::CSVec* Kdd_row = Kdd_.create_or_getRow(rows[i]);
      fei::CSVec* Kdi_row = Kdi_.create_or_getRow(rows[i]);

      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          add_entry(*Kdd_row, cols[j], myvalues[i][j]);
        }
        else {
          add_entry(*Kdi_row, cols[j], myvalues[i][j]);
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
      fei::CSVec* Kid_row = Kid_.create_or_getRow(rows[i]);

      unsigned offset = 0;
      for(int j=0; j<numCols; ++j) {
        if (bool_array_[j]) {
          add_entry(*Kid_row, cols[j], myvalues[i][j]);
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
      }
      else {
        if (!soln_vector) fd_.putEntry(globalIndices[i], values[i]);
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

  if (rhs_vec_counter_ > 600) {
    assembleReducedVector(soln_vector, feivec);
  }

  return(0);
}

void
Reducer::assembleReducedVector(bool soln_vector,
                               fei::Vector& feivec)
{
  if (output_level_ >= fei::BRIEF_LOGS && output_stream_ != NULL) {
    FEI_OSTREAM& os = *output_stream_;
    os << dbgprefix_<<"assembleReducedVector(fei::Vector)"<<FEI_ENDL;
  }

  if (soln_vector) {
    return;
  }

  fei::FillableVec& vec = fd_;

  if (vec.size() > 0) {
    //form tmpVec1 = D^T*vec.
    csvec = vec;
    fei::multiply_trans_CSRMat_CSVec(csrD_, csvec, tmpVec1_);

    vec.clear();

    fei::impl_utils::translate_to_reduced_eqns(*this, tmpVec1_);

    int which_vector = 0;
    feivec.sumIn(tmpVec1_.size(), &(tmpVec1_.indices()[0]),
                 &(tmpVec1_.coefs()[0]), which_vector);
  }

  fei::FillableVec& vec_i = fi_;

  if (vec_i.size() > 0) {
    csvec_i = vec_i;
    fei::impl_utils::translate_to_reduced_eqns(*this, csvec_i);

    int which_vector = 0;
    feivec.sumIn(csvec_i.size(), &(csvec_i.indices()[0]),
                 &(csvec_i.coefs()[0]), which_vector);

    vec_i.clear();
  }

  rhs_vec_counter_ = 0;
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

  tmpVec1_.clear();
  tmpVec2_.clear();
  std::vector<int> reduced_indices;
  std::vector<int> offsets;

  for(int i=0; i<numValues; ++i) {
    if (isSlaveEqn(globalIndices[i])) {
      fei::put_entry(tmpVec1_, globalIndices[i], 0.0);
      offsets.push_back(i);
    }
    else {
      int reduced_idx = translateToReducedEqn(globalIndices[i]);
      feivec.copyOut(1, &reduced_idx, &values[i], vectorIndex);
    }
  }

  if (tmpVec1_.size() > 0) {
    fei::multiply_trans_CSRMat_CSVec(csrD_, tmpVec1_, tmpVec2_);
    int* tmpVec2Indices = &(tmpVec2_.indices()[0]);
    for(size_t i=0; i<tmpVec2_.size(); ++i) {
      reduced_indices.push_back(translateToReducedEqn(tmpVec2Indices[i]));
    }

    feivec.copyOut(tmpVec2_.size(), &reduced_indices[0],
                   &(tmpVec2_.coefs()[0]), vectorIndex);

    fei::multiply_CSRMat_CSVec(csrD_, tmpVec2_, tmpVec1_);

    if (g_nonzero_) {
      int* ginds = &(csg_.indices()[0]);
      double* gcoefs = &(csg_.coefs()[0]);
      for(size_t ii=0; ii<csg_.size(); ++ii) {
        fei::add_entry(tmpVec1_, ginds[ii], -gcoefs[ii]);
      }
    }

    int len = tmpVec1_.size();
    int* indices = &(tmpVec1_.indices()[0]);
    double* coefs = &(tmpVec1_.coefs()[0]);

    for(unsigned ii=0; ii<offsets.size(); ++ii) {
      int index = globalIndices[offsets[ii]];
      int idx = fei::binarySearch(index, indices, len);
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

