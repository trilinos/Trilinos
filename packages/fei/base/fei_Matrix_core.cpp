/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <fei_macros.hpp>

#include <fei_utils.hpp>

#include <fei_ParameterSet.hpp>
#include <fei_LogManager.hpp>
#include <fei_Matrix_core.hpp>
#include <fei_CSRMat.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_CommUtils.hpp>
#include <fei_impl_utils.hpp>

#include <fei_VectorSpace.hpp>
#include <snl_fei_PointBlockMap.hpp>
#include <fei_MatrixGraph.hpp>
#include <snl_fei_Utils.hpp>

#include <fei_MatrixTraits.hpp>
#include <fei_MatrixTraits_FillableMat.hpp>

#undef fei_file
#define fei_file "fei_Matrix_core.cpp"
#include <fei_ErrMacros.hpp>

fei::Matrix_core::Matrix_core(fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                              int numLocalEqns)
  : name_(),
    work_indices_(),
    work_indices2_(),
    work_ints_(),
    work_data1D_(),
    work_data2D_(),
    eqnComm_(),
    rhsVector_(),
    comm_(matrixGraph->getRowSpace()->getCommunicator()),
    localProc_(0),
    numProcs_(1),
    vecSpace_(),
    matrixGraph_(matrixGraph),
    remotelyOwned_(),
    remotelyOwned_last_requested_(NULL),
    sendProcs_(),
    recvProcs_(),
    recv_chars_(),
    send_chars_(),
    sendRecvProcsNeedUpdated_(true),
    proc_last_requested_(-1),
    haveBlockMatrix_(false),
    haveFEMatrix_(false),
    globalOffsets_(),
    firstLocalOffset_(0),
    lastLocalOffset_(0)
{
  if (matrixGraph.get() == NULL) {
    throw std::runtime_error("fei::Matrix_core constructed with NULL fei::MatrixGraph");
  }

  vecSpace_ = matrixGraph->getRowSpace();

  if (vecSpace_->initCompleteAlreadyCalled()) {
    vecSpace_->getGlobalIndexOffsets(globalOffsets_);
    eqnComm_.reset(new fei::EqnComm(comm_, numLocalEqns, globalOffsets_));
  }
  else {
    eqnComm_.reset(new fei::EqnComm(comm_, numLocalEqns));
  }

  localProc_ = fei::localProc(comm_);
  numProcs_ = fei::numProcs(comm_);

  setName("dbg");

  globalOffsets_ = eqnComm_->getGlobalOffsets();

  firstLocalOffset_ = globalOffsets_[localProc_];
  lastLocalOffset_ = globalOffsets_[localProc_+1]-1;
}

std::map<int,fei::FillableMat*>&
fei::Matrix_core::getRemotelyOwnedMatrices()
{
  return remotelyOwned_;
}

void
fei::Matrix_core::putScalar_remotelyOwned(double scalar)
{
  std::map<int,FillableMat*>::iterator
    it = remotelyOwned_.begin(),
    it_end = remotelyOwned_.end();
  for(; it!=it_end; ++it) {
    fei::MatrixTraits<FillableMat>::setValues(it->second, scalar);
  }
}

void
fei::Matrix_core::setEqnComm(fei::SharedPtr<fei::EqnComm> eqnComm)
{
  eqnComm_ = eqnComm;
  globalOffsets_ = eqnComm_->getGlobalOffsets();
  firstLocalOffset_ = globalOffsets_[localProc_];
  lastLocalOffset_ = globalOffsets_[localProc_+1]-1;
  sendRecvProcsNeedUpdated_ = true;
}

fei::Matrix_core::~Matrix_core()
{
  std::map<int,FillableMat*>::iterator
    it = remotelyOwned_.begin(),
    it_end = remotelyOwned_.end();
  for(; it!=it_end; ++it) {
    delete it->second;
  }
}

void
fei::Matrix_core::parameters(const fei::ParameterSet& paramset)
{
  const fei::Param* param = paramset.get("name");
  fei::Param::ParamType ptype = param != NULL ?
    param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    setName(param->getStringValue().c_str());
  }

  param = paramset.get("FEI_OUTPUT_LEVEL");
  ptype = param != NULL ? param->getType() : fei::Param::BAD_TYPE;
  if (ptype == fei::Param::STRING) {
    fei::LogManager& log_manager = fei::LogManager::getLogManager();
    log_manager.setOutputLevel(param->getStringValue().c_str());
    setOutputLevel(fei::utils::string_to_output_level(param->getStringValue()));
  }

}

void fei::Matrix_core::setName(const char* name)
{
  if (name == NULL) return;

  if (name_ == name) return;

  name_ = name;
}

int fei::Matrix_core::getOwnerProc(int globalEqn) const
{
  int len = globalOffsets_.size();
  if (globalEqn > globalOffsets_[len-1]) return(-1);

  for(int p=len-2; p>=0; --p) {
    if (globalEqn >= globalOffsets_[p]) {
      return(p);
    }
  }

  return(-1);
}

void fei::Matrix_core::setRHS(fei::SharedPtr<fei::Vector> rhsvector)
{
  rhsVector_ = rhsvector;
}

void fei::Matrix_core::setCommSizes()
{
#ifndef FEI_SER
  //iterate the remotelyOwned_ map and create a list of processors
  //that we will be sending data to. (processors which own matrix rows
  //that we share.)
  sendProcs_.clear();
  std::map<int,FillableMat*>::iterator
    it = remotelyOwned_.begin(),
    it_end = remotelyOwned_.end();
  for(; it!=it_end; ++it) {
    if (it->first == localProc_) continue;
    if (it->second != NULL) {
      if (it->second->getNumRows() == 0) {
        continue;
      }

      sendProcs_.push_back(it->first);
    }
  }

  //now we can find out which procs we'll be receiving from.
  recvProcs_.clear();
  fei::mirrorProcs(comm_, sendProcs_, recvProcs_);

  int tag1 = 11111;
  std::vector<MPI_Request> mpiReqs(recvProcs_.size());

  //next find out the size of each message that will be recvd or sent.
  std::vector<int> recv_sizes(recvProcs_.size(), 0);
  for(size_t i=0; i<recvProcs_.size(); ++i) {
    MPI_Irecv(&recv_sizes[i], 1, MPI_INT, recvProcs_[i],
              tag1, comm_, &mpiReqs[i]);
  }

  //pack our to-be-sent data into buffers, and send the
  //sizes to the receiving procs:
  send_chars_.resize(sendProcs_.size());
  recv_chars_.resize(recvProcs_.size());

  for(size_t i=0; i<sendProcs_.size(); ++i) {
    int proc = sendProcs_[i];

    int num_bytes = fei::impl_utils::num_bytes_FillableMat(*(remotelyOwned_[proc]));
    send_chars_[i].resize(num_bytes);
    char* buffer = &(send_chars_[i][0]);
    fei::impl_utils::pack_FillableMat(*(remotelyOwned_[proc]), buffer);

    int bsize = send_chars_[i].size();

    MPI_Send(&bsize, 1, MPI_INT, proc, tag1, comm_);
  }

  int numRecvProcs = recvProcs_.size();
  for(size_t i=0; i<recvProcs_.size(); ++i) {
    MPI_Status status;
    int index;
    MPI_Waitany(numRecvProcs, &mpiReqs[0], &index, &status);

    recv_chars_[index].resize(recv_sizes[index]);
  }

  sendRecvProcsNeedUpdated_ = false;
#endif
}

int fei::Matrix_core::gatherFromOverlap(bool accumulate)
{
  if (numProcs() == 1) return(0);

#ifndef FEI_SER
  //this function gathers shared-but-not-owned data onto the
  //owning processors.

  if (sendRecvProcsNeedUpdated_) {
    setCommSizes();
  }

  std::vector<MPI_Request> mpiReqs(recvProcs_.size());

  int tag1 = 11111;
  int numRecvProcs = recvProcs_.size();

  //post the recvs.
  for(size_t i=0; i<recvProcs_.size(); ++i) {
    int bsize = recv_chars_[i].size();

    MPI_Irecv(&(recv_chars_[i][0]), bsize, MPI_CHAR, recvProcs_[i],
              tag1, comm_, &mpiReqs[i]);
  }

  //now pack and send our buffers.
  for(size_t i=0; i<sendProcs_.size(); ++i) {
    int proc = sendProcs_[i];
    fei::FillableMat* remoteMat = remotelyOwned_[proc];
    char* buffer = &(send_chars_[i][0]);
    fei::impl_utils::pack_FillableMat(*remoteMat, buffer);
    remoteMat->setValues(0.0);

    MPI_Send(&(send_chars_[i][0]), send_chars_[i].size(), MPI_CHAR, proc, tag1, comm_);
  }

  for(size_t ir=0; ir<recvProcs_.size(); ++ir) {
    int index;
    MPI_Status status;
    MPI_Waitany(numRecvProcs, &mpiReqs[0], &index, &status);
  }

  //and finally, unpack and store the received buffers.
  CSRMat recvMat;
  for(size_t ir=0; ir<recvProcs_.size(); ++ir) {
    size_t len = recv_chars_[ir].size();
    const char* data = &(recv_chars_[ir][0]);
    bool all_zeros = fei::impl_utils::unpack_CSRMat(data, data+len, recvMat);

    if (all_zeros) continue;

    int nrows = recvMat.getNumRows();

    for(int i=0; i<nrows; ++i) {
      int row = recvMat.getGraph().rowNumbers[i];
      int rowOffset = recvMat.getGraph().rowOffsets[i];
      int rowLen = recvMat.getGraph().rowOffsets[i+1] - rowOffset;

      int* indices = &(recvMat.getGraph().packedColumnIndices[rowOffset]);
      double* vals = &(recvMat.getPackedCoefs()[rowOffset]);
      int err = 0;
      if (haveBlockMatrix()) {
        err = giveToBlockMatrix(1, &row, rowLen, indices, &vals, accumulate);
      }
      else {
        err = giveToUnderlyingMatrix(1, &row, rowLen, indices,
                                     &vals, accumulate, FEI_DENSE_ROW);
      }
      if (err != 0) {
        FEI_COUT << "fei::Matrix_core::gatherFromOverlap ERROR storing "
         << "values for row=" << row<<", recv'd from proc " << recvProcs_[i]
         << FEI_ENDL;
        
        return(err);
      }
    }
  }
#endif

  return(0);
}

void
fei::Matrix_core::copyTransposeToWorkArrays(int numRows, int numCols,
                                            const double*const* values,
                                            std::vector<double>& work_1D,
                                            std::vector<const double*>& work_2D)
{
  //copy the transpose of 'values' into work-arrays.

  int arrayLen = numCols*numRows;
  if ((int)work_1D.size() != arrayLen) {
    work_1D.resize(arrayLen);
  }
  if ((int)work_2D.size() != numCols) {
    work_2D.resize(numCols);
  }

  const double** dataPtr = &work_2D[0];
  double* data1DPtr = &work_1D[0];

  for(int i=0; i<numCols; ++i) {
    dataPtr[i] = data1DPtr;

    for(int j=0; j<numRows; ++j) {
      data1DPtr[j] = values[j][i];
    }

    data1DPtr += numRows;
  }
}


int fei::Matrix_core::convertPtToBlk(int numRows,
				       const int* rows,
				       int numCols,
				       const int* cols,
				       int* blkRows,
				       int* blkRowOffsets,
				       int* blkCols,
				       int* blkColOffsets)
{
  snl_fei::PointBlockMap* pointBlockMap = vecSpace_->getPointBlockMap();

  int i;
  for(i=0; i<numRows; ++i) {
    if (pointBlockMap->getPtEqnInfo(rows[i], blkRows[i], blkRowOffsets[i])!=0){
      ERReturn(-1);
    }
  }
  for(i=0; i<numCols; ++i) {
    if(pointBlockMap->getPtEqnInfo(cols[i], blkCols[i], blkColOffsets[i])!=0){
      ERReturn(-1);
    }
  }

  return(0);
}

int fei::Matrix_core::copyPointRowsToBlockRow(int numPtRows,
						int numPtCols,
						const double*const* ptValues,
						int numBlkCols,
						const int* blkColDims,
						double** blkValues)
{
  int ptColOffset = 0;
  for(int blki=0; blki<numBlkCols; ++blki) {
    int blkvalueOffset = 0;
    double* blkvalues_i = blkValues[blki];
    for(int j=0; j<blkColDims[blki]; ++j) {
      int loc = ptColOffset+j;
      for(int i=0; i<numPtRows; ++i) {
        blkvalues_i[blkvalueOffset++] = ptValues[i][loc];
      }
    }
    ptColOffset += blkColDims[blki];
  }

  return(0);
}

void fei::Matrix_core::setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
{
  matrixGraph_ = matrixGraph;
  vecSpace_ = matrixGraph->getRowSpace();
}

