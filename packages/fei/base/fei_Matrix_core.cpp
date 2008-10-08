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
#include <fei_Matrix_core.hpp>
#include <feiArray.hpp>
#include <fei_SSVec.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_CommUtils.hpp>

#include <fei_VectorSpace.hpp>
#include <snl_fei_PointBlockMap.hpp>
#include <fei_MatrixGraph.hpp>
#include <snl_fei_Utils.hpp>

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
    haveBlockMatrix_(false),
    haveFEMatrix_(false),
    globalOffsets_(),
    firstLocalOffset_(0),
    lastLocalOffset_(0)
{
  if (matrixGraph.get() == NULL) {
    throw std::runtime_error("fei::Matrix_core constructed with NULL fei::MatrixGraph");
  }

  eqnComm_.reset(new fei::EqnComm(comm_, numLocalEqns));

  localProc_ = fei::localProc(comm_);
  numProcs_ = fei::numProcs(comm_);

  vecSpace_ = matrixGraph->getRowSpace();

  remotelyOwned_.resize(numProcs_);
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    remotelyOwned_[i] = new SSMat;
  }

  setName("dbg");

  globalOffsets_ = eqnComm_->getGlobalOffsets();

  firstLocalOffset_ = globalOffsets_[localProc_];
  lastLocalOffset_ = globalOffsets_[localProc_+1]-1;
}

void
fei::Matrix_core::setEqnComm(fei::SharedPtr<fei::EqnComm> eqnComm)
{
  eqnComm_ = eqnComm;
  globalOffsets_ = eqnComm_->getGlobalOffsets();
  firstLocalOffset_ = globalOffsets_[localProc_];
  lastLocalOffset_ = globalOffsets_[localProc_+1]-1;
}

fei::Matrix_core::~Matrix_core()
{
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    delete remotelyOwned_[i];
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
}

void fei::Matrix_core::setName(const char* name)
{
  if (name == NULL) return;

  if (name_ == name) return;

  name_ = name;
}

int fei::Matrix_core::getOwnerProc(int globalEqn)
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

int fei::Matrix_core::gatherFromOverlap(bool accumulate)
{
  if (numProcs() == 1) return(0);

#ifndef FEI_SER
  //this function gathers shared-but-not-owned data onto the
  //owning processors.

  //iterate the remotelyOwned_ array and create a list of processors
  //that we will be sending data to. (processors which own matrix rows
  //that we share.)
  std::vector<int> sendProcs;
  for(unsigned i=0; i<remotelyOwned_.size(); ++i) {
    if ((int)i == localProc_) continue;
    if (remotelyOwned_[i] != NULL) {
      if (remotelyOwned_[i]->numNonzeros() == 0) {
        continue;
      }

      sendProcs.push_back((int)i);
    }
  }

  //now we can find out which procs we'll be receiving from.
  std::vector<int> recvProcs;
  fei::mirrorProcs(comm_, sendProcs, recvProcs);

  //next we'll declare arrays to receive into.
  std::vector<std::vector<int> > recv_ints(recvProcs.size());
  std::vector<std::vector<double> > recv_doubles(recvProcs.size());

  //...and an array for the sizes of the recv buffers:
  std::vector<int> recv_sizes(recvProcs.size()*2);
  std::vector<MPI_Request> mpiReqs(recvProcs.size()*2);
  std::vector<MPI_Status> mpiStatuses(recvProcs.size()*2);

  int tag1 = 20070430;
  int tag2 = 20070460;

  unsigned offset = 0;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    MPI_Irecv(&recv_sizes[offset], 1, MPI_INT, recvProcs[i],
              tag1, comm_, &mpiReqs[offset]);
    ++offset;
    MPI_Irecv(&recv_sizes[offset], 1, MPI_INT, recvProcs[i],
              tag2, comm_, &mpiReqs[offset]);
    ++offset;
  }

  //now we'll pack our to-be-sent data into buffers, and send the
  //sizes to the receiving procs:
  std::vector<std::vector<int> > send_ints(sendProcs.size());
  std::vector<std::vector<double> > send_doubles(sendProcs.size());

  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];

    snl_fei::packSSMat(*(remotelyOwned_[proc]), send_ints[i],
                        send_doubles[i]);

    int isize = send_ints[i].size();
    int dsize = send_doubles[i].size();

    MPI_Send(&isize, 1, MPI_INT, proc, tag1, comm_);
    MPI_Send(&dsize, 1, MPI_INT, proc, tag2, comm_);

    remotelyOwned_[proc]->logicalClear();
  }

  MPI_Waitall(mpiReqs.size(), &mpiReqs[0], &mpiStatuses[0]);


  //now resize our recv buffers, and post the recvs.
  offset = 0;
  unsigned offset2 = 0;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    int intsize = recv_sizes[offset++];
    int doublesize = recv_sizes[offset++];

    recv_ints[i].resize(intsize);
    recv_doubles[i].resize(doublesize);

    MPI_Irecv(&(recv_ints[i][0]), intsize, MPI_INT, recvProcs[i],
              tag1, comm_, &mpiReqs[offset2++]);
    MPI_Irecv(&(recv_doubles[i][0]), doublesize, MPI_DOUBLE, recvProcs[i],
              tag2, comm_, &mpiReqs[offset2++]);
  }

  //now send our packed buffers.
  for(unsigned i=0; i<sendProcs.size(); ++i) {
    int proc = sendProcs[i];

    MPI_Send(&(send_ints[i][0]), send_ints[i].size(), MPI_INT, proc, tag1, comm_);
    MPI_Send(&(send_doubles[i][0]), send_doubles[i].size(), MPI_DOUBLE,
             proc, tag2, comm_);
  }

  MPI_Waitall(mpiReqs.size(), &mpiReqs[0], &mpiStatuses[0]);


  //and finally, unpack and store the received buffers.
  SSMat recvMat;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    snl_fei::unpackIntoSSMat(recv_ints[i], recv_doubles[i], recvMat);

    feiArray<int>& rowNumbers = recvMat.getRowNumbers();
    feiArray<SSVec*>& rows = recvMat.getRows();
    for(int r=0; r<rows.length(); ++r) {
      int row = rowNumbers[r];
      int rowLen = rows[r]->length();
      int* indices = rows[r]->indices().dataPtr();
      double* vals = rows[r]->coefs().dataPtr();
      int err = giveToUnderlyingMatrix(1, &row, rowLen, indices,
                                       &vals, accumulate, FEI_DENSE_ROW);
      if (err != 0) {
        FEI_COUT << "fei::Matrix_Impl::gatherFromOverlap ERROR storing "
         << "values for row=" << row<<", recv'd from proc " << recvProcs[i]
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

