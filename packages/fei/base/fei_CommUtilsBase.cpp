/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_CommUtilsBase.hpp>
#include <fei_mpiTraits.hpp>
#include <fei_ctg_set.hpp>
#include <fei_chk_mpi.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_ErrMacros.hpp>

namespace fei {

CommUtilsBase::CommUtilsBase(MPI_Comm comm)
  : commCore_(fei::CommCore::create()),
    comm_(comm),
    numProcs_(1),
    localProc_(0)
{
  //Figure out how many procs there are, and what our rank is.
#ifndef FEI_SER
  MPI_Comm_rank(comm_, &localProc_);
  MPI_Comm_size(comm_, &numProcs_);
#endif
}

CommUtilsBase::~CommUtilsBase()
{
}

int CommUtilsBase::mirrorProcs(std::vector<int>& toProcs,
                               std::vector<int>& fromProcs)
{
  commCore_->increment_tag();
  fromProcs.resize(0);
#ifdef FEI_SER
  fromProcs.push_back(0);
  return(0);
#else
  std::vector<int>& tmpIntData = commCore_->tmpIntData_;
  tmpIntData.assign(numProcs_*3, 0);

  int* buf = &tmpIntData[0];
  int* recvbuf = buf+numProcs_;

  for(unsigned i=0; i<toProcs.size(); ++i) {
    buf[toProcs[i]] = 1;
  }

  for(int ii=2*numProcs_; ii<3*numProcs_; ++ii) {
    buf[ii] = 1;
  }

  CHK_MPI( MPI_Reduce_scatter(buf, &(buf[numProcs_]), &(buf[2*numProcs_]),
                              MPI_INT, MPI_SUM, comm_) );

  int numRecvProcs = buf[numProcs_];

  int tag = commCore_->get_tag();
  std::vector<MPI_Request>& mpiReqs = commCore_->mpiReqs_;
  mpiReqs.resize(numRecvProcs);

  int offset = 0;
  for(int ii=0; ii<numRecvProcs; ++ii) {
    CHK_MPI( MPI_Irecv(&(recvbuf[ii]), 1, MPI_INT, MPI_ANY_SOURCE, tag,
                       comm_, &(mpiReqs[offset++])) );
  }

  for(unsigned i=0; i<toProcs.size(); ++i) {
    CHK_MPI( MPI_Send(&(toProcs[i]), 1, MPI_INT, toProcs[i], tag, comm_) );
  }

  MPI_Status status;
  for(int ii=0; ii<numRecvProcs; ++ii) {
    int index;
    MPI_Waitany(numRecvProcs, &mpiReqs[0], &index, &status);
    fromProcs.push_back(status.MPI_SOURCE);
  }

  std::sort(fromProcs.begin(), fromProcs.end());

  return(0);
#endif
}

int CommUtilsBase::exchangeIntData(std::vector<int>& sendProcs,
                                   std::vector<int>& sendData,
                                   std::vector<int>& recvProcs,
                                   std::vector<int>& recvData)
{
  commCore_->increment_tag();

  if (sendProcs.size() == 0 && recvProcs.size() == 0) return(0);
  if (sendProcs.size() != sendData.size()) return(-1);
#ifndef FEI_SER
  recvData.resize(recvProcs.size());
  std::vector<MPI_Request>& mpiReqs = commCore_->mpiReqs_;
  mpiReqs.resize(recvProcs.size());

  int tag = commCore_->get_tag();
  MPI_Datatype mpi_dtype = MPI_INT;

  //launch Irecv's for recvData:

  int numRecvProcs = recvProcs.size();
  int req_offset = 0;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    if (recvProcs[i] == localProc_) {--numRecvProcs; continue; }

    CHK_MPI( MPI_Irecv(&(recvData[i]), 1, mpi_dtype, recvProcs[i], tag,
                       comm_, &mpiReqs[req_offset++]) );
  }

  //send the sendData:

  for(unsigned i=0; i<sendProcs.size(); ++i) {
    if (sendProcs[i] == localProc_) continue;

    CHK_MPI( MPI_Send(&(sendData[i]), 1, mpi_dtype,
                      sendProcs[i], tag, comm_) );
  }

  //complete the Irecv's:

  for(int ii=0; ii<numRecvProcs; ++ii) {
    int index;
    MPI_Status status;
    CHK_MPI( MPI_Waitany(numRecvProcs, &mpiReqs[0], &index, &status) );
  }

#endif
  return(0);
}

int CommUtilsBase::mirrorCommPattern(comm_map* inPattern,
                                     comm_map*& outPattern)
{
#ifdef FEI_SER
  (void)inPattern;
  (void)outPattern;
#else
  int localP = localProc();
  int numP  = numProcs();
  MPI_Comm comm = getCommunicator();

  if (numP < 2) return(0);

  std::vector<int> buf(numP*2, 0);

  int numInProcs = inPattern->getMap().size();
  std::vector<int> inProcs(numInProcs);
  fei::copyKeysToVector(inPattern->getMap(), inProcs);

  std::vector<int> outProcs;

  int err = mirrorProcs(inProcs, outProcs);
  if (err != 0) ERReturn(-1);

  std::vector<int> recvbuf(outProcs.size(), 0);

  outPattern = new comm_map(0,1);

  MPI_Datatype mpi_ttype = fei::mpiTraits<int>::mpi_type();

  //now recv a length (the contents of buf[i]) from each "out-proc", which
  //will be the length of the equation data that will also be recvd from that
  //proc.
  std::vector<MPI_Request>& mpiReqs = commCore_->mpiReqs_;
  std::vector<MPI_Status>& mpiStss = commCore_->mpiStatuses_;
  mpiReqs.resize(outProcs.size());
  mpiStss.resize(outProcs.size());
  MPI_Request* requests = &mpiReqs[0];
  MPI_Status* statuses = &mpiStss[0];

  commCore_->increment_tag();
  int firsttag = commCore_->get_tag();
  int offset = 0;
  int* outProcsPtr = &outProcs[0];
  for(unsigned i=0; i<outProcs.size(); ++i) {
    if (MPI_Irecv(&(recvbuf[i]), 1, MPI_INT, outProcsPtr[i], firsttag,
                  comm, &requests[offset++]) != MPI_SUCCESS) ERReturn(-1);
  }

  comm_map::map_type& in_row_map = inPattern->getMap();
  comm_map::map_type::iterator
    in_iter = in_row_map.begin(),
    in_end  = in_row_map.end();
 
  int* inProcsPtr = &inProcs[0];
  for(int ii=0; in_iter!= in_end; ++in_iter, ++ii) {
    comm_map::row_type* in_row = in_iter->second;
    buf[ii] = in_row->size();
    if (MPI_Send(&(buf[ii]), 1, MPI_INT, inProcsPtr[ii], firsttag,
                 comm) != MPI_SUCCESS) ERReturn(-1);
  }

  int numOutProcs = outProcs.size();

  MPI_Waitall(numOutProcs, requests, statuses);
  std::vector<int> lengths(numOutProcs);
  int totalRecvLen = 0;
  offset = 0;
  for(int ii=0; ii<numOutProcs; ++ii) {
    if (recvbuf[ii] > 0) {
      lengths[offset++] = recvbuf[ii];
      totalRecvLen += recvbuf[ii];
    }
  }

  //now we need to create the space into which we'll receive the
  //lists that other procs send to us.
  std::vector<int> recvData(totalRecvLen, 999999);

  commCore_->increment_tag();
  int tag2 = commCore_->get_tag();
  offset = 0;
  for(int ii=0; ii<numOutProcs; ++ii) {
    CHK_MPI(MPI_Irecv(&(recvData[offset]), lengths[ii], mpi_ttype,
                      outProcs[ii], tag2, comm, &requests[ii]) );
    offset += lengths[ii];
  }

  std::vector<int> sendList;

  in_iter = in_row_map.begin();

  for(int ii=0; in_iter != in_end; ++in_iter,++ii) {
    if (inProcs[ii] == localP) {
      continue;
    }
    sendList.resize(in_iter->second->size());
    err = in_iter->second->copy_to_array(sendList.size(), &sendList[0]);
    if (err != 0) return(err);

    CHK_MPI(MPI_Send(&sendList[0], sendList.size(), mpi_ttype,
                     inProcs[ii], tag2, comm) );
  }

  //our final communication operation is to catch the Irecvs we started above.
  for(int ii=0; ii<numOutProcs; ++ii) {
    MPI_Wait(&requests[ii], &statuses[ii]);
  }

  //now we've completed all the communication, so we're ready to put the data
  //we received into the outPattern object.
  offset = 0;
  for(int ii=0; ii<numOutProcs; ii++) {
    outPattern->addIndices(outProcs[ii], lengths[ii],
                           &(recvData[offset]));
    offset += lengths[ii];
  }

#endif
  return(0);
}


int CommUtilsBase::numProcs() const { return( numProcs_ ); }
int CommUtilsBase::localProc() const { return( localProc_); }
MPI_Comm CommUtilsBase::getCommunicator() const { return( comm_ ); }
int CommUtilsBase::Barrier() const
{
#ifdef FEI_SER
  return(0);
#else
  MPI_Comm tmpcomm = comm_;
  return( MPI_Barrier(tmpcomm) );
#endif
}

} //namespace fei

