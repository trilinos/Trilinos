
#include <fei_CommUtils.hpp>
#include <fei_TemplateUtils.hpp>

#undef fei_file
#define fei_file "fei_CommUtils.cpp"

#include <fei_ErrMacros.hpp>

namespace fei {

//------------------------------------------------------------------------
int localProc(MPI_Comm comm)
{
#ifdef FEI_SER
  return 0;
#else
  int local_proc = 0;
  MPI_Comm_rank(comm, &local_proc);
  return local_proc;
#endif
}

//------------------------------------------------------------------------
int numProcs(MPI_Comm comm)
{
#ifdef FEI_SER
  return 1;
#else
  int num_procs = 1;
  MPI_Comm_size(comm, &num_procs);
  return num_procs;
#endif
}

//------------------------------------------------------------------------
void Barrier(MPI_Comm comm)
{
#ifndef FEI_SER
  MPI_Barrier(comm);
#endif
}

//------------------------------------------------------------------------
int mirrorProcs(MPI_Comm comm, std::vector<int>& toProcs, std::vector<int>& fromProcs)
{
  fromProcs.resize(0);
#ifdef FEI_SER
  fromProcs.push_back(0);
  return(0);
#else
  int num_procs = fei::numProcs(comm);
  std::vector<int> tmpIntData(num_procs*3, 0);

  int* buf = &tmpIntData[0];
  int* recvbuf = buf+num_procs;

  for(unsigned i=0; i<toProcs.size(); ++i) {
    buf[toProcs[i]] = 1;
  }

  for(int ii=2*num_procs; ii<3*num_procs; ++ii) {
    buf[ii] = 1;
  }

  CHK_MPI( MPI_Reduce_scatter(buf, &(buf[num_procs]), &(buf[2*num_procs]),
                              MPI_INT, MPI_SUM, comm) );

  int numRecvProcs = buf[num_procs];

  int tag = 11116;
  std::vector<MPI_Request> mpiReqs(numRecvProcs);

  int offset = 0;
  for(int ii=0; ii<numRecvProcs; ++ii) {
    CHK_MPI( MPI_Irecv(&(recvbuf[ii]), 1, MPI_INT, MPI_ANY_SOURCE, tag,
                       comm, &(mpiReqs[offset++])) );
  }

  for(unsigned i=0; i<toProcs.size(); ++i) {
    CHK_MPI( MPI_Send(&(toProcs[i]), 1, MPI_INT, toProcs[i], tag, comm) );
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

//------------------------------------------------------------------------
int mirrorCommPattern(MPI_Comm comm, comm_map* inPattern, comm_map*& outPattern)
{
#ifdef FEI_SER
  (void)inPattern;
  (void)outPattern;
#else
  int localP = localProc(comm);
  int numP  = numProcs(comm);

  if (numP < 2) return(0);

  std::vector<int> buf(numP*2, 0);

  int numInProcs = inPattern->getMap().size();
  std::vector<int> inProcs(numInProcs);
  fei::copyKeysToVector(inPattern->getMap(), inProcs);

  std::vector<int> outProcs;

  int err = mirrorProcs(comm, inProcs, outProcs);
  if (err != 0) ERReturn(-1);

  std::vector<int> recvbuf(outProcs.size(), 0);

  outPattern = new comm_map(0,1);

  MPI_Datatype mpi_ttype = fei::mpiTraits<int>::mpi_type();

  //now recv a length (the contents of buf[i]) from each "out-proc", which
  //will be the length of the equation data that will also be recvd from that
  //proc.
  std::vector<MPI_Request> mpiReqs(outProcs.size());
  std::vector<MPI_Status> mpiStss(outProcs.size());
  MPI_Request* requests = &mpiReqs[0];
  MPI_Status* statuses = &mpiStss[0];

  int firsttag = 11117;
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

  int tag2 = 11118;
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
    fei::copySetToArray(*(in_iter->second), sendList.size(), &sendList[0]);

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


//------------------------------------------------------------------------
int exchangeIntData(MPI_Comm comm,
                    const std::vector<int>& sendProcs,
                    std::vector<int>& sendData,
                    const std::vector<int>& recvProcs,
                    std::vector<int>& recvData)
{
  if (sendProcs.size() == 0 && recvProcs.size() == 0) return(0);
  if (sendProcs.size() != sendData.size()) return(-1);
#ifndef FEI_SER
  recvData.resize(recvProcs.size());
  std::vector<MPI_Request> mpiReqs;
  mpiReqs.resize(recvProcs.size());

  int tag = 11114;
  MPI_Datatype mpi_dtype = MPI_INT;

  //launch Irecv's for recvData:

  int localProc = fei::localProc(comm);
  int numRecvProcs = recvProcs.size();
  int req_offset = 0;
  for(unsigned i=0; i<recvProcs.size(); ++i) {
    if (recvProcs[i] == localProc) {--numRecvProcs; continue; }

    CHK_MPI( MPI_Irecv(&(recvData[i]), 1, mpi_dtype, recvProcs[i], tag,
                       comm, &mpiReqs[req_offset++]) );
  }

  //send the sendData:

  for(unsigned i=0; i<sendProcs.size(); ++i) {
    if (sendProcs[i] == localProc) continue;

    CHK_MPI( MPI_Send(&(sendData[i]), 1, mpi_dtype,
                      sendProcs[i], tag, comm) );
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

//------------------------------------------------------------------------
int Allreduce(MPI_Comm comm, bool localBool, bool& globalBool)
{
#ifndef FEI_SER
  int localInt = localBool ? 1 : 0;
  int globalInt = 0;

  CHK_MPI( MPI_Allreduce(&localInt, &globalInt, 1, MPI_INT, MPI_MAX, comm) );

  globalBool = globalInt==1 ? true : false;
#else
  globalBool = localBool;
#endif

  return(0);
}

}//namespace fei

