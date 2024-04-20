// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef stk_util_parallel_ParallelComm_hpp
#define stk_util_parallel_ParallelComm_hpp

#include "stk_util/stk_config.h"            // for STK_HAS_MPI
#include "stk_util/parallel/Parallel.hpp"   // for MPI_Irecv, MPI_Wait, MPI_Barrier, MPI_Send
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_util/parallel/CouplingVersions.hpp"
#include "stk_util/parallel/MPITagManager.hpp"  // for MPITagManager
#include "stk_util/parallel/DataExchangeUnknownPatternBlocking.hpp"  // for DataExchangeUnknownPatternBlocking
#include "stk_util/parallel/CommBuffer.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg, ThrowRequire
#include <stddef.h>
#include <cstddef>                          // for size_t, ptrdiff_t
#include <map>                              // for map
#include <stdexcept>                        // for runtime_error
#include <string>                           // for string
#include <vector>                           // for vector
#include <sstream>

namespace stk {

class CommBroadcast {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  int             parallel_size() const { return m_size ; }
  int             parallel_rank() const { return m_rank ; }
  int             root_rank() const { return m_root_rank; }

  /** Obtain the message buffer for the root_rank processor */
  CommBuffer & send_buffer();

  /** Obtain the message buffer for the local processor */
  CommBuffer & recv_buffer();

  //----------------------------------------

  CommBroadcast( ParallelMachine , int root_rank );

  void communicate();

  bool allocate_buffer( const bool local_flag = false );

  ~CommBroadcast();

private:

  CommBroadcast();
  CommBroadcast( const CommBroadcast & );
  CommBroadcast & operator = ( const CommBroadcast & );

  ParallelMachine m_comm ;
  int             m_size ;
  int             m_rank ;
  int             m_root_rank ;
  CommBuffer      m_buffer ;
};

template<typename PACK_ALGORITHM>
bool pack_and_communicate(CommBroadcast & comm, const PACK_ALGORITHM & algorithm)
{
  algorithm();
  comm.allocate_buffer();
  algorithm();
  comm.communicate();

  return (comm.parallel_rank() == comm.root_rank() && comm.send_buffer().capacity() > 0)
      || (comm.parallel_rank() != comm.root_rank() && comm.recv_buffer().capacity() > 0);
}

std::vector<int> ComputeReceiveList(std::vector<int>& sendSizeArray, MPI_Comm &mpi_communicator);


//
//  Parallel_Data_Exchange: General object exchange template with unknown comm plan
//
template<typename T>
void parallel_data_exchange_t(std::vector< std::vector<T> > &sendLists,
                              std::vector< std::vector<T> > &recvLists,
                              MPI_Comm &mpiCommunicator) {
#ifdef STK_HAS_MPI
  stk::util::print_unsupported_version_warning(3, __LINE__, __FILE__);

  if (stk::util::get_common_coupling_version() >= 4) {
    DataExchangeUnknownPatternBlocking exchanger(mpiCommunicator, 11242);
    exchanger.execute(sendLists, recvLists);
  } else
  {
    //
    //  Determine the number of processors involved in this communication
    //
    auto msg_tag = get_mpi_tag_manager().get_tag(mpiCommunicator, 10242);
    int num_procs;
    MPI_Comm_size(mpiCommunicator, &num_procs);
    int my_proc;
    MPI_Comm_rank(mpiCommunicator, &my_proc);
    STK_ThrowRequire((unsigned int) num_procs == sendLists.size() && (unsigned int) num_procs == recvLists.size());
    int class_size = sizeof(T);
    //
    //  Determine number of items each other processor will send to the current processor
    //
    std::vector<int> global_number_to_send(num_procs);
    for(int iproc=0; iproc<num_procs; ++iproc) {
      global_number_to_send[iproc] = sendLists[iproc].size();
    }
    std::vector<int> numToRecvFrom = ComputeReceiveList(global_number_to_send, mpiCommunicator);
    //
    //  Send the actual messages as raw byte streams.
    //
    std::vector<MPI_Request> recv_handles(num_procs);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      recvLists[iproc].resize(numToRecvFrom[iproc]);
      if(recvLists[iproc].size() > 0) {
        char* recv_buffer = (char*)recvLists[iproc].data();
        int recv_size = recvLists[iproc].size()*class_size;
        MPI_Irecv(recv_buffer, recv_size, MPI_CHAR, iproc, msg_tag, mpiCommunicator, &recv_handles[iproc]);
      }
    }
    MPI_Barrier(mpiCommunicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(sendLists[iproc].size() > 0) {
        char* send_buffer = (char*)sendLists[iproc].data();
        int send_size = sendLists[iproc].size()*class_size;
        MPI_Send(send_buffer, send_size, MPI_CHAR,
                iproc, msg_tag, mpiCommunicator);
      }
    }
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(recvLists[iproc].size() > 0) {
        MPI_Wait( &recv_handles[iproc], MPI_STATUS_IGNORE );
      }
    }
  }
#endif
}

//
//  Generalized comm plans
//
//  This plan assumes the send and recv lists have identical sizes so no extra sizing communications are needed
//
template<typename T>
void parallel_data_exchange_sym_t(std::vector< std::vector<T> > &send_lists,
                                  std::vector< std::vector<T> > &recv_lists,
                                  const MPI_Comm &mpi_communicator )
{
  //
  //  Determine the number of processors involved in this communication
  //
#if defined( STK_HAS_MPI)
  auto msg_tag = get_mpi_tag_manager().get_tag(mpi_communicator, 10242);
  int num_procs = stk::parallel_machine_size(mpi_communicator);
  int class_size = sizeof(T);

  //
  //  Send the actual messages as raw byte streams.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  std::vector<MPI_Request> send_handles(num_procs);
  int numRecvs = 0;
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    recv_lists[iproc].resize(send_lists[iproc].size());
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)recv_lists[iproc].data();
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR,
                iproc, msg_tag, mpi_communicator, &recv_handles[numRecvs++]);
    }
  }
  if (stk::util::get_common_coupling_version() <= 11) {
    MPI_Barrier(mpi_communicator);
  }

  int numSends = 0;
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)send_lists[iproc].data();
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Isend(send_buffer, send_size, MPI_CHAR,
                iproc, msg_tag, mpi_communicator, &send_handles[numSends++]);
    }
  }

  MPI_Waitall(numSends, send_handles.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(numRecvs, recv_handles.data(), MPI_STATUSES_IGNORE);
#endif
}

//
//The checkInput flag must have the same value on every processor, or this function will hang.
//
template<typename T>
inline
void parallel_data_exchange_nonsym_known_sizes_t(const int* sendOffsets,
                                                 T* sendData,
                                                 const int* recvOffsets,
                                                 T* recvData,
                                                 MPI_Comm mpi_communicator,
                                                 bool checkInput = false)
{
#if defined( STK_HAS_MPI)
  const auto msg_tag = get_mpi_tag_manager().get_tag(mpi_communicator, 10243); //arbitrary tag value, anything less than 32768 is legal
  const int num_procs = stk::parallel_machine_size(mpi_communicator);
  const int bytesPerScalar = sizeof(T);

  if (checkInput && stk::util::get_common_coupling_version() > 11) {
    std::vector<int> reducedSendOffsets(num_procs+1);
    const int result = MPI_Allreduce(sendOffsets, reducedSendOffsets.data(), num_procs+1, MPI_INT, MPI_SUM, mpi_communicator);
    int myProc = -1;
    MPI_Comm_rank(mpi_communicator, &myProc);
    const int totalNumRecvs = reducedSendOffsets[myProc+1] - reducedSendOffsets[myProc];
    const int expectedNumRecvs = recvOffsets[num_procs];
    STK_ThrowRequireMsg(result == MPI_SUCCESS, "parallel_data_exchange_nonsym_known_sizes_t MPI_Allreduce return-code "<<result);
    const int localErr = totalNumRecvs != expectedNumRecvs;
    const int globalErr = stk::get_global_max(mpi_communicator, localErr);
    if (globalErr) {
      std::ostringstream oss;
      if (totalNumRecvs != expectedNumRecvs) {
        oss << "parallel_data_exchange_nonsym_known_sizes_t P"<<myProc<<" expecting to recv "<<expectedNumRecvs<<" values, but senders planning to send "<<totalNumRecvs<<" total values" << std::endl;
      }
      STK_ThrowErrorMsg(oss.str());
    }
  }

  //
  //  Send the actual messages as raw byte streams.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  std::vector<MPI_Request> send_handles;
  if (stk::util::get_common_coupling_version() <= 13) {
    send_handles.resize(num_procs);
  }
  int numRecvs = 0;
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    const int recvSize = recvOffsets[iproc+1]-recvOffsets[iproc];
    if(recvSize > 0) {
      char* recvBuffer = (char*)(&recvData[recvOffsets[iproc]]);
      const int recvSizeBytes = recvSize*bytesPerScalar;
      MPI_Irecv(recvBuffer, recvSizeBytes, MPI_CHAR, iproc, msg_tag, mpi_communicator, &recv_handles[numRecvs++]);
    }
  }

  if (stk::util::get_common_coupling_version() <= 11 || stk::util::get_common_coupling_version() >= 13) {
    MPI_Barrier(mpi_communicator);
  }

  int numSends = 0;
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    const int sendSize = sendOffsets[iproc+1]-sendOffsets[iproc];
    if(sendSize > 0) {
      char* sendBuffer = (char*)(&sendData[sendOffsets[iproc]]);
      const int sendSizeBytes = sendSize*bytesPerScalar;
      if (stk::util::get_common_coupling_version() <= 13) {
        MPI_Isend(sendBuffer, sendSizeBytes, MPI_CHAR, iproc, msg_tag, mpi_communicator, &send_handles[numSends++]);
      }
      else {
        MPI_Send(sendBuffer, sendSizeBytes, MPI_CHAR, iproc, msg_tag, mpi_communicator);
      }
    }
  }

  if (stk::util::get_common_coupling_version() <= 13) {
    MPI_Waitall(numSends, send_handles.data(), MPI_STATUSES_IGNORE);
  }
  MPI_Waitall(numRecvs, recv_handles.data(), MPI_STATUSES_IGNORE);
#endif
}

//
//  This plan assumes the send and recv lists are matched, but that the actual amount of data to send is unknown.
//  A processor knows which other processors it will be receiving data from, but does not know how much data.
//  Thus the comm plan is known from the inputs, but an additional message sizing call must be done.
//
template<typename T>
void parallel_data_exchange_sym_unknown_size_t(std::vector< std::vector<T> > &send_lists,
                                               std::vector< std::vector<T> > &recv_lists,
                                               MPI_Comm &mpi_communicator )
{
#if defined( STK_HAS_MPI)
  const auto msg_tag = get_mpi_tag_manager().get_tag(mpi_communicator, 10242);
  int num_procs = stk::parallel_machine_size(mpi_communicator);
  int class_size = sizeof(T);

  //
  //  Send the message sizes
  //
  std::vector<int> send_msg_sizes(num_procs);
  std::vector<int> recv_msg_sizes(num_procs);
  std::vector<MPI_Request> recv_handles(num_procs);

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    send_msg_sizes[iproc] = send_lists[iproc].size();
  }    
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size()>0) {
      MPI_Irecv(&recv_msg_sizes[iproc], 1, MPI_INT, iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }

  if (stk::util::get_common_coupling_version() <= 11) {
    MPI_Barrier(mpi_communicator);
  }

  int numSends = 0;
  std::vector<MPI_Request> send_handles(num_procs);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size()>0) {
      MPI_Isend(&send_msg_sizes[iproc], 1, MPI_INT, iproc, msg_tag, mpi_communicator, &send_handles[numSends++]);
    }
  }

  MPI_Waitall(numSends, send_handles.data(), MPI_STATUSES_IGNORE);

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Wait( &recv_handles[iproc], MPI_STATUS_IGNORE );
      recv_lists[iproc].resize(recv_msg_sizes[iproc]);
    }
  }
  //
  //  Send the actual messages as raw byte streams.
  //
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)recv_lists[iproc].data();
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR,
                iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }

  if (stk::util::get_common_coupling_version() <= 11) {
    MPI_Barrier(mpi_communicator);
  }

  numSends = 0;
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)send_lists[iproc].data();
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Isend(send_buffer, send_size, MPI_CHAR,
               iproc, msg_tag, mpi_communicator, &send_handles[numSends++]);
    }
  }

  MPI_Waitall(numSends, send_handles.data(), MPI_STATUSES_IGNORE);

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Wait( &recv_handles[iproc], MPI_STATUS_IGNORE );
    }
  }
#endif
}

template<typename T, typename MsgPacker, typename MsgUnpacker>
void parallel_data_exchange_sym_pack_unpack(MPI_Comm mpi_communicator,
                                            const std::vector<int>& comm_procs,
                                            MsgPacker& pack_msg,
                                            MsgUnpacker& unpack_msg,
                                            bool deterministic)
{
#if defined( STK_HAS_MPI)
  const auto msg_tag = get_mpi_tag_manager().get_tag(mpi_communicator, 10242);
  int class_size = sizeof(T);

  int num_comm_procs = comm_procs.size();
  std::vector<std::vector<T> > send_data(num_comm_procs);
  std::vector<std::vector<T> > recv_data(num_comm_procs);
  std::vector<MPI_Request> send_requests(num_comm_procs);
  std::vector<MPI_Request> recv_requests(num_comm_procs);

  for(int i=0; i<num_comm_procs; ++i) {
    int iproc = comm_procs[i];
    pack_msg(iproc, send_data[i]);
    recv_data[i].resize(send_data[i].size());

    char* recv_buffer = (char*)recv_data[i].data();
    int buf_size = recv_data[i].size()*class_size;
    MPI_Irecv(recv_buffer, buf_size, MPI_CHAR, iproc, msg_tag, mpi_communicator, &recv_requests[i]);
    char* send_buffer = (char*)send_data[i].data();
    MPI_Isend(send_buffer, buf_size, MPI_CHAR, iproc, msg_tag, mpi_communicator, &send_requests[i]);
  }

  for(int i = 0; i < num_comm_procs; ++i) {
      int idx = i;
      if (deterministic) {
          MPI_Wait(&recv_requests[i], MPI_STATUS_IGNORE);
      }   
      else {
          MPI_Waitany(num_comm_procs, recv_requests.data(), &idx, MPI_STATUS_IGNORE);
      }   
      unpack_msg(comm_procs[idx], recv_data[idx]);
  }

  MPI_Waitall(num_comm_procs, send_requests.data(), MPI_STATUSES_IGNORE);
#endif
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

