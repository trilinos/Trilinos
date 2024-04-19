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

#include <stk_transfer/TransferUtil.hpp>
#include <stk_transfer/ReducedDependencyCommData.hpp>
#include <stk_util/parallel/CouplingVersions.hpp>

namespace stk {
namespace transfer {
namespace impl {

int get_next_transfer_id()
{
  static int s_transferId = 1;
  return s_transferId++;
}

void exchange_transfer_ids(ReducedDependencyCommData& comm_data)
{
  std::vector<MPI_Request> receiveRequests(comm_data.numToMeshCommunications);
  std::vector<MPI_Request> sendRequests(comm_data.numFromMeshCommunications);

  int sendTag = 0;
  if (stk::util::get_common_coupling_version() >= 10) {
    comm_data.m_transferId = impl::get_next_transfer_id();
    sendTag = comm_data.m_transferId;
  }

  for (int ii = 0; ii < comm_data.numToMeshCommunications; ++ii) {
    int source = comm_data.uniqueToProcVec[ii];
    int recvMessageSize = 0;

    MPI_Irecv(nullptr, recvMessageSize, MPI_BYTE, source,
              MPI_ANY_TAG, comm_data.m_shared_comm, &receiveRequests[ii]);
  }

  for (int ii = 0; ii < comm_data.numFromMeshCommunications; ++ii) {
    int destination = comm_data.uniqueFromProcVec[ii];
    int sendMessageSize = 0;

    MPI_Isend(nullptr, sendMessageSize, MPI_BYTE, destination,
              sendTag, comm_data.m_shared_comm, &sendRequests[ii]);
  }

  std::vector<MPI_Status> receiveStati(receiveRequests.size());
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), receiveStati.data());

  if (stk::util::get_common_coupling_version() >= 10) {
    comm_data.m_otherTransferId.resize(comm_data.numToMeshCommunications);
    for(unsigned i=0; i<receiveStati.size(); ++i) {
      comm_data.m_otherTransferId[i] = receiveStati[i].MPI_TAG;
    }
  }

  MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
}

}
}
}

