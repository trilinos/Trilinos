/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER
*/

#pragma once

#include <mpi.h>

/*! reference alltoall impl
 */
inline void Fake_Alltoall(const void *sendbuf, int sendcount,
                      MPI_Datatype sendtype, void *recvbuf,
                      int recvcount, MPI_Datatype recvtype, 
                      MPI_Comm comm) {

  constexpr int ARBITRARY_TAG = 0;

  // communicator properties
  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  // requests for sends and receives
  std::vector<MPI_Request> sreqs(size, MPI_REQUEST_NULL);
  std::vector<MPI_Request> rreqs(size, MPI_REQUEST_NULL);

  auto rb = reinterpret_cast<char *>(recvbuf);
  auto sb = reinterpret_cast<const char *>(sendbuf);
  
  // get sizes for indexing into sb / rb
  int sendsize, recvsize;
  MPI_Type_size(sendtype, &sendsize);
  MPI_Type_size(recvtype, &recvsize);

  // issue sends & recvs
  for (int source = 0; source < size; ++source) {
    MPI_Irecv(&rb[source * recvcount * recvsize], recvcount, recvtype, source,
              ARBITRARY_TAG, comm, &rreqs[source]);
  }
  for (int dest = 0; dest < size; ++dest) {
    MPI_Isend(&sb[dest * sendcount * sendsize], sendcount, sendtype, dest,
              ARBITRARY_TAG, comm, &sreqs[dest]);
  }

  // wait for communication to finish
  MPI_Waitall(size, sreqs.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(size, rreqs.data(), MPI_STATUSES_IGNORE);
}