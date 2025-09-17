// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once

#include <mpi.h>

/*! reference alltoallv impl
 */
inline void Fake_Alltoallv(const void *sendbuf, const int *sendcounts,
                           const int *sdispls, MPI_Datatype sendtype, void *recvbuf,
                           const int *recvcounts, const int *rdispls,
                           MPI_Datatype recvtype, MPI_Comm comm) {
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
    MPI_Irecv(&rb[rdispls[source] * recvsize], recvcounts[source], recvtype, source,
              ARBITRARY_TAG, comm, &rreqs[source]);
  }
  for (int dest = 0; dest < size; ++dest) {
    MPI_Isend(&sb[sdispls[dest] * sendsize], sendcounts[dest], sendtype, dest,
              ARBITRARY_TAG, comm, &sreqs[dest]);
  }

  // wait for communication to finish
  MPI_Waitall(size, sreqs.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(size, rreqs.data(), MPI_STATUSES_IGNORE);
}