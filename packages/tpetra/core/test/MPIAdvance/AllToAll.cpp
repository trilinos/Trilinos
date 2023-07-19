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

#include <iostream>
#include <random>
#include <vector>

#include <mpi.h>

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_Core.hpp"

#include <mpi_advance.h>

#include "print_vec.hpp"
#include "fake_alltoall.hpp"

/*! \brief all ranks send and receive nothing

    Passes if it does not crash or error out
*/
void test_nothing(MPI_Comm comm, bool nullBufs, bool sameBufs,
                  Teuchos::FancyOStream &out, bool &success) {

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  int sendcount = 0;
  int recvcount = 0;

  char stackByte0, stackByte1;
  void *sbuf, *rbuf;
  if (nullBufs) {
    sbuf = nullptr;
    rbuf = nullptr;
  } else {
    sbuf = &stackByte0;
    rbuf = &stackByte1;
  }

  if (sameBufs) {
    sbuf = rbuf;
  }

  // create MPIX communicator
  MPIX_Comm *mpixComm = nullptr;
  MPIX_Comm_init(&mpixComm, comm);

  // reference implementation should be okay
  Fake_Alltoall(sbuf, sendcount, MPI_BYTE, rbuf,
                recvcount, MPI_BYTE, comm);

  // MPI advance implementation
  MPIX_Alltoall(sbuf, sendcount, MPI_BYTE, rbuf,
                recvcount, MPI_BYTE, mpixComm);

  MPIX_Comm_free(mpixComm);

  // we just require that we got this far
  success = true;
}

/*! \brief all ranks send and receive some
 */
void test_random(MPI_Comm comm, int seed, Teuchos::FancyOStream &out,
                 bool &success) {

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  int sendcount, recvcount; // alltoall

  // seeded rng so all ranks generate the same count
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> dist(1, 100);  
  sendcount = dist(rng);
  recvcount = sendcount;

  // create MPIX communicator
  MPIX_Comm *mpixComm = nullptr;
  MPIX_Comm_init(&mpixComm, comm);

  // allocate send/recv bufs
  std::vector<char> sbuf(sendcount * size), exp(recvcount * size), act(recvcount * size);

  // fill send buf
  std::iota(sbuf.begin(), sbuf.end(), 0); // 0, 1, 2, ...

  // Use reference and MPI_Advance implementation to fill buffers
  Fake_Alltoall(sbuf.data(), sendcount, MPI_BYTE,
                   exp.data(), recvcount, MPI_BYTE, comm);


  MPIX_Alltoall(sbuf.data(), sendcount, MPI_BYTE,
                 act.data(), recvcount, MPI_BYTE, mpixComm);

  MPIX_Comm_free(mpixComm);

  // two recv buffers should be the same
  for (int i = 0; i < recvcount * size; ++i) {
    TEST_ASSERT(exp[i] == act[i]);
  }
}

static MPI_Comm tpetra_default_comm_as_mpi_comm() {
  Teuchos::RCP<const Teuchos::Comm<int>> teuchosComm = Tpetra::getDefaultComm();

  const Teuchos::MpiComm<int> *teuchosMpiComm =
      dynamic_cast<const Teuchos::MpiComm<int> *>(teuchosComm.get());

  if (teuchosMpiComm) {
    return *(teuchosMpiComm->getRawMpiComm());
  } else {
    throw std::runtime_error(
        "couldn't convert tpetra::getDefaultComm() into MPI_Comm");
  }
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_nothing) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, false, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_nothing_null) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, true, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_nothing_same) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, false, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_nothing_nullsame) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, true, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_random) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_random(comm, 42, out, success);
}
