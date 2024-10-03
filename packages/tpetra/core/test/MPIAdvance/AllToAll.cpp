// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <random>
#include <vector>

#include <mpi.h>

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Kokkos_Core.hpp"

#include "Tpetra_Core.hpp"

#include <mpi_advance.h>

#include "print_vec.hpp"
#include "fake_alltoall.hpp"

/*! \brief all ranks send and receive nothing

    Passes if it does not crash or error out
*/
template <typename Device>
void test_nothing(MPI_Comm comm, bool nullBufs, bool sameBufs,
                  Teuchos::FancyOStream &out, bool &success) {
  static_assert(Kokkos::is_device_v<Device>, "");

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
template <typename Device>
void test_random(MPI_Comm comm, int seed, Teuchos::FancyOStream &out,
                 bool &success) {
  static_assert(Kokkos::is_device_v<Device>, "");

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
  Kokkos::View<char *, typename Device::memory_space> 
    sbuf("sbuf", sendcount * size),
    exp("exp", recvcount * size),
    act("act", recvcount * size);

  // fill send buf
  Kokkos::parallel_for(sbuf.size(), KOKKOS_LAMBDA (size_t i) {sbuf(i) = i;});

  // Use reference and MPI_Advance implementation to fill buffers
  Fake_Alltoall(sbuf.data(), sendcount, MPI_BYTE,
                exp.data(), recvcount, MPI_BYTE, comm);


  MPIX_Alltoall(sbuf.data(), sendcount, MPI_BYTE,
                act.data(), recvcount, MPI_BYTE, mpixComm);

  MPIX_Comm_free(mpixComm);

  // two recv buffers should be the same
  auto exp_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), exp);
  auto act_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), act);
  TEST_ASSERT(exp_h.size() == act_h.size());
  for (size_t i = 0; i < exp_h.size(); ++i) {
    TEST_ASSERT(exp_h(i) == act_h(i));
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
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing<device_type>(comm, false, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_nothing_null) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing<device_type>(comm, true, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_nothing_same) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing<device_type>(comm, false, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_nothing_nullsame) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing<device_type>(comm, true, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAll_random) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_random<device_type>(comm, 42, out, success);
}

// Let Tpetra initialize Kokkos
// We define this because we don't also include ${TEUCHOS_STD_UNIT_TEST_MAIN}
// in the CMakeLists.txt
int main(int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}