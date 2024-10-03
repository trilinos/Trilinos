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
#include "fake_alltoallv.hpp"

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

  std::vector<int> sendcounts(size, 0);
  std::vector<int> recvcounts(size, 0);

  std::vector<int> senddispls(size, 0);
  std::vector<int> recvdispls(size, 0);

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
  MPIX_Dist_graph_create_adjacent(
      comm, 0, /*indegree*/
      nullptr, /*sources*/
      nullptr, /*sourceweights*/
      0,       /*outdegree*/
      nullptr /*destinations*/, nullptr /*destweights*/, MPI_INFO_NULL /*info*/,
      0 /*reorder*/, &mpixComm);

  // reference implementation should be okay
  Fake_Alltoallv(sbuf, sendcounts.data(), senddispls.data(), MPI_BYTE, rbuf,
                   recvcounts.data(), recvdispls.data(), MPI_BYTE, comm);

  // MPI advance implementation
  MPIX_Neighbor_part_locality_alltoallv(sbuf, sendcounts.data(), senddispls.data(), MPI_BYTE,
                                          rbuf, recvcounts.data(), recvdispls.data(), MPI_BYTE,
                                          mpixComm);

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

  // seeded rng so all ranks generate the same communication plan
  // plan [i * size + j] is count i sends to j
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> dist(1, 100);
  std::vector<int> plan;
  for (int i = 0; i < size * size; ++i) {
    plan.push_back(dist(rng));
  }

  // read my part of the plan
  std::vector<int> sendcounts, recvcounts, senddispls, recvdispls; // alltoallv
  std::vector<int> nbrsendcounts, nbrrecvcounts, nbrsenddispls, nbrrecvdispls; // neighbor part locality alltoallv
  std::vector<int> sources, sourceweights, destinations, destweights; // communicator

  std::uniform_int_distribution<int> soffsetdist(0, 150);
  rng.seed(seed + rank); // different seed -> different displs per rank
  int initsdispl = soffsetdist(rng);
  int sdispl = initsdispl;
  int nbrsdispl = initsdispl;
  for (int dest = 0; dest < size; ++dest) {
    senddispls.push_back(sdispl);
    int count = plan[rank * size + dest];
    sendcounts.push_back(count);
    int offset = soffsetdist(rng);
    sdispl += count + offset;

    if (count > 0) {
      destinations.push_back(dest);
      destweights.push_back(count);
      nbrsendcounts.push_back(count);
      nbrsenddispls.push_back(nbrsdispl);
      nbrsdispl += count + offset;
    }
  }

  std::uniform_int_distribution<int> roffsetdist(0, 150);
  rng.seed(seed + size + rank); // different seed -> different displs (also per rank)
  int initrdispl = roffsetdist(rng);
  int rdispl = initrdispl;
  int nbrrdispl = initrdispl;  
  for (int source = 0; source < size; ++source) {
    recvdispls.push_back(rdispl);
    int count = plan[source * size + rank];
    recvcounts.push_back(count);
    int offset = roffsetdist(rng);
    rdispl += count + offset;

    if (count > 0) {
      sources.push_back(source);
      sourceweights.push_back(count);
      nbrrecvcounts.push_back(count);
      nbrrecvdispls.push_back(nbrrdispl);
      nbrrdispl += count + offset;
    }
  }


  print_vec(rank, "sources", sources);
  print_vec(rank, "sourceweights", sourceweights);
  print_vec(rank, "destinations", destinations);
  print_vec(rank, "destweights", destweights);

  // create MPIX communicator
  MPIX_Comm *mpixComm = nullptr;
  MPIX_Dist_graph_create_adjacent(
      comm, sources.size(), /*indegree*/
      sources.data(),       /*sources*/
      sourceweights.data(), /*sourceweights*/
      destinations.size(),  /*outdegree*/
      destinations.data() /*destinations*/, destweights.data() /*destweights*/,
      MPI_INFO_NULL /*info*/, 0 /*reorder*/, &mpixComm);

  // allocate send/recv bufs
  // displs are in elements, so the displs are correct since MPI_BYTE 
  // matches type in bufs, alltoallv calls as calculated above
  Kokkos::View<char *, typename Device::memory_space>
    sbuf("sbuf", sdispl), exp("exp", rdispl), act("act", nbrrdispl);

  // fill send buf
  Kokkos::parallel_for(sbuf.size(), KOKKOS_LAMBDA (size_t i) {sbuf(i) = i;});

  // Use reference and MPI_Advance implementation to fill buffers
  Fake_Alltoallv(sbuf.data(), sendcounts.data(), senddispls.data(), MPI_BYTE,
                   exp.data(), recvcounts.data(), recvdispls.data(), MPI_BYTE,
                   comm);

  print_vec(rank, "nbrsendcounts", nbrsendcounts);
  print_vec(rank, "nbrsenddispls", nbrsenddispls);
  print_vec(rank, "nbrrecvcounts", nbrrecvcounts);
  print_vec(rank, "nbrrecvdispls", nbrrecvdispls);


  MPIX_Neighbor_part_locality_alltoallv(sbuf.data(), nbrsendcounts.data(), nbrsenddispls.data(), MPI_BYTE,
                                          act.data(), nbrrecvcounts.data(), nbrrecvdispls.data(), MPI_BYTE,
                                          mpixComm);

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

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborPartLocalityAllToAllV_nothing) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing<device_type>(comm, false, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborPartLocalityAllToAllV_nothing_null) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing<device_type>(comm, true, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborPartLocalityAllToAllV_nothing_same) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing<device_type>(comm, false, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborPartLocalityAllToAllV_nothing_nullsame) {
  using execution_space = Kokkos::DefaultExecutionSpace;
  using memory_space = execution_space::memory_space;
  using device_type = Kokkos::Device<execution_space, memory_space>;
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing<device_type>(comm, true, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, NeighborPartLocalityAllToAllV_random) {
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