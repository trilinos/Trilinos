#include <iostream>
#include <random>
#include <vector>

#include <mpi.h>

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_UnitTestHarness.hpp"

#include "Tpetra_Core.hpp"

#include <mpi_advance.h>

#if defined(HAVE_TPETRACORE_MPI_ADVANCE)
#warning HAVE_TPETRACORE_MPI_ADVANCE is defined!

/*! reference alltoallv impl
 */
void Tpetra_Alltoallv(const void *sendbuf, const int *sendcounts,
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

  // issue sends & recvs
  for (int source = 0; source < size; ++source) {
    MPI_Irecv(&rb[rdispls[source]], recvcounts[source], recvtype, source,
              ARBITRARY_TAG, comm, &rreqs[source]);
  }
  for (int dest = 0; dest < size; ++dest) {
    MPI_Isend(&sb[sdispls[dest]], sendcounts[dest], sendtype, dest,
              ARBITRARY_TAG, comm, &sreqs[dest]);
  }

  // wait for communication to finish
  MPI_Waitall(size, sreqs.data(), MPI_STATUSES_IGNORE);
  MPI_Waitall(size, rreqs.data(), MPI_STATUSES_IGNORE);
}

/*! \brief all ranks send and receive nothing

    Passes if it does not crash or error out
*/
void test_nothing(MPI_Comm comm, bool nullBufs, bool sameBufs,
                  Teuchos::FancyOStream &out, bool &success) {

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
  Tpetra_Alltoallv(sbuf, sendcounts.data(), senddispls.data(), MPI_BYTE, rbuf,
                   recvcounts.data(), recvdispls.data(), MPI_BYTE, comm);

  // MPI advance implementation
  MPIX_Neighbor_alltoallv(sbuf, sendcounts.data(), senddispls.data(), MPI_BYTE,
                          rbuf, recvcounts.data(), recvdispls.data(), MPI_BYTE,
                          mpixComm);

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

  // seeded rng so all ranks generate the same communication plan
  // plan [i * size + j] is count i sends to j
  std::mt19937 rng(seed);
  std::uniform_int_distribution<int> dist(1, 100);
  std::vector<int> plan;
  for (int i = 0; i < size * size; ++i) {
    plan.push_back(dist(rng));
  }

  // read my part of the plan
  std::vector<int> sendcounts, recvcounts, senddispls, recvdispls;

  int sdispl = 0;
  for (int dest = 0; dest < size; ++dest) {
    senddispls.push_back(sdispl);
    int count = plan[rank * size + dest];
    sendcounts.push_back(count);
    sdispl += count;
  }

  int rdispl = 0;
  for (int source = 0; source < size; ++source) {
    recvdispls.push_back(rdispl);
    int count = plan[source * size + rank];
    recvcounts.push_back(count);
    rdispl += count;
  }

  // allocate send/recv bufs
  std::vector<char> sbuf(sdispl), exp(rdispl), act(rdispl);

  // fill send buf
  std::iota(sbuf.begin(), sbuf.end(), 0); // 0, 1, 2, ...

  // Use reference and MPI_Advance implementation to fill buffers
  Tpetra_Alltoallv(sbuf.data(), sendcounts.data(), senddispls.data(), MPI_BYTE,
                   exp.data(), recvcounts.data(), recvdispls.data(), MPI_BYTE,
                   comm);
  MPI_Alltoallv(sbuf.data(), sendcounts.data(), senddispls.data(), MPI_BYTE,
                act.data(), recvcounts.data(), recvdispls.data(), MPI_BYTE,
                comm);

  // two recv buffers should be the s ame
  for (int i = 0; i < rdispl; ++i) {
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

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAllV_nothing) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, false, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAllV_nothing_null) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, true, false, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAllV_nothing_same) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, false, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAllV_nothing_nullsame) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_nothing(comm, true, true, out, success);
}

TEUCHOS_UNIT_TEST(MpiAdvance, AllToAllV_random) {
  MPI_Comm comm = tpetra_default_comm_as_mpi_comm();
  test_random(comm, 42, out, success);
}

#endif // HAVE_TPETRACORE_MPI_ADVANCE