// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2026 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_Ialltofewv.hpp"
#include "Tpetra_TestingUtilities.hpp"

#include <Teuchos_UnitTestHarness.hpp>

#include <array>
#include <string_view>
#include <vector>

namespace {

using Tpetra::Details::Ialltofewv;

int cachedAllocations = 0;

void countCachedAllocation(Kokkos::Tools::SpaceHandle, const char *label, const void *, uint64_t) {
  const std::string_view name(label);
  if (name == "rootBufDev" || name == "rootBufHost" ||
      name == "aggBufDev" || name == "argsDev" || name == "argsHost") {
    ++cachedAllocations;
  }
}

struct CachedAllocationCounter {
  CachedAllocationCounter()
    : previous_(Kokkos::Tools::Experimental::get_callbacks().allocate_data) {
    cachedAllocations = 0;
    Kokkos::Tools::Experimental::set_allocate_data_callback(countCachedAllocation);
  }

  ~CachedAllocationCounter() {
    Kokkos::Tools::Experimental::set_allocate_data_callback(previous_);
  }

  int globalCount() const {
    int result = 0;
    MPI_Allreduce(&cachedAllocations, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return result;
  }

 private:
  Kokkos::Tools::allocateDataFunction previous_;
};

struct Fixture {
  explicit Fixture(int value, int count = 1)
    : send(count, value)
    , sendcounts{count} {
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (rank == root) {
      recv.resize(size * count, -1);
      recvcounts.resize(size, count);
      rdispls.resize(size);
      for (int i = 0; i < size; ++i) rdispls[i] = i * count;
    }
  }

  void post(Ialltofewv &impl, int aggTag, int rootTag) {
    const int err = impl.post<false>(send.data(), sendcounts.data(), sdispls.data(), MPI_INT,
                                     recv.data(), recvcounts.data(), rdispls.data(),
                                     roots.data(), int(roots.size()), MPI_INT, aggTag, rootTag,
                                     MPI_COMM_WORLD, &req);
    TEUCHOS_ASSERT_EQUALITY(err, MPI_SUCCESS);
  }

  int rank = 0;
  int size = 0;
  std::vector<int> send;
  int root = 0;
  std::array<int, 1> sendcounts{1};
  std::array<int, 1> sdispls{0};
  std::array<int, 1> roots{root};
  std::vector<int> recv;
  std::vector<int> recvcounts;
  std::vector<int> rdispls;
  Ialltofewv::Req req;
};

void checkResult(const Fixture &fixture, int offset,
                 Teuchos::FancyOStream &out, bool &success) {
  if (fixture.rank != fixture.root) return;
  for (int rank = 0; rank < fixture.size; ++rank) {
    TEST_EQUALITY(fixture.recv[rank], offset + rank);
  }
}

TEUCHOS_UNIT_TEST(Ialltofewv, differentWaitOrder) {
  // Test that we don't deadlock when different ranks wait on different Ialltofewv
  // in different orders.

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Ialltofewv firstImpl;
  Ialltofewv secondImpl;
  Fixture first(rank);
  Fixture second(100 + rank);
  first.post(firstImpl, 101, 102);
  second.post(secondImpl, 103, 104);

  if (rank % 2 == 0) {
    TEST_EQUALITY(firstImpl.wait(first.req), MPI_SUCCESS);
    TEST_EQUALITY(secondImpl.wait(second.req), MPI_SUCCESS);
  } else {
    TEST_EQUALITY(secondImpl.wait(second.req), MPI_SUCCESS);
    TEST_EQUALITY(firstImpl.wait(first.req), MPI_SUCCESS);
  }

  checkResult(first, 0, out, success);
  checkResult(second, 100, out, success);
}

TEUCHOS_UNIT_TEST(Ialltofewv, pollingOneProgressesAll) {
  // Check that get_status() on a request makes progress on other Ialltofewvs

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Ialltofewv firstImpl;
  Ialltofewv secondImpl;
  Fixture first(rank);
  Fixture second(100 + rank);
  first.post(firstImpl, 105, 106);
  second.post(secondImpl, 107, 108);

  int firstComplete  = 0;
  int secondComplete = 0;
  while (!firstComplete || !secondComplete) {
    TEST_EQUALITY(firstImpl.get_status(first.req, &firstComplete, MPI_STATUS_IGNORE), MPI_SUCCESS);
    secondComplete = second.req.completed;
  }

  TEST_EQUALITY(firstImpl.wait(first.req), MPI_SUCCESS);
  TEST_EQUALITY(secondImpl.wait(second.req), MPI_SUCCESS);
  checkResult(first, 0, out, success);
  checkResult(second, 100, out, success);
}

TEUCHOS_UNIT_TEST(Ialltofewv, zeroCountsMultipleRoots) {
  // Mutliple source/root pairs send nothing
  int rank = 0;
  int size = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::array<int, 2> send{rank, 100 + rank};
  std::array<int, 2> sendcounts{1, rank % 2};
  std::array<int, 2> sdispls{0, 1};
  std::array<int, 2> roots{0, 1};
  std::vector<int> recvcounts;
  std::vector<int> rdispls;
  std::vector<int> recv;
  if (rank == roots[0] || rank == roots[1]) {
    recvcounts.resize(size);
    rdispls.resize(size);
    int total = 0;
    for (int src = 0; src < size; ++src) {
      recvcounts[src] = rank == roots[0] ? 1 : src % 2;
      rdispls[src]    = total;
      total += recvcounts[src];
    }
    recv.resize(total, -1);
  }

  Ialltofewv impl;
  Ialltofewv::Req req;
  const int postErr = impl.post<false>(send.data(), sendcounts.data(), sdispls.data(), MPI_INT,
                                       recv.data(), recvcounts.data(), rdispls.data(), roots.data(),
                                       int(roots.size()), MPI_INT, 109, 110, MPI_COMM_WORLD, &req);
  TEST_EQUALITY(postErr, MPI_SUCCESS);
  TEST_EQUALITY(impl.wait(req), MPI_SUCCESS);

  if (rank == roots[0]) {
    for (int src = 0; src < size; ++src) TEST_EQUALITY(recv[src], src);
  } else if (rank == roots[1]) {
    int offset = 0;
    for (int src = 0; src < size; ++src) {
      if (src % 2) {
        TEST_EQUALITY(recv[offset], 100 + src);
        ++offset;
      }
    }
  }
}

TEUCHOS_UNIT_TEST(Ialltofewv, reusesCachedViews) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Ialltofewv impl;
  Fixture first(rank, 2);
  first.post(impl, 111, 112);
  TEST_EQUALITY(impl.wait(first.req), MPI_SUCCESS);

  CachedAllocationCounter counter;
  Fixture second(100 + rank);
  second.post(impl, 113, 114);
  TEST_EQUALITY(impl.wait(second.req), MPI_SUCCESS);
  TEST_EQUALITY(counter.globalCount(), 0);
  checkResult(second, 100, out, success);
}

TEUCHOS_UNIT_TEST(Ialltofewv, copyHasIndependentCache) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  Ialltofewv original;
  Fixture first(rank, 2);
  first.post(original, 115, 116);
  TEST_EQUALITY(original.wait(first.req), MPI_SUCCESS);

  Ialltofewv copy(original);
  CachedAllocationCounter counter;
  Fixture second(100 + rank);
  second.post(copy, 117, 118);
  TEST_EQUALITY(copy.wait(second.req), MPI_SUCCESS);
  TEST_ASSERT(counter.globalCount() > 0);
  checkResult(second, 100, out, success);
}

}  // namespace

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard scopeGuard(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
