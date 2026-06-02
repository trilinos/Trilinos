// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_HeuristicBlockPermutation.hpp"
#include <Teuchos_UnitTestHarness.hpp>
#include "Teko_Utilities.hpp"
#include <random>
#include <numeric>

namespace Teko {

using Teuchos::ParameterList;
using Teuchos::RCP;
using TpetraComm = RCP<const Teuchos::Comm<int>>;
using Teko::GO;
using Teko::LO;
using Teko::NT;
using Teko::ST;

std::random_device r;
std::default_random_engine random_engine(r());

namespace {

using BlockGIDSMappingType   = std::vector<std::vector<GO>>;
using MatrixAndBlockGIDSPair = std::pair<Teuchos::RCP<Tpetra::CrsMatrix<>>, BlockGIDSMappingType>;
using ConnectionsType        = std::map<int, std::set<int>>;
MatrixAndBlockGIDSPair generate_random_matrix(const TpetraComm& comm,
                                              const std::vector<int>& blockSizes,
                                              const ConnectionsType& couplings) {
  std::vector<int> blockStarts(blockSizes.size());

  // use of std::exclusive_scan causes errors on NVIDIA builds
  blockStarts[0] = 0;
  for (auto blockNum = 1U; blockNum < blockSizes.size(); ++blockNum) {
    blockStarts[blockNum] = blockStarts[blockNum - 1] + blockSizes[blockNum - 1];
  }

  GO gid = 0;
  BlockGIDSMappingType blockGids(blockSizes.size());
  for (auto bid = 0U; bid < blockSizes.size(); ++bid) {
    blockGids[bid].resize(blockSizes[bid]);
    std::iota(blockGids[bid].begin(), blockGids[bid].end(), gid);
    gid = blockGids[bid].back() + 1;
  }

  const auto n                    = std::accumulate(blockSizes.begin(), blockSizes.end(), 0);
  RCP<const Tpetra::Map<>> rowMap = rcp(new const Tpetra::Map<>(n, 0, comm));

  std::vector<GO> indices(n);
  std::iota(indices.begin(), indices.end(), 0);

  std::vector<ST> randomRowValues(n);

  // small 'noise' to represent coupling between coupling between diffusion blocks
  std::uniform_real_distribution<ST> uniform_dist(-0.1, 0.0);

  auto mat = Teuchos::rcp(new Tpetra::CrsMatrix<>(rowMap, rowMap, n));

  GO row = 0;
  for (auto blockNum = 0U; blockNum < blockSizes.size(); ++blockNum) {
    const auto blockSize = blockSizes[blockNum];
    for (auto blockRow = 0; blockRow < blockSize; ++blockRow) {
      std::fill(randomRowValues.begin(), randomRowValues.end(), 0.0);
      if (couplings.count(blockNum) > 0) {
        for (const auto coupledBlock : couplings.at(blockNum)) {
          std::generate_n(randomRowValues.begin() + blockStarts[coupledBlock],
                          blockSizes[coupledBlock], [&]() { return uniform_dist(random_engine); });
        }
      }

      randomRowValues[row] = 2.0;
      if (blockRow - 1 >= 0) randomRowValues[row - 1] = -1.0;
      if (blockRow + 1 < blockSize) randomRowValues[row + 1] = -1.0;
      mat->insertGlobalValues(row, Teuchos::ArrayView<GO>(indices.data(), n),
                              Teuchos::ArrayView<ST>(randomRowValues.data(), n));
      row++;
    }
  }
  mat->fillComplete();
  return std::make_pair(mat, blockGids);
}

}  // namespace

class BlockReorderingFixture {
 public:
  BlockReorderingFixture() { SetUp(); }

  void SetUp() {
    this->comm = Tpetra::getDefaultComm();

    std::uniform_int_distribution<int> uniform_dist(1, 10);
    blockSizes.resize(nBlockRows);
    std::generate(blockSizes.begin(), blockSizes.end(),
                  [&]() { return uniform_dist(random_engine); });

    ConnectionsType all_to_all_connections;
    for (int row = 0; row < nBlockRows; ++row) {
      for (int connectedRow = 0; connectedRow < nBlockRows; ++connectedRow) {
        if (row == connectedRow) continue;
        all_to_all_connections[row].insert(connectedRow);
      }
    }

    std::tie(this->A, this->blockGids) =
        generate_random_matrix(comm, blockSizes, all_to_all_connections);

    this->tpetra_A_blk =
        Teuchos::make_rcp<Teko::TpetraHelpers::BlockedTpetraOperator>(this->blockGids, this->A);
  }

  void make_block_upper_triangular() {
    ConnectionsType block_ut_connections;
    for (int row = 0; row < nBlockRows; ++row) {
      for (int connectedRow = 0; connectedRow < nBlockRows; ++connectedRow) {
        if (row == connectedRow) continue;
        if (row > connectedRow) continue;
        block_ut_connections[row].insert(connectedRow);
      }
    }

    this->A.reset();
    this->tpetra_A_blk.reset();
    std::tie(this->A, this->blockGids) =
        generate_random_matrix(this->comm, this->blockSizes, block_ut_connections);

    this->tpetra_A_blk =
        Teuchos::make_rcp<Teko::TpetraHelpers::BlockedTpetraOperator>(this->blockGids, this->A);
  }

  void make_connections(const ConnectionsType& connections) {
    this->A.reset();
    this->tpetra_A_blk.reset();
    std::tie(this->A, this->blockGids) =
        generate_random_matrix(this->comm, this->blockSizes, connections);
    this->tpetra_A_blk =
        Teuchos::make_rcp<Teko::TpetraHelpers::BlockedTpetraOperator>(this->blockGids, this->A);
  }

  const int nBlockRows = 10;
  Teuchos::RCP<const Teuchos::Comm<int>> comm;
  std::vector<int> blockSizes;
  Teuchos::RCP<Tpetra::CrsMatrix<>> A;
  BlockGIDSMappingType blockGids;
  Teuchos::RCP<Teko::TpetraHelpers::BlockedTpetraOperator> tpetra_A_blk;
};

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation, LossMinimizationUTOptimality) {
  BlockReorderingFixture test_fixture{};
  test_fixture.make_block_upper_triangular();

  Teuchos::ParameterList p;
  p.set("Heuristic Method", "Loss Minimizing Reordering");
  p.set("Block Inverse Type", "Block Gauss-Seidel");
  p.set("Use Block Upper Triangle", true);

  // allow sufficient time to find optimal configuration for testing purposes
  const double maxWalltime = 1e-2;
  p.set("Max Heuristic Walltime", maxWalltime);

  auto A_blk                = Teko::toBlockedLinearOp(test_fixture.tpetra_A_blk->getThyraOp());
  auto [permutation, score] = generate_heuristic_permutation(A_blk, Teuchos::rcpFromRef(p));
  const auto tol            = std::numeric_limits<double>::epsilon();
  TEST_FLOATING_EQUALITY(score, 0.0, tol);

  // Initial ordering is already optimal for UT
  for (const auto& [row, permutedRow] : permutation) {
    TEST_EQUALITY(row, permutedRow);
  }
}

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation, LossMinimizationLTOptimality) {
  BlockReorderingFixture test_fixture{};
  test_fixture.make_block_upper_triangular();

  Teuchos::ParameterList p;
  p.set("Heuristic Method", "Loss Minimizing Reordering");
  p.set("Block Inverse Type", "Block Gauss-Seidel");
  p.set("Use Block Upper Triangle", false);

  // allow sufficient time to find optimal configuration for testing purposes
  const double maxWalltime = 1e-2;
  p.set("Max Heuristic Walltime", maxWalltime);

  auto A_blk                = Teko::toBlockedLinearOp(test_fixture.tpetra_A_blk->getThyraOp());
  auto [permutation, score] = generate_heuristic_permutation(A_blk, Teuchos::rcpFromRef(p));
  const auto tol            = std::numeric_limits<double>::epsilon();
  TEST_FLOATING_EQUALITY(score, 0.0, tol);
  const int nBlockRows = Teko::blockRowCount(A_blk);

  // Optimal ordering for LT is reversed
  for (const auto& [row, permutedRow] : permutation) {
    TEST_EQUALITY(permutedRow, nBlockRows - row - 1);
  }
}

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation, GreedyBlockMergingFuseIntoSingleBlock) {
  BlockReorderingFixture test_fixture{};

  Teuchos::ParameterList p;
  p.set("Heuristic Method", "Greedy Block Merging Heuristic");
  p.set("Block Inverse Type", "Block Gauss-Seidel");
  p.set("Use Block Upper Triangle", true);

  // allow sufficient time to find optimal configuration for testing purposes
  const double maxWalltime = 1e-2;
  p.set("Max Heuristic Walltime", maxWalltime);

  // norm loss is set such that we are guaranteed only one, monolithic sub-block
  const double targetNormLoss = 0.0;
  p.set("Target Norm Loss", targetNormLoss);

  auto A_blk                = Teko::toBlockedLinearOp(test_fixture.tpetra_A_blk->getThyraOp());
  auto [permutation, score] = generate_heuristic_permutation(A_blk, Teuchos::rcpFromRef(p));
  const auto tol            = std::numeric_limits<double>::epsilon();
  TEST_FLOATING_EQUALITY(score, 0.0, tol);

  for (const auto& [row, permutedRow] : permutation) {
    TEST_EQUALITY(permutedRow, 0);
  }
}

bool run_sparse_connectivity_test(const ConnectionsType& connections, int expectedBlockFuses) {
  bool test_passes = true;
  BlockReorderingFixture test_fixture{};

  std::set<std::pair<int, int>> fullyConnectedBlocks;
  for (const auto& [block, connectedBlocks] : connections) {
    for (const auto& connectedBlock : connectedBlocks) {
      if (connections.count(connectedBlock) == 0) continue;

      for (const auto& opposingConnectedBlock : connections.at(connectedBlock)) {
        if (opposingConnectedBlock == block) {
          const auto minBlock = std::min(block, connectedBlock);
          const auto maxBlock = std::max(block, connectedBlock);
          fullyConnectedBlocks.emplace(std::make_pair(minBlock, maxBlock));
        }
      }
    }
  }

  test_fixture.make_connections(connections);

  Teuchos::ParameterList p;
  p.set("Heuristic Method", "Greedy Block Merging Heuristic");
  p.set("Block Inverse Type", "Block Gauss-Seidel");
  p.set("Use Block Upper Triangle", true);

  // allow sufficient time to find optimal configuration for testing purposes
  const double maxWalltime = 1e-2;
  p.set("Max Heuristic Walltime", maxWalltime);

  // norm loss is set such that we are guaranteed only one, monolithic sub-block
  const double targetNormLoss = 0.0;
  p.set("Target Norm Loss", targetNormLoss);

  auto A_blk           = Teko::toBlockedLinearOp(test_fixture.tpetra_A_blk->getThyraOp());
  const int nBlockRows = Teko::blockRowCount(A_blk);

  auto [permutation, score] = generate_heuristic_permutation(A_blk, Teuchos::rcpFromRef(p));
  auto inverseMapping       = generate_inverse_permutation(permutation);

  // check that doubly-connected blocks were merged
  for (const auto& [block, connectedBlock] : fullyConnectedBlocks) {
    test_passes &= permutation.at(block) == permutation.at(connectedBlock);
  }

  const int nBlockFinal       = inverseMapping.size();
  const int expectedNumBlocks = nBlockRows - expectedBlockFuses;
  test_passes &= nBlockFinal == expectedNumBlocks;
  return test_passes;
}

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation, GreedyBlockMergingSingleConnectedBlock) {
  ConnectionsType connections{{1, {0}}, {0, {1}}};
  TEST_EQUALITY(run_sparse_connectivity_test(connections, 1), true);
}

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation,
                  GreedyBlockMergingSingleConnectedBlockWithTrivialConnection) {
  ConnectionsType connections{{1, {0, 2}}, {0, {1}}};
  TEST_EQUALITY(run_sparse_connectivity_test(connections, 1), true);
}

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation, GreedyBlockMergingTwoConnectedBlocks) {
  ConnectionsType connections{
      {1, {0}},
      {0, {1}},
      {3, {4}},
      {4, {3}},
  };
  TEST_EQUALITY(run_sparse_connectivity_test(connections, 2), true);
}

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation,
                  GreedyBlockMergingTwoConnectedBlocksWithTrivialConnection) {
  ConnectionsType connections{
      {1, {0}}, {0, {1}}, {0, {2}}, {3, {4}}, {4, {3}},
  };
  TEST_EQUALITY(run_sparse_connectivity_test(connections, 2), true);
}

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation,
                  GreedyBlockMergingTwoConnectedBlocksWithSeveralTrivialConnection) {
  ConnectionsType connections{
      {1, {0}}, {1, {9}}, {0, {1}}, {0, {2}}, {0, {7}},
      {0, {8}}, {3, {4}}, {4, {3}}, {5, {6}}, {6, {8}},
  };
  TEST_EQUALITY(run_sparse_connectivity_test(connections, 2), true);
}

TEUCHOS_UNIT_TEST(tHeuristicBlockPermutation, GreedyBlockMergingThreeFullyConnectedBlocks) {
  ConnectionsType connections{{2, {0, 1}}, {1, {0, 2}}, {0, {1, 2}}};
  TEST_EQUALITY(run_sparse_connectivity_test(connections, 2), true);
}

}  // namespace Teko
