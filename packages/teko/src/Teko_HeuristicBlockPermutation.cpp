// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teko_BlockedReordering.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Teko_Utilities.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teko_HeuristicBlockPermutation.hpp"
#include "Teko_AdaptivePreconditionerFactory.hpp"
#include <algorithm>
#include <numeric>
#include <queue>
#include <utility>
#include <random>

using Teuchos::RCP;

namespace Teko {

namespace {

constexpr bool is_gpu_build() { return Teko::NT::is_gpu; }

Teuchos::RCP<const Teuchos::Comm<int>> extract_communicator(const Teko::BlockedLinearOp& blocked) {
  bool transposed = false;
  Teko::ST scalar = 0.0;
  auto A_00 =
      Teko::TpetraHelpers::getTpetraCrsMatrix(Teko::getBlock(0, 0, blocked), &scalar, &transposed);
  auto rowMap_ij = A_00->getRowMap();
  return rowMap_ij->getComm();
}

double upperTriangularLossScore(const BlockNormsViewType& blockNorms,
                                const std::vector<int>& tentativePermutation) {
  double loss = 0.0;
  const int N = blockNorms.extent(0);
  for (int row = 0; row < N; row++) {
    const auto row_idx = tentativePermutation[row];
    for (int col = 0; col < row; col++) {
      const auto col_idx = tentativePermutation[col];
      const auto v       = blockNorms(row_idx, col_idx);
      loss += v * v;
    }
  }
  return loss;
}

double lowerTriangularLossScore(const BlockNormsViewType& blockNorms,
                                const std::vector<int>& tentativePermutation) {
  double loss = 0.0;
  const int N = blockNorms.extent(0);
  for (int row = 0; row < N; row++) {
    const auto row_idx = tentativePermutation[row];
    for (int col = row + 1; col < N; col++) {
      const auto col_idx = tentativePermutation[col];
      const auto v       = blockNorms(row_idx, col_idx);
      loss += v * v;
    }
  }
  return loss;
}

struct HeuristicSettings {
  std::string heuristicMethod{"Greedy Block Merging Heuristic"};
  std::string blockInverseType{"Block Gauss-Seidel"};
  bool gsUseUpperTriangle{true};
  double maxHeuristicWalltime{1e-2};
  double targetNormLoss{1e-3};

  HeuristicMethodFunctionType userDefinedPermutation;
};

HeuristicSettings parse_parameters_from_pl(Teuchos::RCP<Teuchos::ParameterList> pl) {
  HeuristicSettings params;
  if (!pl) return params;

  if (pl->isParameter("Heuristic Method")) {
    params.heuristicMethod = pl->get<std::string>("Heuristic Method");
  }

  if (pl->isParameter("Max Heuristic Walltime")) {
    params.maxHeuristicWalltime = pl->get<double>("Max Heuristic Walltime");
  }

  if (pl->isParameter("Target Norm Loss")) {
    params.targetNormLoss = pl->get<double>("Target Norm Loss");
  }

  if (pl->isParameter("Block Inverse Type")) {
    params.blockInverseType = pl->get<std::string>("Block Inverse Type");
  }

  if (pl->isParameter("Use Block Upper Triangle")) {
    params.gsUseUpperTriangle = pl->get<bool>("Use Block Upper Triangle");
  }

  if (!pl->isSublist("user data")) return params;

  const auto& user_data = pl->sublist("user data");
  if (user_data.isParameter("User Defined Permutation Function")) {
    params.userDefinedPermutation =
        user_data.get<HeuristicMethodFunctionType>("User Defined Permutation Function");
  }

  return params;
}

PermutationScoreType lossMinimizationBlockPermutation(
    const BlockNormsViewType& blockNorms,
    const Teuchos::RCP<const Teuchos::Comm<int>>& communicator, const HeuristicSettings& settings) {
  const int N = blockNorms.extent(0);

  auto objective_function = [&](const std::vector<int>& indices) {
    return settings.gsUseUpperTriangle ? upperTriangularLossScore(blockNorms, indices)
                                       : lowerTriangularLossScore(blockNorms, indices);
  };

  auto [ordering, obj_f] = LinearOrdering::compute_min_ordering(
      communicator, settings.maxHeuristicWalltime, blockNorms, settings.gsUseUpperTriangle,
      objective_function, Teuchos::null);

  PermutationType reordering;
  for (int i = 0; i < N; ++i) {
    reordering[i] = ordering[i];
  }

  return std::make_pair(reordering, obj_f);
}

using BlockSet         = std::set<int>;
using Index            = std::pair<BlockSet, BlockSet>;
using BlockOrderMap    = std::map<BlockSet, int>;
using BlockNormMapType = std::map<Index, double>;

Index computeMaxOffDiagNorm(const BlockNormMapType& norms, const BlockOrderMap& ordering,
                            bool upper_triangular) {
  double maxOffDiagNorm = 0.0;
  Index maxIndices;
  for (auto&& [indices, norm] : norms) {
    const auto row = ordering.at(indices.first);
    const auto col = ordering.at(indices.second);

    // want largest off-diagonal sub-block norm in the section of the matrix _not_ encapsulated by
    // the preconditioner
    if (col >= row && upper_triangular) continue;
    if (row >= col && !upper_triangular) continue;

    if (norm > maxOffDiagNorm) {
      maxIndices     = indices;
      maxOffDiagNorm = norm;
    }
  }

  return maxIndices;
}

std::pair<BlockNormMapType, BlockOrderMap> combineBlocks(const Index& indexStar,
                                                         const BlockNormMapType& norms,
                                                         const BlockOrderMap& orderMap) {
  const auto& [istar, jstar] = indexStar;
  BlockNormMapType newNorms;
  BlockOrderMap newOrderMap;

  const auto index_istar = orderMap.at(istar);
  const auto index_jstar = orderMap.at(jstar);
  int min_index          = std::min(index_istar, index_jstar);
  int max_index          = std::max(index_istar, index_jstar);

  for (const auto& [indexSet, index] : orderMap) {
    if (indexSet == jstar || indexSet == istar) {
      continue;
    }
    const auto newIndex   = index > max_index ? index - 1 : index;
    newOrderMap[indexSet] = newIndex;
  }

  // populate norms for non-merged entries
  for (const auto& [indices, norm] : norms) {
    const auto& [i, j] = indices;
    if (i == istar) continue;
    if (i == jstar) continue;
    if (j == istar) continue;
    if (j == jstar) continue;
    newNorms[indices] = norm;
  }

  // merged diagonal block
  BlockSet new_entry;
  std::set_union(istar.cbegin(), istar.cend(), jstar.cbegin(), jstar.cend(),
                 std::inserter(new_entry, new_entry.end()));
  newOrderMap[new_entry] = min_index;

  const auto norm_A_ii = norms.at(std::make_pair(istar, istar));
  const auto norm_A_ij = norms.at(std::make_pair(istar, jstar));
  const auto norm_A_ji = norms.at(std::make_pair(jstar, istar));
  const auto norm_A_jj = norms.at(std::make_pair(jstar, jstar));

  double blockFrobeniusNorm = 0.0;
  blockFrobeniusNorm += norm_A_ii * norm_A_ii;
  blockFrobeniusNorm += norm_A_ij * norm_A_ij;
  blockFrobeniusNorm += norm_A_ji * norm_A_ji;
  blockFrobeniusNorm += norm_A_jj * norm_A_jj;
  newNorms[std::make_pair(new_entry, new_entry)] = std::sqrt(blockFrobeniusNorm);

  // cross-terms
  for (const auto& [entry, index] : orderMap) {
    if (entry == istar) continue;
    if (entry == jstar) continue;

    newNorms[std::make_pair(entry, new_entry)] =
        std::hypot(norms.at(std::make_pair(entry, istar)), norms.at(std::make_pair(entry, jstar)));

    newNorms[std::make_pair(new_entry, entry)] =
        std::hypot(norms.at(std::make_pair(istar, entry)), norms.at(std::make_pair(jstar, entry)));
  }

  return std::make_pair(newNorms, newOrderMap);
}

BlockNormsViewType blockNormsFromNormMap(const BlockNormMapType& normMap,
                                         const BlockOrderMap& ordering) {
  const auto numBlocksFinal = ordering.size();
  BlockNormsViewType fusedBlockNorms("fusedBlockNorms", numBlocksFinal, numBlocksFinal);

  for (const auto& [i_index, i_blk] : ordering) {
    for (const auto& [j_index, j_blk] : ordering) {
      const auto v                  = normMap.at(std::make_pair(i_index, j_index));
      fusedBlockNorms(i_blk, j_blk) = v;
    }
  }

  return fusedBlockNorms;
}

BlockOrderMap construct_reordered_map(const std::map<int, int>& newOrdering,
                                      BlockOrderMap& orderMap) {
  BlockOrderMap newOrderMap;
  for (const auto& [indexSet, index] : orderMap) {
    newOrderMap[indexSet] = newOrdering.at(index);
  }
  return newOrderMap;
}

double compute_norm_A(const BlockNormsViewType& blockNorms) {
  double norm = 0.0;
  const int N = blockNorms.extent(0);
  for (int row = 0; row < N; ++row) {
    for (int col = 0; col < N; ++col) {
      const auto v = blockNorms(row, col);
      norm += v * v;
    }
  }
  return std::sqrt(norm);
}

PermutationScoreType greedyBlockMergingBlockPermutation(
    const BlockNormsViewType& blockNorms,
    const Teuchos::RCP<const Teuchos::Comm<int>>& communicator,
    const HeuristicSettings& heuristicSettings) {
  const int N             = blockNorms.extent(0);
  const auto tMaxWalltime = heuristicSettings.maxHeuristicWalltime / N;

  auto reorderingMinimizationSettings                 = heuristicSettings;
  reorderingMinimizationSettings.maxHeuristicWalltime = tMaxWalltime;
  auto [reordering, initialScore] =
      lossMinimizationBlockPermutation(blockNorms, communicator, reorderingMinimizationSettings);

  const auto norm_A = compute_norm_A(blockNorms);
  initialScore *= norm_A;

  if (initialScore == 0) {  // lucky convergence, i.e., one-way coupled problem
    return std::make_pair(reordering, initialScore);
  }

  PermutationType permutation;

  BlockNormMapType normMap;
  BlockOrderMap orderMap;
  for (int row = 0; row < N; ++row) {
    orderMap[{row}] = reordering[row];
    for (int col = 0; col < N; ++col) {
      normMap[{{row}, {col}}] = blockNorms(row, col);
    }
  }

  double currentScore = initialScore;

  auto indexStar = computeMaxOffDiagNorm(normMap, orderMap, heuristicSettings.gsUseUpperTriangle);

  while (normMap.size() > 1 && currentScore > heuristicSettings.targetNormLoss) {
    std::tie(normMap, orderMap)    = combineBlocks(indexStar, normMap, orderMap);
    auto newBlockNorms             = blockNormsFromNormMap(normMap, orderMap);
    auto [orderPermutation, score] = lossMinimizationBlockPermutation(
        newBlockNorms, communicator, reorderingMinimizationSettings);
    currentScore = score * norm_A;
    orderMap     = construct_reordered_map(orderPermutation, orderMap);
    indexStar    = computeMaxOffDiagNorm(normMap, orderMap, heuristicSettings.gsUseUpperTriangle);
  }

  const auto finalScore = currentScore;

  for (const auto& [indexSet, newIndex] : orderMap) {
    for (auto&& index : indexSet) {
      permutation[index] = newIndex;
    }
  }

  return std::make_pair(permutation, finalScore);
}

}  // namespace

InversePermutationType generate_inverse_permutation(const PermutationType& permutation) {
  InversePermutationType inverseMapping;
  for (auto&& [row, permutedRow] : permutation) {
    inverseMapping[permutedRow].insert(row);
  }
  return inverseMapping;
}

PermutationScoreType generate_heuristic_permutation(
    const Teko::BlockedLinearOp& dof_blocked_matrix,
    Teuchos::RCP<Teuchos::ParameterList> parameters) {
  int blockRows = Teko::blockRowCount(dof_blocked_matrix);
  int blockCols = Teko::blockColCount(dof_blocked_matrix);
  TEUCHOS_ASSERT(blockRows == blockCols);

  BlockNormsViewType blockNorms("Block Norms", blockRows, blockCols);
  for (int blockRow = 0; blockRow < blockRows; ++blockRow) {
    for (int blockCol = 0; blockCol < blockCols; ++blockCol) {
      blockNorms(blockRow, blockCol) = 0.0;
      auto A_ij                      = Teko::getBlock(blockRow, blockCol, dof_blocked_matrix);
      if (A_ij != Teuchos::null) {
        blockNorms(blockRow, blockCol) = Teko::frobeniusNorm(A_ij);
      }
    }
  }

  auto comm = extract_communicator(dof_blocked_matrix);

  return generate_heuristic_permutation(blockNorms, comm, parameters);
}

PermutationScoreType generate_heuristic_permutation(
    const Teuchos::RCP<Teko::TpetraHelpers::BlockedTpetraOperator>& dof_blocked_matrix,
    Teuchos::RCP<Teuchos::ParameterList> parameters) {
  using Teko::GO;
  using Teko::LO;
  using Teko::NT;
  using Teko::ST;

  int blockRows = dof_blocked_matrix->GetBlockRowCount();
  int blockCols = dof_blocked_matrix->GetBlockColCount();
  TEUCHOS_ASSERT(blockRows == blockCols);
  BlockNormsViewType blockNorms("Block Norms", blockRows, blockCols);

  for (int blockRow = 0; blockRow < blockRows; ++blockRow) {
    for (int blockCol = 0; blockCol < blockCols; ++blockCol) {
      blockNorms(blockRow, blockCol) = 0.0;
      const auto A_ij = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<ST, LO, GO, NT>>(
          dof_blocked_matrix->GetBlock(blockRow, blockCol), true);
      if (A_ij) {
        blockNorms(blockRow, blockCol) = A_ij->getFrobeniusNorm();
      }
    }
  }

  return generate_heuristic_permutation(blockNorms, Teuchos::rcpFromRef(dof_blocked_matrix->Comm()),
                                        parameters);
}

namespace {
auto generate_method_map(const HeuristicSettings& settings) {
  auto greedy_merging_func = [settings](
                                 const BlockNormsViewType& blockNorms,
                                 const Teuchos::RCP<const Teuchos::Comm<int>>& communicator) {
    return greedyBlockMergingBlockPermutation(blockNorms, communicator, settings);
  };

  auto loss_minimizing_reordering_func =
      [settings](const BlockNormsViewType& blockNorms,
                 const Teuchos::RCP<const Teuchos::Comm<int>>& communicator) {
        return lossMinimizationBlockPermutation(blockNorms, communicator, settings);
      };

  auto methods = std::map<std::string, HeuristicMethodFunctionType>{
      {"Greedy Block Merging Heuristic", greedy_merging_func},
      {"Loss Minimizing Reordering", loss_minimizing_reordering_func}};

  if (settings.userDefinedPermutation) {
    methods["User Defined Heuristic"] = settings.userDefinedPermutation;
  }
  return methods;
}
}  // namespace

PermutationScoreType generate_heuristic_permutation(
    const BlockNormsViewType& blockNorms, Teuchos::RCP<const Teuchos::Comm<int>> communicator,
    Teuchos::RCP<Teuchos::ParameterList> parameters) {
  auto heuristicSettings = parse_parameters_from_pl(parameters);
  auto methodMap         = generate_method_map(heuristicSettings);
  if (methodMap.count(heuristicSettings.heuristicMethod) == 0) {
    std::ostringstream ss;
    ss << "Error: could not find method named " << heuristicSettings.heuristicMethod
       << " in methods available to generate_heuristic_permutation.\n";
    ss << "Available methods are:\n";
    for (const auto& [methodName, _] : methodMap) {
      ss << "\t" << methodName << "\n";
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidArgument, ss.str());
  }
  auto method = methodMap[heuristicSettings.heuristicMethod];
  return method(blockNorms, communicator);
}

std::vector<std::vector<Teko::GO>> construct_block_gids_from_permutation(
    const PermutationType& permutation, const std::vector<std::vector<Teko::GO>>& dof_gids) {
  auto inverse_permutation = generate_inverse_permutation(permutation);
  return construct_block_gids_from_permutation(inverse_permutation, dof_gids);
}

std::vector<std::vector<Teko::GO>> construct_block_gids_from_permutation(
    const InversePermutationType& inverse_permutation,
    const std::vector<std::vector<Teko::GO>>& dof_gids) {
  const int numBlocksFinal = inverse_permutation.size();
  std::vector<std::vector<Teko::GO>> block_gids;
  block_gids.resize(numBlocksFinal);
  for (int block = 0; block < numBlocksFinal; ++block) {
    size_t numGIDBlock = 0;
    for (auto&& dof : inverse_permutation.at(block)) {
      numGIDBlock += dof_gids[dof].size();
    }
    block_gids[block].reserve(numGIDBlock);

    for (auto&& dof : inverse_permutation.at(block)) {
      block_gids[block].insert(block_gids[block].end(), dof_gids[dof].begin(), dof_gids[dof].end());
    }
    TEUCHOS_ASSERT(block_gids[block].size() == numGIDBlock);
  }

  return block_gids;
}

namespace {
bool is_inverse_type(const std::string& name) {
  bool is_inverse = false;
  is_inverse |= name.find("Inverse") != std::string::npos;
  is_inverse |= name.find("Preconditioner") != std::string::npos;
  is_inverse &= name.find("Type") != std::string::npos;
  return is_inverse;
}

bool is_package_name(const std::string& name) {
  bool is_package = false;
  is_package |= name.find("Amesos") != std::string::npos;
  is_package |= name.find("Amesos2") != std::string::npos;
  is_package |= name.find("Ifpack") != std::string::npos;
  is_package |= name.find("Ifpack2") != std::string::npos;
  is_package |= name.find("MueLu") != std::string::npos;
  is_package |= name.find("Belos") != std::string::npos;
  return is_package;
}

Teuchos::RCP<Teuchos::ParameterList> mangle_subblock_inverses(const Teuchos::ParameterList& params,
                                                              std::string suffix) {
  auto newParams = Teuchos::make_rcp<Teuchos::ParameterList>();

  for (auto& param : params) {
    if (!params.isSublist(param.key)) continue;
    const auto& sublist = params.sublist(param.key);
    if (!sublist.isParameter("Type")) continue;

    const auto newSublistName = param.key + suffix;
    auto& newSublist          = newParams->sublist(newSublistName);
    newSublist.setParameters(sublist);
  }

  for (auto& param : *newParams) {
    if (!newParams->isSublist(param.key)) continue;
    auto& sublist = newParams->sublist(param.key);
    if (!sublist.isParameter("Type")) continue;
    for (auto&& sublistParam : sublist) {
      if (!is_inverse_type(sublistParam.key)) continue;
      const auto inverseMethod = sublist.get<std::string>(sublistParam.key);
      if (is_package_name(inverseMethod)) continue;
      sublist.set(sublistParam.key, inverseMethod + suffix);
    }
  }

  return newParams;
}

std::string block_suffix(int block) { return "__" + std::to_string(block + 1); }

Teuchos::RCP<Teuchos::ParameterList> generate_parameters_from_permutation(
    const PermutationType& permutation, const InversePermutationType& inverse_permutation,
    std::string inverseName, Teuchos::RCP<Teuchos::ParameterList> heuristicParameters,
    SubblockParameters subblockParameters) {
  const auto heuristicSettings = parse_parameters_from_pl(heuristicParameters);

  Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::make_rcp<Teuchos::ParameterList>();
  auto& blockInverseParams                          = solverParams->sublist(inverseName);
  blockInverseParams.set("Type", heuristicSettings.blockInverseType);
  if (heuristicSettings.blockInverseType == "Block Gauss-Seidel") {
    blockInverseParams.set("Use Upper Triangle", heuristicSettings.gsUseUpperTriangle);
  }

  for (const auto& [dofBlockNumber, originalDofBlockNumbers] : inverse_permutation) {
    auto subblockInverseParams        = subblockParameters.mergedSolverParameters;
    auto subblockPreconditionerParams = subblockParameters.mergedPreconditionerParameters;
    auto subblockInverseName          = subblockParameters.mergedSolverName;
    auto subblockPreconditionerName   = subblockParameters.mergedPreconditionerName;

    if (originalDofBlockNumbers.size() == 1) {
      const auto originalBlock = *originalDofBlockNumbers.begin();
      subblockInverseParams    = subblockParameters.subblockSolverParameters[originalBlock];
      subblockInverseName      = subblockParameters.subblockSolverNames[originalBlock];
      if (!subblockParameters.subblockPreconditionerParameters.empty()) {
        subblockPreconditionerParams =
            subblockParameters.subblockPreconditionerParameters[originalBlock];
        subblockPreconditionerName = subblockParameters.subblockPreconditionerNames[originalBlock];
      }
    }

    const auto suffix = block_suffix(dofBlockNumber);
    blockInverseParams.set("Inverse Type " + std::to_string(dofBlockNumber + 1),
                           subblockInverseName + suffix);

    auto mangledSubblockSolverParams = mangle_subblock_inverses(*subblockInverseParams, suffix);
    solverParams->setParameters(*mangledSubblockSolverParams);

    if (subblockPreconditionerParams) {
      blockInverseParams.set("Preconditioner Type " + std::to_string(dofBlockNumber + 1),
                             subblockPreconditionerName + suffix);
      auto mangledSubblockPreconditionerParams =
          mangle_subblock_inverses(*subblockPreconditionerParams, suffix);
      solverParams->setParameters(*mangledSubblockPreconditionerParams);
    }
  }

  return solverParams;
}

}  // namespace

std::tuple<Teuchos::RCP<Teuchos::ParameterList>, std::string> default_merged_solver_parameters() {
  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::make_rcp<Teuchos::ParameterList>();

  const std::string subblockSolverName{"default_subblock_solver"};

  const std::string jacobiType      = "jacobi_pc";
  const std::string ddiluType       = "ddilu_pc";
  const std::string gmresType       = "gmres";
  const std::string robustDDILUType = "robust_ddilu_pc";

  // Penultimate before resorting to direct solver
  const std::string reallyRobustDDILUType = "really_robust_ddilu_pc";

  // Solver of last resort
  const std::string directSolverType = "Amesos2";

  const double targetResidualReduction = 0.1;

  // default sub-block solver
  {
    auto& tekoBlockSolverParams = params->sublist(subblockSolverName, false);
    tekoBlockSolverParams.set("Type", "Adaptive");
    tekoBlockSolverParams.set("Target Residual Reduction", targetResidualReduction);
    tekoBlockSolverParams.set("Inverse Type 1", gmresType);
    tekoBlockSolverParams.set("Preconditioner Type 1", jacobiType);
    tekoBlockSolverParams.set("Inverse Type 2", gmresType);
    tekoBlockSolverParams.set("Preconditioner Type 2", ddiluType);
    tekoBlockSolverParams.set("Inverse Type 3", gmresType);
    tekoBlockSolverParams.set("Preconditioner Type 3", robustDDILUType);
    tekoBlockSolverParams.set("Inverse Type 4", gmresType);
    tekoBlockSolverParams.set("Preconditioner Type 4", reallyRobustDDILUType);

    // Direct sub-block solver must be limited to relatively small problem sizes.
    // Otherwise, several Amesos2 factorization failures may kill the simulation.
    using SizeOrdinalType                         = AdaptivePreconditionerFactory::SizeOrdinalType;
    const SizeOrdinalType maximumDirectSolverSize = 100000;
    tekoBlockSolverParams.set("Maximum Size 5", maximumDirectSolverSize);
    tekoBlockSolverParams.set("Inverse Type 5", directSolverType);
  }

  {
    auto& jacobiParams = params->sublist(jacobiType, false);
    jacobiParams.set("Type", "Ifpack2");
    jacobiParams.set("Prec Type", "relaxation");
    auto& ifpack2Params = jacobiParams.sublist("Ifpack2 Settings", false);
    ifpack2Params.set("relaxation: type", "Jacobi");
    ifpack2Params.set("relaxation: sweeps", 1);
    ifpack2Params.set("relaxation: damping factor", 1.0);
  }

  auto gen_schwarz_settings = [](int fillLevel, int overlapLevel) {
    Teuchos::ParameterList schwarzParams;
    schwarzParams.set("Type", "Ifpack2");
    schwarzParams.set("Prec Type", "schwarz");
    auto& ifpack2Settings = schwarzParams.sublist("Ifpack2 Settings", false);
    ifpack2Settings.set("subdomain solver name", "RILUK");
    auto& subdomainParams = ifpack2Settings.sublist("subdomain solver parameters", false);
    subdomainParams.set("fact: iluk level-of-fill", fillLevel);
    ifpack2Settings.set("schwarz: overlap level", overlapLevel);
    ifpack2Settings.set("schwarz: use reordering", true);
    auto& reorderingParams = ifpack2Settings.sublist("schwarz: reordering list", false);
    reorderingParams.set("order_method", "rcm");
    return schwarzParams;
  };

  {
    const int fill    = 0;
    const int overlap = 0;
    auto& ddiluParams = params->sublist(ddiluType, false);
    ddiluParams.setParameters(gen_schwarz_settings(fill, overlap));
  }

  {
    const int fill    = 1;
    const int overlap = 1;
    auto& ddiluParams = params->sublist(robustDDILUType, false);
    ddiluParams.setParameters(gen_schwarz_settings(fill, overlap));
  }

  {
    const int fill    = 2;
    const int overlap = 2;
    auto& ddiluParams = params->sublist(reallyRobustDDILUType, false);
    ddiluParams.setParameters(gen_schwarz_settings(fill, overlap));
  }

  const auto gmresTol = 0.9 * targetResidualReduction;
  {
    int maxIterations = 30;
    if (is_gpu_build()) {
      // Sub-block preconditioner construction costs can become dominant on the GPU,
      // especially for DD-ILU based approaches.
      // Bumping up the maximum numer of GMRES iterations is therefore warranted.
      maxIterations = 50;
    }

    auto& gmresParams = params->sublist(gmresType, false);
    gmresParams.set("Type", "Belos");
    gmresParams.set("Solver Type", "Pseudo Block GMRES");
    auto& settings =
        gmresParams.sublist("Solver Types", false).sublist("Pseudo Block GMRES", false);
    settings.set("Maximum Iterations", maxIterations);
    settings.set("Num Blocks", maxIterations);
    settings.set("Maximum Restarts", 1);
    settings.set("Orthogonalization", "ICGS");
    settings.set("Convergence Tolerance", gmresTol);
    settings.set("Implicit Residual Scaling", "Norm of RHS");
    settings.set("Explicit Residual Scaling", "Norm of RHS");
  }

  return std::make_tuple(params, subblockSolverName);
}

Teuchos::RCP<Teuchos::ParameterList> generate_parameters_from_permutation(
    const PermutationType& permutation, std::string inverseName,
    Teuchos::RCP<Teuchos::ParameterList> heuristicParameters,
    SubblockParameters subblockParameters) {
  const auto inversePermutation = generate_inverse_permutation(permutation);

  auto [defaultSolverParams, defaultSolverName] = default_merged_solver_parameters();
  if (!subblockParameters.mergedSolverParameters) {
    subblockParameters.mergedSolverName       = defaultSolverName;
    subblockParameters.mergedSolverParameters = defaultSolverParams;
  }

  if (subblockParameters.subblockSolverParameters.empty()) {
    const auto nBlock = permutation.size();
    subblockParameters.subblockSolverNames.resize(nBlock);
    subblockParameters.subblockSolverParameters.resize(nBlock);
  }

  // assign default solver for any subblocks that have yet to have one specified
  for (auto block = 0U; block < subblockParameters.subblockSolverParameters.size(); ++block) {
    if (subblockParameters.subblockSolverParameters[block]) continue;
    subblockParameters.subblockSolverNames[block]      = defaultSolverName;
    subblockParameters.subblockSolverParameters[block] = defaultSolverParams;
  }

  return generate_parameters_from_permutation(permutation, inversePermutation, inverseName,
                                              heuristicParameters, subblockParameters);
}

}  // namespace Teko
