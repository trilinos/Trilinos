// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_HeuristicBlockPermutation_hpp__
#define __Teko_HeuristicBlockPermutation_hpp__

#include "Kokkos_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Teko_BlockedTpetraOperator.hpp"
#include "Teko_Utilities.hpp"
#include <map>
#include <set>
#include <vector>

namespace Teko {

/// \brief A (dense) matrix type representing the multiphysics system via the coupling strength
/// between each physics block.
///
/// Typically, each entry corresponds to the Frobenius norm of the sub-block matrix:
/// \f\[
///   H_{i,j} := \|A_{i,j}\|_F
/// \f\].
using BlockNormsViewType = Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::HostSpace>;

/// \brief A mapping between the sub-block number correspondng to each physics and the sub-block
/// number in the (block) preconditioner. Note that a permutation need not be one-to-one. For
/// example, two sub-blocks could be "merged" into the same output sub-block.
using PermutationType = std::map<int, int>;

/// \brief Return type from generate_parameters_from_permutation.
///  The first entry is a PertmuationType specifying the mapping from sub-block physics number to
///  sub-block number in the preconditioner. The second entry is the score, defined as
///  \f$\|M-A\|_F\f$, where \f$A\f$ is the (block) matrix, and \f$M\f$ is the forward preconditioner
///  operator.
using PermutationScoreType = std::pair<PermutationType, double>;

/// \brief The inverse mapping from sub-block number in the preconditioner and the set of sub-block
/// physics.
using InversePermutationType = std::map<int, std::set<int>>;

/// \brief The inverse mapping from sub-block number in the preconditioner and the set of sub-block
/// physics. A user may provide an object of this type to use their own heuristic method of
/// generating a block permutation. See the documentation for generate_parameters_from_permutation
/// below for additional information.
using HeuristicMethodFunctionType = std::function<PermutationScoreType(
    const BlockNormsViewType &, const Teuchos::RCP<const Teuchos::Comm<int>> &)>;

/// \brief Given a dense matrix of sub-block norms with
/// \f\[
///   H_{i,j} := \|A_{i,j}\|_F
/// \f\],
/// compute a block permutation that seeks to improve the preconditioner quality.
/// By default, this means a block permutation that minimizes \f$\|M-A\|_F\f$ through a combination
/// of block merging and block reordering operations for block Gauss-Seidel.
///
/// NOTE: this must be called collectively by every rank in the communicator.
///
/// \param dof_blocked_matrix [in]  Blocked operator based on one-to-one physics-to-block mapping.
/// \param communicator  [in] Input communicator. This function must be called collectively by every
/// rank in the communicator. \param parameters  [in] Optional parameters specifying heuristic
/// algorithm settings.
///    These settings are:
///      - "Heuristic Method" : string representing the heuristic method. The currently accepted
///      values are:
///        1. "Greedy Block Merging Heuristic" (default).
///           This greedily merges blocks until \f$\|M-A\|_F \|A\|_F \leq \tau\f$.
///           \f$\tau\f$ may be specified through the "Target Norm Loss" option below.
///        2. "Loss Minimizing Reordering".
///           Find the (simultaneous) row and column permutation minimizing \f$\|M-A\|_F\f$
///        3. "User Defined Heuristic".
///           If a user defined heuristic function is supplied, use that.
///      - "Max Heuristic Walltime" : double representing how long (in seconds) to allow the
///      heuristic method to run.
///        By default, this is 1e-2 seconds.
///      - "Target Norm Loss" : double representing user-defined threshold for merging blocks.
///         See "Heuristic Method" -> "Greedy Block Merging Heuristic" above for additional
///         information.
///      - "Block Inverse Type" : string representing block splitting approach.
///        Currently, only "Block Gauss-Seidel" is supported.
///      - "Use Block Upper Triangle" : bool on whether to use upper triangular or lower triangular
///      Block Gauss-Seidel.
///      - "user data" : sublist with user-defined data, including:
///        -# "User Defined Permutation Function" : HeuristicMethodFunctionType-type for
///        user-defined permutation.
///
/// \return The permutation and score satisfying the heuristic settings from the parameters input.
PermutationScoreType generate_heuristic_permutation(
    const BlockNormsViewType &blockNorms, Teuchos::RCP<const Teuchos::Comm<int>> communicator,
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::null);

/// \brief This overload is the same as generate_heuristic_permutation above, but computes the dense
/// matrix of sub-block norms \f\[
///   H_{i,j} := \|A_{i,j}\|_F
/// \f\]
/// based on the input operator.
///
/// NOTE: this must be called collectively by every rank in the input operator.
///
/// \param dof_blocked_matrix [in]  Blocked operator based on one-to-one physics-to-block mapping.
/// \param parameters  [in] Optional parameters specifying heuristic algorithm settings.
///    For additional information, please consult the overload of generate_heuristic_permutation
///    above.
/// \return The permutation and score satisfying the heuristic settings from the parameters input.
PermutationScoreType generate_heuristic_permutation(
    const Teko::BlockedLinearOp &dof_blocked_matrix,
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::null);

/// \brief This overload is the same as generate_heuristic_permutation above, but computes the dense
/// matrix of sub-block norms \f\[
///   H_{i,j} := \|A_{i,j}\|_F
/// \f\]
/// based on the input operator.
///
/// NOTE: this must be called collectively by every rank in the input operator.
///
/// \param dof_blocked_matrix [in]  Blocked operator based on one-to-one physics-to-block mapping.
/// \param parameters  [in] Optional parameters specifying heuristic algorithm settings.
///    For additional information, please consult the overload of generate_heuristic_permutation
///    above.
/// \return The permutation and score satisfying the heuristic settings from the parameters input.
PermutationScoreType generate_heuristic_permutation(
    const Teuchos::RCP<Teko::TpetraHelpers::BlockedTpetraOperator> &dof_blocked_matrix,
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::null);

/// \brief Given the <tt>std::vector<std::vector<Teko::GO>></tt> representing the mapping of global
/// ordinals to individual physics,
///   construct the <tt>std::vector<std::vector<Teko::GO>></tt> representing the mapping of global
///   ordinals to blocks in the preconditioner.
///
/// \param permutation [in] Physics-to-sub-block mapping, as one may generate from calling
/// generate_heuristic_permutation above. \param parameters  [in] Mapping of global ordinals to
/// individual physics. \return Mapping of global ordinals to sub-blocks in the preconditioner.
std::vector<std::vector<Teko::GO>> construct_block_gids_from_permutation(
    const PermutationType &permutation, const std::vector<std::vector<Teko::GO>> &dof_gids);

/// \brief This overload is the same as construct_block_gids_from_permutation above, but uses the
/// inverse permutation.
///
/// \param permutation [in] Physics-to-sub-block mapping, as one may generate from calling
/// generate_heuristic_permutation above. \param parameters  [in] Mapping of global ordinals to
/// individual physics. \return Mapping of global ordinals to sub-blocks in the preconditioner.
std::vector<std::vector<Teko::GO>> construct_block_gids_from_permutation(
    const InversePermutationType &inverse_permutation,
    const std::vector<std::vector<Teko::GO>> &dof_gids);

/// \brief Given a permutation, generate the inverse permutation.
/// \param permutation  [in] Physics-to-sub-block mapping, as one may generate from calling
/// generate_heuristic_permutation above. \return Inverse mapping of sub-blocks to individual
/// physics.
InversePermutationType generate_inverse_permutation(const PermutationType &permutation);

/// \class SubblockParameters
/// \brief Specifies the solvers/preconditioners for each sub-block, as used in
/// generate_parameters_from_permutation. Note that the SubblockParametersBuilder provided below is
/// the preferred way of constructing this object.
struct SubblockParameters {
  //! Solver parameters to used for merged sub-blocks
  Teuchos::RCP<Teuchos::ParameterList> mergedSolverParameters{Teuchos::null};

  //! Name of merged sub-block solver
  std::string mergedSolverName{""};

  //! Preconditioner parameters to used for merged sub-blocks
  Teuchos::RCP<Teuchos::ParameterList> mergedPreconditionerParameters{Teuchos::null};

  //! Name of merged sub-block preconditioner
  std::string mergedPreconditionerName{""};

  //! Individual physics sub-block solver.
  //  If an entry is null, the mergedSolverParameters are used by default.
  //  Note this is only used for sub-blocks that are not merged.
  std::vector<Teuchos::RCP<Teuchos::ParameterList>> subblockSolverParameters{};

  //! Names of individual sub-block solvers
  std::vector<std::string> subblockSolverNames{};

  //! Individual physics sub-block preconditioner.
  //  If an entry is null, the mergedPreconditionerParameters are used by default.
  //  Note this is only used for sub-blocks that are not merged.
  std::vector<Teuchos::RCP<Teuchos::ParameterList>> subblockPreconditionerParameters{};

  //! Names of individual sub-block preconditioners
  std::vector<std::string> subblockPreconditionerNames{};
};

/// \class SubblockParametersBuilder
/// \brief Builder class to construct SubblockParameters.
/// Note that this class supports chaining, e.g.
/// <tt>
///   auto subblockParams = SubblockParametersBuilder{}
///     .set_merged_solver_params(myMergedSolverParams, myMergedSolverName)
///     .set_subblock_solver_params(mySubblockSolverParams, mySubblockSolverNames)
///     .data();
/// </tt>
struct SubblockParametersBuilder {
  /// \brief Set the parameters for the merged sub-block solver.
  ///
  /// \param params [in] Solver parameters to use for merged sub-blocks.
  /// \param solverName  [in] Name to use for merged sub-block solver.
  /// \return Another SubblockParametersBuilder, but with the merged sub-block solver set.
  SubblockParametersBuilder &set_merged_solver_params(Teuchos::RCP<Teuchos::ParameterList> params,
                                                      std::string solverName) {
    data.mergedSolverParameters = params;
    data.mergedSolverName       = solverName;
    return *this;
  };

  /// \brief Set the parameters for the merged sub-block preconditioner.
  ///
  /// \param params [in] Preconditioner parameters to use for merged sub-blocks.
  /// \param preconditionerName  [in] Name to use for merged sub-block preconditioner.
  /// \return Another SubblockParametersBuilder, but with the merged sub-block preconditioner set.
  SubblockParametersBuilder &set_merged_preconditioner_params(
      Teuchos::RCP<Teuchos::ParameterList> params, std::string preconditionerName) {
    data.mergedPreconditionerParameters = params;
    data.mergedPreconditionerName       = preconditionerName;
    return *this;
  };

  /// \brief Set the parameters for each sub-block solver.
  ///  Note that for merged blocks, the merged solver parameters take precedence.
  ///
  /// \param params [in] Solver parameters to use for each sub-blocks, if unmerged.
  /// \param preconditionerName  [in] Name to use for each sub-block solver, if unmerged.
  /// \return Another SubblockParametersBuilder, but with the sub-block solvers set.
  SubblockParametersBuilder &set_subblock_solver_params(
      std::vector<Teuchos::RCP<Teuchos::ParameterList>> params,
      std::vector<std::string> solverNames) {
    data.subblockSolverParameters = params;
    data.subblockSolverNames      = solverNames;
    return *this;
  };

  /// \brief Set the parameters for each sub-block preconditioner.
  ///  Note that for merged blocks, the merged preconditioner parameters take precedence.
  ///
  /// \param params [in] preconditioner parameters to use for each sub-blocks, if unmerged.
  /// \param preconditionerName  [in] Name to use for each sub-block preconditioner, if unmerged.
  /// \return Another SubblockParametersBuilder, but with the sub-block preconditioners set.
  SubblockParametersBuilder &set_subblock_preconditioner_params(
      std::vector<Teuchos::RCP<Teuchos::ParameterList>> params,
      std::vector<std::string> preconditionerNames) {
    data.subblockPreconditionerParameters = params;
    data.subblockPreconditionerNames      = preconditionerNames;
    return *this;
  };

  /// \brief Construct the SubblockParameters.
  /// \return SubblockParameters object to be passed to generate_parameters_from_permutation, for
  /// example.
  SubblockParameters build() const { return data; };

 private:
  SubblockParameters data;
};

/// \brief The default parameters and solver names used in generate_parameters_from_permutation.
/// This is an adaptive sub-block solver that utilizes a schedule of solvers to attempt to solve the
/// sub-block problem. On CPU systems, this is:
///  1. GMRES(30) + Jacobi
///  2. GMRES(30) + Ifpack2 additive Schwarz (overlap=0) with ILU(0)
///  3. GMRES(30) + Ifpack2 additive Schwarz (overlap=1) with ILU(1)
///  4. GMRES(30) + Ifpack2 additive Schwarz (overlap=2) with ILU(2)
///  5. Amesos2/KLU2 direct solver, provided the problem size is less than 100,000 rows.
/// On GPU systems, this is:
///  1. GMRES(50) + Jacobi
///  2. GMRES(50) + Ifpack2 additive Schwarz (overlap=0) with ILU(0)
///  3. GMRES(50) + Ifpack2 additive Schwarz (overlap=1) with ILU(1)
///  4. GMRES(50) + Ifpack2 additive Schwarz (overlap=2) with ILU(2)
///  5. Amesos2/KLU2 direct solver, provided the problem size is less than 100,000 rows.
std::tuple<Teuchos::RCP<Teuchos::ParameterList>, std::string> default_merged_solver_parameters();

/// \brief Given a permutation, generate preconditioner parameters to pass to Teko.
///
/// \param permutation [in] Physics-to-sub-block mapping, as one may generate from calling
/// generate_heuristic_permutation above. \param inverseName [in] Name of the block inverse method.
/// Later, a user may retrieve the InverseFactory associated with the preconditioner via this name.
/// \param heuristicParameters [in] Optional parameters used to construct permutation in
/// generate_heuristic_permutation. If this was provided in generate_heuristic_permutation, the same
/// settings should be provided here. \param subblockParameters  [in] Optional SubblockParameters
/// specifying the solvers/preconditioners to use for each sub-block.
///   For additional information, consult the documentation for SubblockParameters above.
///   By default, the solver specification in described in default_merged_solver_parameters is used.
/// \return Preconditioner parameters as a ParameterList, to be passed to
/// InverseLibrary::buildFromParameterList to construct the InverseLibrary.
Teuchos::RCP<Teuchos::ParameterList> generate_parameters_from_permutation(
    const PermutationType &permutation, std::string inverseName,
    Teuchos::RCP<Teuchos::ParameterList> heuristicParameters = Teuchos::null,
    SubblockParameters subblockParameters                    = {});

namespace LinearOrdering {
/// \brief Compute minimum ordering for block Gauss-Seidel.
/// NOTE: this is an implementation detail. A user should not depend on this.
/// NOTE: This must be called collectively from every rank in the input communicator.
///
/// \param communicator [in] Communicator
/// \param tMaxWalltime [in] Maximum allowable time to run heuristic algorithm.
/// \param blockNorms [in] (Dense) matrix representing multi-physics coupling
/// \param upperTriangular [in] Whether to use upper triangular block Gauss-Seidel
/// \param objectiveFunction  [in] Objective function to minimize.
/// \param out  [in] Optional output stream. If provided, output information from linear ordering
/// solver. \return Pair: first entry is the permutation minimizing the score. Second entry is the
/// score.
std::pair<std::vector<int>, double> compute_min_ordering(
    const Teuchos::RCP<const Teuchos::Comm<int>> &communicator, double tMaxWalltime,
    const BlockNormsViewType &blockNorms, bool upperTriangular,
    std::function<double(const std::vector<int> &)> objectiveFunction,
    Teuchos::RCP<Teuchos::FancyOStream> out);
}  // namespace LinearOrdering

}  // namespace Teko

#endif