// Teko_KrylovSurrogate.hpp — the Belos adaptive hook that drives the
// Krylov-informed block-preconditioner reconfiguration loop.
//
// This header is the orchestration only. The pieces it coordinates live in
// three sibling headers it includes:
//   * Teko_KrylovReducedModel.hpp — the surrogate math (C_hat / b_hat from the
//     Krylov state via per-block Gram-matrix thin SVDs) and the shared type
//     aliases.
//   * Teko_KrylovReconfigIO.hpp   — the JSON file formats (request / reconfig /
//     conv / solved) and the watcher handshake.
//   * Teko_KrylovReconfigPrec.hpp — assembling the reconfigured ("as-if-
//     original") flat blocked system and its preconditioner from an ordering.
//
// Usage (transparent to callers):
//   adaptiveLoop() is registered as the BelosAdaptiveHook and fires after
//   every flexible GMRES solve on a blocked operator (converged or stalled).
//   The request/response/convergence files live in TEKO_RECONFIG_REQUESTS_DIR
//   if set, otherwise kDefaultRequestsDir below.
//
// adaptiveLoop():
//   1. Computes C_hat and writes s<N>_request.json to requests_dir, where
//      N = nextRequestNumber(requests_dir) (0, 1, 2, ... — never reused).
//   2. Polls for s<N>_reconfig.json in the same directory (left in place after
//      being read — its presence is how wait_for_request.py recognizes an
//      already-answered request).
//   3. Optionally sweeps/times each candidate ordering on a fresh flat blocked
//      operator (as if the application had been called with that grouping
//      originally), selecting per selection_mode, and re-solves.
//   4. Overwrites the LHS with the selected solve's result and writes
//      s<N>_conv.json (and s<N>_solved.json for a sweep).
//
// Depends on: Belos, Thyra, Tpetra, Teko, Stratimikos, Teuchos.
// No external JSON library — JSON is written/parsed manually.
#pragma once

// ── Standard library ──────────────────────────────────────────────────────
#include <chrono>
#include <cstdlib>    // getenv
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

// ── Teuchos ───────────────────────────────────────────────────────────────
#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// ── Belos ─────────────────────────────────────────────────────────────────
#include "BelosAdaptiveHook.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosThyraAdapter.hpp"
#include "BelosTypes.hpp"

// ── Thyra ─────────────────────────────────────────────────────────────────
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

// ── Teko ──────────────────────────────────────────────────────────────────
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_FactorTimeRegistry.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_KrylovReconfigIO.hpp"     // JSON file formats + record structs
#include "Teko_KrylovReconfigPrec.hpp"   // reconfigured system + preconditioner
#include "Teko_KrylovReducedModel.hpp"   // type aliases + computeCHat (surrogate math)
#include "Teko_Utilities.hpp"

namespace Teko {
namespace KrylovSurrogate {

// (comm, rank) from block 0 of a product multivector — used to gate the JSON
// file I/O in adaptiveLoop() to rank 0. (getBlockTpetraMV is defined in
// Teko_KrylovReducedModel.hpp.)
inline std::pair<Teuchos::RCP<const Teuchos::Comm<int>>, int>
getCommAndRank(Teuchos::RCP<const MV> mv)
{
    auto comm = getBlockTpetraMV(mv, 0)->getMap()->getComm();
    return {comm, comm->getRank()};
}

// Default location for s<N>_request.json / s<N>_reconfig.json / s<N>_conv.json
// when TEKO_RECONFIG_REQUESTS_DIR is unset. All file types live together in
// this single directory. Hardcoded for now, relative to the Trilinos checkout
// this package lives in.
constexpr const char* kDefaultRequestsDir = "/home/node/codespace/Trilinos/teko-reconfig/requests";

// adaptiveLoop()'s Phase 4 below runs a second flexible FGMRES solve, which
// (if it converges) would re-invoke this same hook recursively on the
// already-reconfigured system — at best redundant work, at worst an unbounded
// reconfigure/solve recursion. g_inAdaptiveLoop / ScopedGuard
// below prevent this: the outermost call sets the flag before Phase 4, so a
// recursive invocation from solver2.solve() sees it set and returns
// immediately at the top of adaptiveLoop(). thread_local since Belos solves
// (and therefore this hook) could in principle run on different threads.
namespace detail {
thread_local bool g_inAdaptiveLoop = false;
}

// RAII guard around g_inAdaptiveLoop, reset on scope exit even if
// solver2.solve() throws (e.g. BlockGmresSolMgrOrthoFailure).
class ScopedAdaptiveLoopGuard {
public:
    ScopedAdaptiveLoopGuard()  { detail::g_inAdaptiveLoop = true; }
    ~ScopedAdaptiveLoopGuard() { detail::g_inAdaptiveLoop = false; }
    ScopedAdaptiveLoopGuard(const ScopedAdaptiveLoopGuard&) = delete;
    ScopedAdaptiveLoopGuard& operator=(const ScopedAdaptiveLoopGuard&) = delete;
};

// Build the as-if-original flat system for `ordering`, solve it, and return
// timing/convergence. Factor and iterate wall times are max-reduced across
// ranks (critical path) so every rank gets identical numbers — which makes
// any time-based selection deterministic across ranks. If writeToLHS is true
// and the solve converged, the solution is unpacked into problem->getLHS().
// The hook-recursion guard must already be held by the caller (this runs a
// flexible solve, which would otherwise re-enter the hook).
//
// Preconditioner construction is guarded: if buildReconfiguredPrec throws on
// ANY rank (e.g. a singular merged block), the failure flag is max-reduced and
// every rank returns a non-converged result without entering the collective
// solve — so one bad candidate in a sweep can't deadlock the run (other ranks
// would otherwise hang in solver.solve()).
inline OrderingResult solveOrdering(
    const std::vector<int>&                            ordering,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>> A_blocked,
    int                                                nb,
    const std::string&                                 method,
    Teuchos::RCP<Teko::InverseFactory>                 invFact,
    Teuchos::RCP<Belos::AdaptiveHook::Problem>         problem,
    Teuchos::RCP<const Teuchos::ParameterList>         orig_params,
    Teuchos::RCP<const Teuchos::Comm<int>>             comm,
    double                                             b_norm,
    bool                                               writeToLHS)
{
    OrderingResult r;
    r.ordering         = ordering;
    r.initial_residual = b_norm;

    ReconfiguredSystem recon;
    int local_err = 0;
    try {
        recon = buildReconfiguredPrec(ordering, A_blocked, nb, method, invFact);
    } catch (const std::exception&) {
        local_err = 1;
    }
    int glob_err = local_err;
    if (comm->getSize() > 1)
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &local_err, &glob_err);
    if (glob_err) {
        // Build failed somewhere; do not enter the collective solve. r is left
        // as not-converged with zero timings.
        return r;
    }

    auto newProblem = Teuchos::rcp(new Problem());
    newProblem->setOperator(Teuchos::rcp_dynamic_cast<const OP>(recon.flatOp));
    newProblem->setRightPrec(Teuchos::rcp_dynamic_cast<const OP>(recon.precOp));

    auto x_new = Thyra::createMembers(recon.flatOp->domain(), 1);
    Thyra::assign(x_new.ptr(), SC(0));
    auto b_new = Thyra::createMembers(recon.flatOp->range(), 1);
    copyOriginalToGrouped(problem->getRHS(), b_new, recon.groups);
    newProblem->setLHS(x_new);
    newProblem->setRHS(b_new);
    newProblem->setProblem();

    auto solverParams = Teuchos::rcp(new Teuchos::ParameterList(*orig_params));
    Belos::BlockGmresSolMgr<SC, MV, OP> solver(newProblem, solverParams);

    const auto t0 = std::chrono::steady_clock::now();
    Belos::ReturnType ret = solver.solve();
    double iterate_sec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    double factor_sec = recon.factor_wall_time_sec;
    if (comm->getSize() > 1) {
        double loc_f = factor_sec, loc_i = iterate_sec;
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &loc_f, &factor_sec);
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &loc_i, &iterate_sec);
    }

    r.iterations            = solver.getNumIters();
    r.converged             = (ret == Belos::Converged);
    r.final_residual        = solver.achievedTol();
    r.factor_wall_time_sec  = factor_sec;
    r.iterate_wall_time_sec = iterate_sec;
    r.total_wall_time_sec   = factor_sec + iterate_sec;

    if (writeToLHS && r.converged)
        copyGroupedToOriginal(x_new, problem->getLHS(), recon.groups);
    return r;
}

// This is the function registered as Belos::AdaptiveHook::HookFn.
// It is called by BelosBlockGmresSolMgr::solve() after a flexible GMRES
// solve converges.
//
// On entry: problem->getLHS() contains the first solve's result.
// On exit:  problem->getLHS() is overwritten with the selected ordering's
//           result (if one was attempted and converged).
// Returns a HookResult: when the selected re-solve converged, it carries that
// solve's iteration count / achieved tolerance so the solver manager can
// report the (possibly rescued) outcome. Fires on both converged and stalled
// first solves — a stall (curDim > 0) takes the full reconfigure path too, so
// a short, deliberately-low-max-iter first solve can be rescued.
inline Belos::AdaptiveHook::HookResult adaptiveLoop(
    const Belos::AdaptiveHook::State&                         state,
    Teuchos::RCP<Belos::AdaptiveHook::Problem>                problem,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>>        A_blocked,
    Teuchos::RCP<const Teuchos::ParameterList>                orig_params,
    const Belos::AdaptiveHook::SolveMetrics&                  solve1_metrics)
{
    // Bail out immediately if this call is the recursive re-entry from
    // Phase 4's solver2.solve() below (see ScopedAdaptiveLoopGuard).
    if (detail::g_inAdaptiveLoop) return {};

    // Opt-in gate. The hook is registered globally at libteko load, so it would
    // otherwise fire on *every* flexible blocked GMRES solve in any process
    // that links Teko — hijacking unrelated solves and blocking them on the
    // reconfig handshake. Do nothing unless TEKO_ADAPTIVE_RECONFIG is set to a
    // truthy value. (A per-solve ParameterList flag is not used because Belos
    // validates the solver parameter list and rejects unknown keys.)
    {
        const char* en = std::getenv("TEKO_ADAPTIVE_RECONFIG");
        const std::string v = en ? std::string(en) : "";
        const bool enabled = !(v.empty() || v == "0" || v == "false" || v == "FALSE");
        if (!enabled) return {};
    }

    // requests_dir comes from TEKO_RECONFIG_REQUESTS_DIR if set, otherwise
    // kDefaultRequestsDir.
    const char* reconfig_dir_env = std::getenv("TEKO_RECONFIG_REQUESTS_DIR");
    const std::string requests_dir = reconfig_dir_env ? std::string(reconfig_dir_env) : kDefaultRequestsDir;

    if (A_blocked.is_null()) {
        std::cerr << "[TekoAdaptive] operator is not blocked; skipping.\n";
        return {};
    }
    if (state.curDim == 0) {
        std::cerr << "[TekoAdaptive] curDim==0; nothing to compute.\n";
        return {};
    }

    // All file I/O (request/reconfig/convergence JSON) is rank-0-only; the
    // numerics (computeCHat, buildReconfiguredPrec, second solve) are
    // collective and run identically on every rank.
    auto [comm, rank] = getCommAndRank(state.V);

    if (rank == 0) fs::create_directories(requests_dir);
    // All request/reconfig/convergence files live together in requests_dir.
    const std::string convergence_dir = requests_dir;

    // Determine block structure from A_blocked.
    const int nb = A_blocked->productRange()->numBlocks();
    std::vector<int> block_sizes(nb);
    for (int j = 0; j < nb; ++j)
        block_sizes[j] =
            static_cast<int>(A_blocked->productRange()->getBlock(j)->dim());

    if (rank == 0)
        std::cout << "[TekoAdaptive] first solve done. curDim=" << state.curDim
                  << "  nb=" << nb << "\n";

    // ||b||: shared initial residual for both solves. Both start from x0 = 0,
    // and reordering (Phase 3) only regroups vector blocks, which preserves
    // the 2-norm.
    double b_norm = 0.0;
    {
        Teuchos::Array<SC> norms(1);
        Thyra::norms_2(*problem->getRHS(), norms());
        b_norm = norms[0];
    }

    // Setup (factorization) time of the first solve: everything the
    // application spent inside Teko::buildInverse / Teko::rebuildInverse since
    // the previous hook invocation — i.e. building the preconditioner that
    // the solve which just converged used. Drain the registry now (Phase 3's
    // buildReconfiguredPrec will bump it again; it is timed separately), and
    // re-drain on every exit path below so solve 2's factorizations never
    // leak into the *next* solve's setup figure. Factorization is rank-local
    // work, so report the max across ranks (the wall-clock critical path).
    double solve1_factor_sec = Teko::FactorTimeRegistry::read();
    Teko::FactorTimeRegistry::reset();
    if (comm->getSize() > 1) {
        double local = solve1_factor_sec;
        Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &local, &solve1_factor_sec);
    }
    struct RegistryResetGuard {
        ~RegistryResetGuard() { Teko::FactorTimeRegistry::reset(); }
    } registry_reset_guard;

    SolveStats s1{solve1_metrics.num_iters, b_norm, solve1_metrics.achieved_tol,
                  solve1_factor_sec, solve1_metrics.wall_time_sec,
                  solve1_factor_sec + solve1_metrics.wall_time_sec};

    // A stalled first solve (hit max iterations / lost accuracy) is NOT a dead
    // end — it is exactly the case reconfiguration should rescue, and the
    // surrogate builds fine from its Krylov data (curDim == max iters). So it
    // falls through to the same surrogate / watcher / sweep / re-solve path as
    // a converged solve; s1 still records the stall (iterations == max iters,
    // time-to-stall), and if the selected re-solve converges this hook reports
    // it back so the manager returns success instead of the original stall.
    if (rank == 0 && !solve1_metrics.converged)
        std::cerr << "[TekoAdaptive] first solve stalled ("
                  << solve1_metrics.num_iters
                  << " iters); attempting reconfiguration rescue.\n";

    // ── Phase 1: compute C_hat, b_hat and write s<N>_request.json ──────────
    // computeCHat is collective (mvTransMv reductions + distributed applies);
    // its result is replicated identically on every rank, so only rank 0
    // writes the request file.
    CHatData chat = computeCHat(state, A_blocked, problem->getRHS(), block_sizes);
    int request_id = -1;
    if (rank == 0) request_id = writeRequestJson(chat, requests_dir);

    // ── Phase 2: wait for s<N>_reconfig.json ───────────────────────────────
    // The response arrives from an external process via the filesystem, so
    // rank 0 polls and the result is broadcast to all ranks.
    ReconfigResponse resp;
    if (rank == 0) resp = waitForReconfig(requests_dir, request_id);

    // Broadcast use_ordering.
    int uo_size = static_cast<int>(resp.use_ordering.size());
    Teuchos::broadcast(*comm, 0, 1, &uo_size);
    resp.use_ordering.resize(uo_size);
    if (uo_size > 0)
        Teuchos::broadcast(*comm, 0, uo_size, resp.use_ordering.data());

    // Broadcast selection_mode as an int code (0 chosen, 1 best_conv, 2 best_time).
    int mode_code = 0;
    if (rank == 0) {
        if      (resp.selection_mode == "best_conv") mode_code = 1;
        else if (resp.selection_mode == "best_time") mode_code = 2;
    }
    Teuchos::broadcast(*comm, 0, 1, &mode_code);
    const char* mode_names[3] = {"chosen", "best_conv", "best_time"};
    const std::string selection_mode = mode_names[mode_code];

    // Broadcast test_orderings (count, then each ordering's length + data).
    int num_tests = static_cast<int>(resp.test_orderings.size());
    Teuchos::broadcast(*comm, 0, 1, &num_tests);
    if (rank != 0) resp.test_orderings.assign(num_tests, std::vector<int>());
    for (int t = 0; t < num_tests; ++t) {
        int len = (rank == 0) ? static_cast<int>(resp.test_orderings[t].size()) : 0;
        Teuchos::broadcast(*comm, 0, 1, &len);
        resp.test_orderings[t].resize(len);
        if (len > 0)
            Teuchos::broadcast(*comm, 0, len, resp.test_orderings[t].data());
    }

    if (resp.use_ordering.empty()) {
        if (rank == 0) {
            std::cerr << "[TekoAdaptive] no ordering received; returning first "
                         "solve result.\n";
            writeConvergenceJson(convergence_dir, request_id, s1, nullptr);
        }
        return {};
    }

    // ── Phase 3: determine method, then sweep / select / re-solve ──────────
    // Determine method from the existing preconditioner type.
    // Heuristic: if it dynamic-casts to BlockLowerTriInverseOp → "gs",
    //            otherwise "jacobi".
    std::string method = "jacobi";
    {
        auto lti = Teuchos::rcp_dynamic_cast<
            const Teko::BlockLowerTriInverseOp>(problem->getRightPrec());
        if (!lti.is_null()) method = "gs";
    }

    // Each candidate solve runs a flexible GMRES that would otherwise re-enter
    // this hook; hold the guard across the whole sweep + final solve.
    ScopedAdaptiveLoopGuard recursion_guard;

    // One Ifpack2 inverse factory, reused for every candidate build below
    // (factory construction is untimed setup; sharing it avoids rebuilding the
    // Stratimikos inverse library once per ordering).
    auto invFact = makeIfpack2InverseFactory();

    // Warm-up before a timed sweep. The first build+solve inside this hook
    // pays one-time costs (kernel first-touch, Stratimikos/Ifpack2 setup
    // paths) that the factorization warm-up does not cover; left unaddressed
    // they would be charged entirely to the first swept ordering and bias
    // best_time. When there is a sweep to time, run one throwaway solve on the
    // identity (all-separate) ordering first so every recorded candidate is
    // measured warm. The result is discarded (writeToLHS = false).
    if (!resp.test_orderings.empty()) {
        std::vector<int> identity(nb);
        std::iota(identity.begin(), identity.end(), 0);
        solveOrdering(identity, A_blocked, nb, method, invFact, problem,
                      orig_params, comm, b_norm, /*writeToLHS=*/false);
    }

    // Sweep: build, solve, and time every candidate ordering on the freshly
    // assembled flat (as-if-original) system, so all candidates are timed
    // comparably to the first solve. Skipped when no test_orderings were given.
    std::vector<OrderingResult> sweep;
    sweep.reserve(resp.test_orderings.size());
    for (const auto& ord : resp.test_orderings)
        sweep.push_back(solveOrdering(ord, A_blocked, nb, method, invFact, problem,
                                      orig_params, comm, b_norm, /*writeToLHS=*/false));

    // Selection. "chosen" returns use_ordering; "best_conv"/"best_time" pick
    // from the sweep. Comparators operate only on cross-rank-reduced,
    // deterministic quantities, so every rank selects the same ordering.
    auto convBetter = [](const OrderingResult& a, const OrderingResult& b) {
        if (a.converged != b.converged) return a.converged;          // converged first
        if (a.iterations != b.iterations) return a.iterations < b.iterations;
        return a.final_residual < b.final_residual;
    };
    auto timeBetter = [&](const OrderingResult& a, const OrderingResult& b) {
        if (a.converged != b.converged) return a.converged;
        if (a.converged && b.converged)
            return a.total_wall_time_sec < b.total_wall_time_sec;
        return convBetter(a, b);  // neither converged: least-bad by convergence
    };

    std::vector<int> selected = resp.use_ordering;
    int selected_index = -1;
    if (!sweep.empty() && mode_code != 0) {
        int best = 0;
        for (int i = 1; i < static_cast<int>(sweep.size()); ++i) {
            const bool better = (mode_code == 1) ? convBetter(sweep[i], sweep[best])
                                                 : timeBetter(sweep[i], sweep[best]);
            if (better) best = i;
        }
        selected = sweep[best].ordering;
        selected_index = best;
    } else {
        // "chosen" (or no candidates): record where use_ordering sits, if present.
        for (int i = 0; i < static_cast<int>(sweep.size()); ++i)
            if (sweep[i].ordering == resp.use_ordering) { selected_index = i; break; }
    }

    if (rank == 0 && !sweep.empty())
        writeSolvedJson(convergence_dir, request_id, selection_mode,
                        sweep, selected_index, selected);

    // ── Phase 4: authoritative solve of the selected ordering ──────────────
    // Re-solve the selected ordering (no cached solution is kept) and write
    // its result into the LHS the application reads back.
    OrderingResult finalRes = solveOrdering(selected, A_blocked, nb, method, invFact, problem,
                                            orig_params, comm, b_norm, /*writeToLHS=*/true);

    SolveStats s2{finalRes.iterations, b_norm, finalRes.final_residual,
                  finalRes.factor_wall_time_sec, finalRes.iterate_wall_time_sec,
                  finalRes.total_wall_time_sec};
    if (rank == 0) writeConvergenceJson(convergence_dir, request_id, s1, &s2);

    if (rank == 0) {
        if (!finalRes.converged)
            std::cerr << "[TekoAdaptive] selected solve did not converge ("
                      << finalRes.iterations << " iters). Using first result.\n";
        else
            std::cout << "[TekoAdaptive] selected solve (" << selection_mode
                      << ") converged (" << finalRes.iterations << " iters).\n";
    }

    // Report the re-solve back to the manager. Only when it converged (and so
    // wrote its result into the LHS) does this override the manager's result —
    // turning a rescued stall into a reported success. If it did not converge,
    // converged stays false and the manager keeps its own (first-solve) result.
    Belos::AdaptiveHook::HookResult result;
    result.converged    = finalRes.converged;
    result.num_iters    = finalRes.iterations;
    result.achieved_tol = finalRes.final_residual;
    return result;
}

} // namespace KrylovSurrogate
} // namespace Teko
