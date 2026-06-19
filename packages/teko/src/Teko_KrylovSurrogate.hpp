// Teko_KrylovSurrogate.hpp — Krylov-informed surrogate for adaptive block
// preconditioning.
//
// After a flexible GMRES solve the Krylov state (V, Z) allows building a
// cheap reduced operator C_hat (equation 3.11 in teko_equations.tex):
//
//   C_hat_{ij} = Q_{i,ri}^T  A_{ij}  Z_j  W_{j,rj}  Sigma_{j,rj}^{-1}
//
// where Q_j / W_j / Sigma_j come from the thin SVD of V_j (the j-th block
// row of the Krylov basis V).  Diagonal blocks are identity by construction.
//
// Usage (transparent to callers):
//   adaptiveLoop() is registered as the BelosAdaptiveHook and fires after
//   every converged flexible GMRES solve on a blocked operator. The
//   request/response/convergence files live in TEKO_RECONFIG_REQUESTS_DIR if
//   set, otherwise kDefaultRequestsDir below.
//
// adaptiveLoop():
//   1. Computes C_hat and writes s<N>_request.json to requests_dir, where
//      N = nextRequestNumber(requests_dir) (0, 1, 2, ... — never reused).
//   2. Polls for s<N>_reconfig.json in the same directory. This file is left
//      in place after being read (not deleted) — its presence is how
//      wait_for_request.py recognizes an already-answered request.
//   3. Reads the ordering vector and assembles the reconfigured system as a
//      fresh flat blocked operator (merged groups become single monolithic
//      blocks), exactly as if the application had been called with that
//      grouping originally — so the second solve's factor/iterate times are
//      directly comparable to the first solve's.
//   4. Runs a second FGMRES solve and overwrites the LHS with the result.
//   5. Writes s<N>_conv.json (in requests_dir, alongside the request and
//      reconfig files) with iteration count, initial/final residual, and wall
//      time for both solves (solve2 is null if step 2 timed out).
//
// Depends on: Belos, Thyra, Tpetra, Teko, Stratimikos, Teuchos.
// No external JSON library — JSON is written/parsed manually.
#pragma once

// ── Standard library ──────────────────────────────────────────────────────
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdlib>    // getenv
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Filesystem (C++17)
#include <filesystem>
namespace fs = std::filesystem;

// ── Teuchos ───────────────────────────────────────────────────────────────
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

// ── Belos ─────────────────────────────────────────────────────────────────
#include "BelosAdaptiveHook.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosThyraAdapter.hpp"
#include "BelosTypes.hpp"

// ── Thyra ─────────────────────────────────────────────────────────────────
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearOpBase.hpp"

// ── Thyra-Tpetra adapter ──────────────────────────────────────────────────
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraMultiVector.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

// ── Tpetra ────────────────────────────────────────────────────────────────
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"

// ── Teko ──────────────────────────────────────────────────────────────────
#include "Teko_BlockDiagonalInverseOp.hpp"
#include "Teko_BlockLowerTriInverseOp.hpp"
#include "Teko_FactorTimeRegistry.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_Utilities.hpp"

// ── Stratimikos ───────────────────────────────────────────────────────────
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

namespace Teko {
namespace KrylovSurrogate {

// ── Type aliases matching teko_ext.cpp ────────────────────────────────────
using SC   = double;
using LO   = int;
using GO   = long long;
using Node = Tpetra::KokkosCompat::KokkosOpenMPWrapperNode;

using SDM      = Teuchos::SerialDenseMatrix<int, SC>;
using TpetraMV = Tpetra::MultiVector<SC, LO, GO, Node>;
using MV       = Thyra::MultiVectorBase<SC>;
using OP       = Thyra::LinearOpBase<SC>;
using Problem  = Belos::LinearProblem<SC, MV, OP>;
using State    = Belos::GmresIterationState<SC, MV>;

// ═══════════════════════════════════════════════════════════════════════════
// Section 1: distributed helpers
// ═══════════════════════════════════════════════════════════════════════════

// Return block j of a Thyra ProductMultiVector as its underlying (possibly
// distributed) Tpetra::MultiVector. No data is copied or gathered.
inline Teuchos::RCP<const TpetraMV> getBlockTpetraMV(
    Teuchos::RCP<const MV> mv,
    int j)
{
    auto pmv  = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC>>(mv);
    TEUCHOS_TEST_FOR_EXCEPTION(pmv.is_null(), std::runtime_error,
        "Teko::KrylovSurrogate::getBlockTpetraMV: mv is not a ProductMultiVectorBase");

    auto blk  = pmv->getMultiVectorBlock(j);
    auto tmv  = Teuchos::rcp_dynamic_cast<const Thyra::TpetraMultiVector<SC, LO, GO, Node>>(blk);
    TEUCHOS_TEST_FOR_EXCEPTION(tmv.is_null(), std::runtime_error,
        "Teko::KrylovSurrogate::getBlockTpetraMV: block " + std::to_string(j) +
        " is not a TpetraMultiVector");

    return tmv->getConstTpetraMultiVector();
}

// (comm, rank) from block 0 of a product multivector — used to gate the JSON
// file I/O in adaptiveLoop() to rank 0.
inline std::pair<Teuchos::RCP<const Teuchos::Comm<int>>, int>
getCommAndRank(Teuchos::RCP<const MV> mv)
{
    auto comm = getBlockTpetraMV(mv, 0)->getMap()->getComm();
    return {comm, comm->getRank()};
}

// C := A^T B as a replicated SerialDenseMatrix (numVecs(A) × numVecs(B)).
// A and B may be distributed arbitrarily (they must share a row distribution).
// Internally Tpetra computes the local GEMM and sum-reduces the small result
// across ranks (Tpetra::MultiVector::multiply with a locally-replicated
// target), so every rank gets the bit-identical matrix — no broadcasts needed
// downstream. This is the only cross-rank reduction in computeCHat(); nothing
// problem-sized is ever gathered.
inline SDM mvTransMv(const TpetraMV& A, const TpetraMV& B)
{
    const int ka = static_cast<int>(A.getNumVectors());
    const int kb = static_cast<int>(B.getNumVectors());

    auto repMap = Teuchos::rcp(new Tpetra::Map<LO, GO, Node>(
        static_cast<Tpetra::global_size_t>(ka), 0,
        A.getMap()->getComm(), Tpetra::LocallyReplicated));
    TpetraMV C(repMap, kb);
    C.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, A, B, 0.0);

    SDM out(ka, kb);
    for (int c = 0; c < kb; ++c) {
        auto view = C.getVector(c)->getData();
        for (int r = 0; r < ka; ++r)
            out(r, c) = view[r];
    }
    return out;
}

// ═══════════════════════════════════════════════════════════════════════════
// Section 2: CHat data structure and computation
// ═══════════════════════════════════════════════════════════════════════════

struct CHatData {
    int                nb;          // number of blocks
    int                krylov_dim;  // curDim used
    std::vector<int>   block_sizes;
    std::vector<int>   ranks;       // r_j per block

    // blocks[i][j] is (ranks[i] × ranks[j]) — reduced operator
    std::vector<std::vector<SDM>>  blocks;
    // b_hat[j] is (ranks[j] × 1) = Q[j]^T * b_j — reduced right-hand side
    std::vector<SDM>               b_hat;
};

// Compute C_hat from the final FGMRES state.
//   state      — from BlockFGmresIter::getState() after convergence
//   A_blocked  — the blocked system operator
//   b          — the right-hand side of the linear system (for b_hat = Q^T b)
//   block_sizes — sizes of each block
//   svd_tol    — truncation threshold relative to leading singular value
//
// Multi-rank (Gram-matrix) formulation. The mathematical definitions are
// unchanged from the GESVD version:
//
//   V_j = Q_j Sigma_j W_j^T  (thin SVD),
//   C_hat_{ij} = Q_{i,ri}^T A_{ij} Z_j W_{j,rj} Sigma_{j,rj}^{-1},
//   b_hat_j    = Q_{j,rj}^T b_j,
//
// but Q_j (problem-sized, distributed) is eliminated via
// Q_j = V_j W_j Sigma_j^{-1}, leaving only small replicated quantities:
//
//   G_j   := V_j^T V_j = W_j Sigma_j^2 W_j^T   (curDim × curDim, SPD)
//     → SYEV eigendecomposition gives W_j and sigma_{j,k} = sqrt(lambda_{j,k})
//   b_hat_j    = Sigma_{j,rj}^{-1} W_{j,rj}^T (V_j^T b_j)
//   C_hat_{ij} = Sigma_{i,ri}^{-1} W_{i,ri}^T (V_i^T A_{ij} Z_j) W_{j,rj} Sigma_{j,rj}^{-1}
//
// Every V^T(...) product is a mvTransMv() reduction (local GEMM + Allreduce of
// a curDim-sized matrix); the sparse A_{ij} Z_j apply stays fully distributed.
// All ranks therefore compute identical C_hat/b_hat with no gathers or
// broadcasts. This is a collective: every rank must call it.
//
// Numerical caveat: truncating sigma_k >= svd_tol * sigma_max on G_j's
// eigenvalues means lambda_k >= svd_tol^2 * lambda_max (= 1e-16*lambda_max at
// the default), right at double roundoff. For ill-conditioned V_j
// (cond >~ 1e4) ranks[j] may differ from the GESVD result. SYEV's eigenvector
// sign convention may also flip individual rows/columns of C_hat/b_hat
// relative to GESVD — magnitudes and norms are unaffected.
inline CHatData computeCHat(
    const State&                                       state,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>> A_blocked,
    Teuchos::RCP<const MV>                             b,
    const std::vector<int>&                            block_sizes,
    double                                             svd_tol = 1e-8)
{
    const int nb       = static_cast<int>(block_sizes.size());
    const int curDim   = state.curDim;

    CHatData chat;
    chat.nb          = nb;
    chat.krylov_dim  = curDim;
    chat.block_sizes = block_sizes;
    chat.ranks.resize(nb, 0);
    chat.b_hat.resize(nb);
    chat.blocks.assign(nb, std::vector<SDM>(nb));

    Teuchos::LAPACK<int, SC> lapack;

    const Teuchos::Range1D krylovCols(0, curDim - 1);

    // Per-block distributed views (first curDim columns; no copies).
    std::vector<Teuchos::RCP<const TpetraMV>> Vt(nb), Zt(nb);
    for (int j = 0; j < nb; ++j) {
        Vt[j] = getBlockTpetraMV(state.V, j)->subView(krylovCols);
        Zt[j] = getBlockTpetraMV(state.Z, j)->subView(krylovCols);
    }

    // ── Per-block: eigendecomposition of the Gram matrix G_j = V_j^T V_j ──
    std::vector<SDM> WTrunc(nb);            // WTrunc[j] = curDim × r_j
    std::vector<std::vector<SC>> Sigma(nb); // Sigma[j]  = leading r_j sing. values

    for (int j = 0; j < nb; ++j) {
        // G_j (curDim × curDim, replicated identically on every rank).
        SDM Gj = mvTransMv(*Vt[j], *Vt[j]);

        // SYEV: eigenvalues ascending; eigenvectors overwrite Gj's columns.
        std::vector<SC> eig(curDim);
        int lwork = -1, info = 0;
        SC work_query;
        lapack.SYEV('V', 'U', curDim, Gj.values(), Gj.stride(),
                    eig.data(), &work_query, lwork, &info);
        lwork = static_cast<int>(work_query);
        std::vector<SC> work(lwork);
        lapack.SYEV('V', 'U', curDim, Gj.values(), Gj.stride(),
                    eig.data(), work.data(), lwork, &info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
            "Teko::KrylovSurrogate::computeCHat: SYEV failed (info=" +
            std::to_string(info) + ") for block " + std::to_string(j));

        // Truncate against the leading eigenvalue. sigma_k >= svd_tol*sigma_max
        // on singular values is lambda_k >= svd_tol^2*lambda_max on G_j's
        // eigenvalues (lambda_k = sigma_k^2). The threshold is floored at a
        // few ulps of lambda_max: eigenvalues of G_j below that are roundoff
        // noise (directions with sigma <~ sqrt(eps)*sigma_max cannot be
        // resolved through the Gram matrix), and keeping one would inject an
        // O(1/sigma) garbage row/column into C_hat.
        const SC eps        = std::numeric_limits<SC>::epsilon();
        const SC lambda_max = eig[curDim - 1];
        const SC threshold  = (lambda_max > 0.0
            ? std::max(svd_tol * svd_tol, 8.0 * eps) * lambda_max : 0.0);
        int rj = 0;
        while (rj < curDim && eig[curDim - 1 - rj] >= threshold) ++rj;
        if (rj == 0) rj = 1;  // keep at least one direction

        chat.ranks[j] = rj;

        // Descending order: r-th largest eigenpair is SYEV column curDim-1-r.
        Sigma[j].resize(rj);
        WTrunc[j].shape(curDim, rj);
        for (int r = 0; r < rj; ++r) {
            const int src = curDim - 1 - r;
            Sigma[j][r] = std::sqrt(std::max(eig[src], 0.0));
            for (int kd = 0; kd < curDim; ++kd)
                WTrunc[j](kd, r) = Gj(kd, src);
        }

        // b_hat[j] = Sigma_j^{-1} W_j^T (V_j^T b_j)   (rj × 1)
        SDM Vtb = mvTransMv(*Vt[j], *getBlockTpetraMV(b, j));   // curDim × 1
        chat.b_hat[j].shape(rj, 1);
        chat.b_hat[j].multiply(Teuchos::TRANS, Teuchos::NO_TRANS,
                               1.0, WTrunc[j], Vtb, 0.0);
        for (int r = 0; r < rj; ++r) {
            const SC inv_s = (std::abs(Sigma[j][r]) > 1e-14)
                             ? 1.0 / Sigma[j][r] : 0.0;
            chat.b_hat[j](r, 0) *= inv_s;
        }
    }

    // ── Per (i,j) block of C_hat ──────────────────────────────────────────
    // C_hat[i][j] = Sigma_i^{-1} W_i^T (V_i^T A_ij Z_j) W_j Sigma_j^{-1}.
    // A_ij Z_j is the existing distributed sparse apply; V_i^T(...) is a
    // mvTransMv reduction.

    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nb; ++j) {
            const int ri = chat.ranks[i];
            const int rj = chat.ranks[j];
            chat.blocks[i][j].shape(ri, rj);

            if (i == j) {
                // Diagonal block is identity by construction.
                for (int d = 0; d < std::min(ri, rj); ++d)
                    chat.blocks[i][j](d, d) = 1.0;
                continue;
            }

            auto A_ij = A_blocked->getBlock(i, j);
            if (A_ij.is_null()) {
                // Zero operator block → C_hat[i][j] = 0 (already zero-initialized)
                continue;
            }

            // Extract A_ij as Tpetra CrsMatrix for direct apply.
            auto A_ij_tpetra_lo =
                Teuchos::rcp_dynamic_cast<const Thyra::TpetraLinearOp<SC, LO, GO, Node>>(A_ij);
            if (A_ij_tpetra_lo.is_null()) {
                // Fall back: try applying via Thyra (slower but general)
                // For now, just zero this block and continue.
                continue;
            }
            auto A_crs = Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<SC, LO, GO, Node>>(
                A_ij_tpetra_lo->getConstTpetraOperator());
            if (A_crs.is_null()) continue;

            // AijZj = A_ij * Z_j — fully distributed sparse apply, directly on
            // block j's Tpetra view (no extract/repack).
            auto AijZj_tmv = Teuchos::rcp(new TpetraMV(A_crs->getRangeMap(), curDim));
            A_crs->apply(*Zt[j], *AijZj_tmv, Teuchos::NO_TRANS, 1.0, 0.0);

            // M_ij = V_i^T (A_ij Z_j)   (curDim × curDim, replicated)
            SDM Mij = mvTransMv(*Vt[i], *AijZj_tmv);

            // T = M_ij * W_j  (curDim × rj), columns scaled by 1/sigma_{j,c}
            SDM T(curDim, rj);
            T.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                       1.0, Mij, WTrunc[j], 0.0);
            for (int c = 0; c < rj; ++c) {
                const SC inv_s = (std::abs(Sigma[j][c]) > 1e-14)
                                 ? 1.0 / Sigma[j][c] : 0.0;
                for (int r = 0; r < curDim; ++r)
                    T(r, c) *= inv_s;
            }

            // C_hat[i][j] = W_i^T * T  (ri × rj), rows scaled by 1/sigma_{i,r}
            chat.blocks[i][j].multiply(Teuchos::TRANS, Teuchos::NO_TRANS,
                                       1.0, WTrunc[i], T, 0.0);
            for (int r = 0; r < ri; ++r) {
                const SC inv_s = (std::abs(Sigma[i][r]) > 1e-14)
                                 ? 1.0 / Sigma[i][r] : 0.0;
                for (int c = 0; c < rj; ++c)
                    chat.blocks[i][j](r, c) *= inv_s;
            }
        }
    }

    return chat;
}

// ═══════════════════════════════════════════════════════════════════════════
// Section 3: JSON I/O
// ═══════════════════════════════════════════════════════════════════════════

// Write a dense matrix as a JSON array-of-arrays.
inline void writeDenseJson(std::ostream& os, const SDM& M, int indent)
{
    const std::string pad(indent, ' ');
    const std::string inner(indent + 2, ' ');
    os << pad << "[\n";
    for (int r = 0; r < M.numRows(); ++r) {
        os << inner << "[";
        for (int c = 0; c < M.numCols(); ++c) {
            if (c) os << ", ";
            os << std::setprecision(17) << M(r, c);
        }
        os << "]";
        if (r + 1 < M.numRows()) os << ",";
        os << "\n";
    }
    os << pad << "]";
}

// Write a (k × 1) column vector as a flat JSON array.
inline void writeVectorJson(std::ostream& os, const SDM& v)
{
    os << "[";
    for (int r = 0; r < v.numRows(); ++r) {
        if (r) os << ", ";
        os << std::setprecision(17) << v(r, 0);
    }
    os << "]";
}

// Scan requests_dir for existing s<N>_request.json files and return
// max(N) + 1, or 0 if none exist. This is how the C++ side picks the id for
// the next request, so request ids are unique and monotonically increasing
// across the lifetime of requests_dir (they are never reused, even if old
// s<N>_request.json files are deleted).
inline int nextRequestNumber(const std::string& requests_dir)
{
    int next = 0;
    if (!fs::exists(requests_dir)) return next;

    // Consider both s<N>_request.json and s<N>_conv.json. A non-converged
    // (record-only) solve writes a s<N>_conv.json with no matching
    // request.json; counting conv files too keeps that N reserved so a later
    // request can't reuse it and overwrite the record.
    const std::string prefix = "s";
    const std::array<std::string, 3> suffixes{"_request.json", "_conv.json", "_solved.json"};
    for (const auto& entry : fs::directory_iterator(requests_dir)) {
        if (!entry.is_regular_file()) continue;
        const std::string name = entry.path().filename().string();
        if (name.compare(0, prefix.size(), prefix) != 0) continue;
        for (const auto& suffix : suffixes) {
            if (name.size() <= prefix.size() + suffix.size()) continue;
            if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0) continue;

            const std::string digits = name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
            if (digits.empty() || !std::all_of(digits.begin(), digits.end(),
                                                [](char c){ return std::isdigit(static_cast<unsigned char>(c)); })) continue;

            const int n = std::stoi(digits);
            if (n + 1 > next) next = n + 1;
        }
    }
    return next;
}

// equation_ends[k] = sum(ranks[0..k-1]) for k = 0..nb, so equation_ends[0]==0
// and equation_ends[nb] == R (the total rank, R = sum(ranks)). Each
// equation_ends[k+1] is one-past-the-last row/col index of block k in the
// assembled C_hat matrix / b_hat vector below.
inline std::vector<int> computeEquationEnds(const CHatData& chat)
{
    std::vector<int> ends(chat.nb + 1, 0);
    for (int j = 0; j < chat.nb; ++j)
        ends[j + 1] = ends[j] + chat.ranks[j];
    return ends;
}

// Assemble the full R x R C_hat matrix (R = sum(ranks)) from chat.blocks,
// placing block (i,j) — size ranks[i] x ranks[j] — at rows
// [equation_ends[i], equation_ends[i+1]) and cols
// [equation_ends[j], equation_ends[j+1]).
inline SDM assembleCHat(const CHatData& chat, const std::vector<int>& equation_ends)
{
    const int R = equation_ends.back();
    SDM C(R, R);
    for (int i = 0; i < chat.nb; ++i) {
        for (int j = 0; j < chat.nb; ++j) {
            const SDM& blk = chat.blocks[i][j];
            for (int r = 0; r < blk.numRows(); ++r)
                for (int c = 0; c < blk.numCols(); ++c)
                    C(equation_ends[i] + r, equation_ends[j] + c) = blk(r, c);
        }
    }
    return C;
}

// Assemble the full R x 1 b_hat vector (R = sum(ranks)) by stacking
// chat.b_hat[j] — size ranks[j] x 1 — at rows
// [equation_ends[j], equation_ends[j+1]).
inline SDM assembleBHat(const CHatData& chat, const std::vector<int>& equation_ends)
{
    const int R = equation_ends.back();
    SDM b(R, 1);
    for (int j = 0; j < chat.nb; ++j) {
        const SDM& bj = chat.b_hat[j];
        for (int r = 0; r < bj.numRows(); ++r)
            b(equation_ends[j] + r, 0) = bj(r, 0);
    }
    return b;
}

// Write s<N>_request.json to requests_dir, where N = nextRequestNumber(requests_dir).
// Returns N so the caller can wait for the matching s<N>_reconfig.json
// and tag s<N>_conv.json with the same id.
//
// Format:
// {
//   "n_blocks": nb,
//   "block_sizes": [...],
//   "krylov_dim": m,
//   "ranks": [...],
//   "equation_ends": [0, r0, r0+r1, ..., R],   // length nb+1, R = sum(ranks)
//   "C_hat": [ [ ... ], ... ],                  // single R x R matrix
//   "b_hat": [ ... ]                            // single length-R vector
// }
//
// C_hat is the full R x R surrogate operator (R = sum(ranks)) and b_hat the
// full length-R reduced RHS, both assembled block-by-block. "equation_ends"
// gives the block boundaries within them: equation_ends[k] is the first
// row/col index of block k and equation_ends[k+1] one past its last, so
// block k occupies [equation_ends[k], equation_ends[k+1]).
//
// Every entry scales with nb (number of blocks) and ranks[j] (Krylov-derived
// rank, <= krylov_dim) only — nothing here scales with the global problem
// size n, so this stays small even for very large systems.
inline int writeRequestJson(const CHatData& chat, const std::string& requests_dir)
{
    const int request_id = nextRequestNumber(requests_dir);

    const std::vector<int> equation_ends = computeEquationEnds(chat);
    const SDM C = assembleCHat(chat, equation_ends);
    const SDM b = assembleBHat(chat, equation_ends);

    // Write to a temp file and atomically rename into place. Without this, a
    // reader's fs::exists(path) check (in wait_for_request.py) can observe
    // the file mid-write and parse incomplete JSON. A same-directory
    // rename() is atomic on POSIX, so readers only ever see the file fully
    // written or not at all.
    const std::string path     = requests_dir + "/s" + std::to_string(request_id) + "_request.json";
    const std::string tmp_path = path + ".tmp";
    std::ofstream ofs(tmp_path);
    TEUCHOS_TEST_FOR_EXCEPTION(!ofs.is_open(), std::runtime_error,
        "Teko::KrylovSurrogate::writeRequestJson: cannot open " + tmp_path);

    ofs << "{\n";
    ofs << "  \"n_blocks\": " << chat.nb << ",\n";

    ofs << "  \"block_sizes\": [";
    for (int j = 0; j < chat.nb; ++j) {
        if (j) ofs << ", ";
        ofs << chat.block_sizes[j];
    }
    ofs << "],\n";

    ofs << "  \"krylov_dim\": " << chat.krylov_dim << ",\n";

    ofs << "  \"ranks\": [";
    for (int j = 0; j < chat.nb; ++j) {
        if (j) ofs << ", ";
        ofs << chat.ranks[j];
    }
    ofs << "],\n";

    ofs << "  \"equation_ends\": [";
    for (size_t k = 0; k < equation_ends.size(); ++k) {
        if (k) ofs << ", ";
        ofs << equation_ends[k];
    }
    ofs << "],\n";

    // C_hat: single R x R matrix
    ofs << "  \"C_hat\":\n";
    writeDenseJson(ofs, C, 2);
    ofs << ",\n";

    // b_hat: single length-R vector
    ofs << "  \"b_hat\": ";
    writeVectorJson(ofs, b);
    ofs << "\n";

    ofs << "}\n";
    ofs.close();

    fs::rename(tmp_path, path);

    std::cout << "[TekoAdaptive] wrote " << path << "\n";
    return request_id;
}

// ── Minimal strict parsing for the reconfig response ──────────────────────
// s<N>_reconfig.json is machine-generated by wait_for_request.py in a fixed
// shape, e.g.
//   {"selection_mode": "chosen",
//    "use_ordering": [0,0,1],
//    "opt_ordering": [0,0,1],
//    "exh_opt_ordering": [],
//    "test_orderings": [[0,1,2],[0,0,1],[0,0,0]]}
// These helpers extract exactly the fields the C++ side consumes — a string
// enum, an int array, and an array-of-int-arrays — keyed off field names and
// rejecting malformed numeric tokens. It is deliberately a small targeted
// parser, not a general JSON reader, but keying off the field names makes it
// tolerant of extra fields in any position.

// Parse one integer from a token, rejecting anything std::stoi would silently
// truncate (e.g. "0.9" -> 0).
inline int parseStrictInt(std::string tok, const std::string& ctx)
{
    tok.erase(std::remove_if(tok.begin(), tok.end(),
                             [](char c){ return std::isspace(static_cast<unsigned char>(c)); }),
              tok.end());
    std::size_t consumed = 0;
    int value = 0;
    try {
        value = std::stoi(tok, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error(ctx + ": non-integer entry \"" + tok + "\"");
    }
    if (consumed != tok.size())
        throw std::runtime_error(ctx + ": non-integer entry \"" + tok + "\"");
    return value;
}

// Index of the ']' matching the '[' at openPos (balanced). npos if unbalanced.
inline std::size_t matchBracket(const std::string& s, std::size_t openPos)
{
    int depth = 0;
    for (std::size_t i = openPos; i < s.size(); ++i) {
        if (s[i] == '[') ++depth;
        else if (s[i] == ']') { if (--depth == 0) return i; }
    }
    return std::string::npos;
}

// Parse "a, b, c" (comma-separated, no brackets) into ints.
inline std::vector<int> parseIntCsv(const std::string& body, const std::string& ctx)
{
    std::vector<int> out;
    std::istringstream ss(body);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        const bool nonspace = std::any_of(tok.begin(), tok.end(),
            [](char c){ return !std::isspace(static_cast<unsigned char>(c)); });
        if (nonspace) out.push_back(parseStrictInt(tok, ctx));
    }
    return out;
}

// Optional "key": "value" string field. Returns false if the key is absent.
inline bool parseStringField(const std::string& content, const std::string& key,
                             std::string& out)
{
    const auto kp = content.find("\"" + key + "\"");
    if (kp == std::string::npos) return false;
    const auto colon = content.find(':', kp + key.size() + 2);
    if (colon == std::string::npos) return false;
    const auto q1 = content.find('"', colon);
    if (q1 == std::string::npos) return false;
    const auto q2 = content.find('"', q1 + 1);
    if (q2 == std::string::npos) return false;
    out = content.substr(q1 + 1, q2 - q1 - 1);
    return true;
}

// Optional "key": [ints] field. Returns false if the key is absent.
inline bool parseIntArrayField(const std::string& content, const std::string& key,
                               std::vector<int>& out, const std::string& ctx)
{
    const auto kp = content.find("\"" + key + "\"");
    if (kp == std::string::npos) return false;
    const auto lb = content.find('[', kp);
    if (lb == std::string::npos)
        throw std::runtime_error(ctx + ": \"" + key + "\" has no array");
    const auto rb = matchBracket(content, lb);
    if (rb == std::string::npos)
        throw std::runtime_error(ctx + ": \"" + key + "\" array not closed");
    out = parseIntCsv(content.substr(lb + 1, rb - lb - 1), ctx + " \"" + key + "\"");
    return true;
}

// Optional "key": [[ints],[ints],...] field. Returns false if the key is absent.
inline bool parseListOfIntArraysField(const std::string& content, const std::string& key,
                                      std::vector<std::vector<int>>& out,
                                      const std::string& ctx)
{
    const auto kp = content.find("\"" + key + "\"");
    if (kp == std::string::npos) return false;
    const auto outerL = content.find('[', kp);
    if (outerL == std::string::npos)
        throw std::runtime_error(ctx + ": \"" + key + "\" has no array");
    const auto outerR = matchBracket(content, outerL);
    if (outerR == std::string::npos)
        throw std::runtime_error(ctx + ": \"" + key + "\" array not closed");
    out.clear();
    std::size_t i = outerL + 1;
    while (true) {
        const auto innerL = content.find('[', i);
        if (innerL == std::string::npos || innerL > outerR) break;
        const auto innerR = matchBracket(content, innerL);
        if (innerR == std::string::npos || innerR > outerR)
            throw std::runtime_error(ctx + ": \"" + key + "\" inner array not closed");
        out.push_back(parseIntCsv(content.substr(innerL + 1, innerR - innerL - 1),
                                  ctx + " \"" + key + "\""));
        i = innerR + 1;
    }
    return true;
}

// Parsed reconfig response. use_ordering empty => no response (timeout).
struct ReconfigResponse {
    std::string                   selection_mode;  // "chosen" | "best_conv" | "best_time"
    std::vector<int>              use_ordering;     // ordering applied in "chosen" mode
    std::vector<std::vector<int>> test_orderings;   // candidates to sweep/time (may be empty)
};

// Poll for s<request_id>_reconfig.json and parse it. The file is left in place
// after being read (NOT deleted) — its presence is how wait_for_request.py
// recognizes an already-answered request. Returns a response with empty
// use_ordering on timeout.
inline ReconfigResponse waitForReconfig(
    const std::string& requests_dir,
    int                request_id,
    int                timeout_s = 600)
{
    const std::string path = requests_dir + "/s" + std::to_string(request_id) + "_reconfig.json";
    const auto deadline =
        std::chrono::steady_clock::now() + std::chrono::seconds(timeout_s);

    std::cout << "[TekoAdaptive] waiting for " << path << " ...\n";
    std::cout.flush();

    while (std::chrono::steady_clock::now() < deadline) {
        if (fs::exists(path)) {
            std::ifstream ifs(path);
            std::string content((std::istreambuf_iterator<char>(ifs)),
                                 std::istreambuf_iterator<char>());
            ifs.close();

            const std::string ctx = "[TekoAdaptive] " + path;
            ReconfigResponse r;
            if (!parseIntArrayField(content, "use_ordering", r.use_ordering, ctx))
                throw std::runtime_error(ctx + ": no \"use_ordering\" field found");
            if (!parseStringField(content, "selection_mode", r.selection_mode))
                r.selection_mode = "chosen";
            parseListOfIntArraysField(content, "test_orderings", r.test_orderings, ctx);

            std::cout << "[TekoAdaptive] reconfig: mode=" << r.selection_mode
                      << "  use_ordering=[";
            for (int k = 0; k < (int)r.use_ordering.size(); ++k) {
                if (k) std::cout << ", ";
                std::cout << r.use_ordering[k];
            }
            std::cout << "]  test_orderings=" << r.test_orderings.size() << "\n";
            return r;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }

    std::cerr << "[TekoAdaptive] timed out waiting for " << path << "\n";
    return {};
}

// ═══════════════════════════════════════════════════════════════════════════
// Section 4: Preconditioner rebuild from ordering
// ═══════════════════════════════════════════════════════════════════════════

// Extract A_blocked(i,j) as a Tpetra CrsMatrix (null if the block is null or
// not a Tpetra CRS operator).
inline Teuchos::RCP<const Tpetra::CrsMatrix<SC, LO, GO, Node>> getBlockCrs(
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>> A_blocked, int i, int j)
{
    auto blk = A_blocked->getBlock(i, j);
    if (blk.is_null()) return Teuchos::null;
    auto lo = Teuchos::rcp_dynamic_cast<
        const Thyra::TpetraLinearOp<SC, LO, GO, Node>>(blk);
    if (lo.is_null()) return Teuchos::null;
    return Teuchos::rcp_dynamic_cast<const Tpetra::CrsMatrix<SC, LO, GO, Node>>(
        lo->getConstTpetraOperator());
}

// Row map of original block k: row map of the first non-null Tpetra block in
// block-row k. This is the distribution of block k's vector entries.
inline Teuchos::RCP<const Tpetra::Map<LO, GO, Node>> blockRowMap(
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>> A_blocked, int k, int nb)
{
    for (int j = 0; j < nb; ++j) {
        auto crs = getBlockCrs(A_blocked, k, j);
        if (!crs.is_null()) return crs->getRowMap();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error,
        "Teko::KrylovSurrogate::blockRowMap: block-row " + std::to_string(k) +
        " has no Tpetra blocks.");
    return Teuchos::null;  // unreachable
}

// Map of a group: each rank's GIDs are the members' local rows (as global
// indices within the group, i.e. shifted by the member's within-group offset)
// concatenated in member order. The vector layout this induces is exactly
// "member-local parts back to back", which is what the local repacking in
// copyOriginalToGrouped/copyGroupedToOriginal assumes. For a singleton group
// this is just the member's own row map.
inline Teuchos::RCP<const Tpetra::Map<LO, GO, Node>> buildGroupMap(
    const std::vector<int>&                                           members,
    const std::vector<Teuchos::RCP<const Tpetra::Map<LO, GO, Node>>>& rowMapOf,
    const std::vector<int>&                                           groupOffsets)
{
    using TpMap = Tpetra::Map<LO, GO, Node>;
    if (members.size() == 1) return rowMapOf[members[0]];

    std::vector<GO> myGIDs;
    for (size_t s = 0; s < members.size(); ++s) {
        const auto&  mmap = rowMapOf[members[s]];
        const size_t nloc = mmap->getLocalNumElements();
        for (size_t lr = 0; lr < nloc; ++lr)
            myGIDs.push_back(static_cast<GO>(groupOffsets[s]) +
                             mmap->getGlobalElement(static_cast<LO>(lr)));
    }
    return Teuchos::rcp(new TpMap(
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
        Teuchos::ArrayView<const GO>(myGIDs.data(), myGIDs.size()), 0,
        rowMapOf[members[0]]->getComm()));
}

// Assemble the (rowMembers × colMembers) super-block of A_blocked into a
// single monolithic Tpetra::CrsMatrix on the given group row/domain maps —
// the matrix the application would have built directly had it been called
// with this grouping in the first place. Off-diagonal coupling between
// members is included. Returns null if every member sub-block is null
// (a structurally zero super-block).
//
// Multi-rank: each member block keeps its own row distribution; only
// locally-owned rows are inserted. Column indices go through the member's
// column map to recover the member-global column (the local index alone is
// not the within-block column index once the block is distributed).
inline Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, Node>> assembleGroupBlock(
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>>                A_blocked,
    const std::vector<int>&                                           rowMembers,
    const std::vector<int>&                                           colMembers,
    const std::vector<int>&                                           rowOffsets,
    const std::vector<int>&                                           colOffsets,
    const std::vector<Teuchos::RCP<const Tpetra::Map<LO, GO, Node>>>& rowMapOf,
    Teuchos::RCP<const Tpetra::Map<LO, GO, Node>>                     groupRowMap,
    Teuchos::RCP<const Tpetra::Map<LO, GO, Node>>                     groupDomainMap)
{
    using CrsMatrix = Tpetra::CrsMatrix<SC, LO, GO, Node>;

    const int mr = static_cast<int>(rowMembers.size());
    const int mc = static_cast<int>(colMembers.size());

    std::vector<std::vector<Teuchos::RCP<const CrsMatrix>>> sub(
        mr, std::vector<Teuchos::RCP<const CrsMatrix>>(mc));
    bool anyNonNull = false;
    for (int si = 0; si < mr; ++si)
        for (int sj = 0; sj < mc; ++sj) {
            sub[si][sj] = getBlockCrs(A_blocked, rowMembers[si], colMembers[sj]);
            anyNonNull  = anyNonNull || !sub[si][sj].is_null();
        }
    if (!anyNonNull) return Teuchos::null;

    // Pass 1: count nonzeros per local group row (groupRowMap's local
    // ordering = members' local rows concatenated in member order).
    Teuchos::Array<size_t> nnzPerRow(groupRowMap->getLocalNumElements(), 0);
    {
        size_t lrow = 0;
        for (int si = 0; si < mr; ++si) {
            const size_t nloc = rowMapOf[rowMembers[si]]->getLocalNumElements();
            for (size_t lr = 0; lr < nloc; ++lr, ++lrow)
                for (int sj = 0; sj < mc; ++sj)
                    if (!sub[si][sj].is_null())
                        nnzPerRow[lrow] += sub[si][sj]->getNumEntriesInLocalRow(static_cast<LO>(lr));
        }
    }

    auto M = Teuchos::rcp(new CrsMatrix(groupRowMap, nnzPerRow()));

    // Pass 2: copy each member block's locally-owned entries in with row/col
    // offsets.
    for (int si = 0; si < mr; ++si) {
        const auto&  rmap = rowMapOf[rowMembers[si]];
        const size_t nloc = rmap->getLocalNumElements();
        for (int sj = 0; sj < mc; ++sj) {
            const auto& crs = sub[si][sj];
            if (crs.is_null()) continue;
            auto colMap = crs->getColMap();
            for (size_t lr = 0; lr < nloc; ++lr) {
                typename CrsMatrix::local_inds_host_view_type lcols;
                typename CrsMatrix::values_host_view_type     lvals;
                crs->getLocalRowView(static_cast<LO>(lr), lcols, lvals);
                const int nnz = static_cast<int>(lcols.extent(0));
                if (nnz == 0) continue;
                Teuchos::Array<GO> gcols(nnz);
                Teuchos::Array<SC> vals(nnz);
                for (int k = 0; k < nnz; ++k) {
                    gcols[k] = static_cast<GO>(colOffsets[sj]) + colMap->getGlobalElement(lcols(k));
                    vals[k]  = lvals(k);
                }
                const GO grow = static_cast<GO>(rowOffsets[si]) +
                                rmap->getGlobalElement(static_cast<LO>(lr));
                M->insertGlobalValues(grow, gcols(), vals());
            }
        }
    }
    M->fillComplete(groupDomainMap, groupRowMap);
    return M;
}

// Copy a 1-column product multivector from the original nb-block layout into
// the grouped nb_new-block layout (group block = members' local parts
// concatenated in member order — the layout induced by buildGroupMap).
// Purely local data movement, no communication.
inline void copyOriginalToGrouped(
    Teuchos::RCP<const MV>               orig,
    Teuchos::RCP<MV>                     grouped,
    const std::vector<std::vector<int>>& groups)
{
    using TMV = Thyra::TpetraMultiVector<SC, LO, GO, Node>;
    auto o  = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC>>(orig, true);
    auto gp = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<SC>>(grouped, true);
    for (size_t g = 0; g < groups.size(); ++g) {
        auto gT = Teuchos::rcp_dynamic_cast<TMV>(
                      gp->getNonconstMultiVectorBlock(static_cast<int>(g)), true)
                      ->getTpetraMultiVector();
        auto gview = gT->getVectorNonConst(0)->getDataNonConst();
        size_t off = 0;
        for (int k : groups[g]) {
            auto mT = Teuchos::rcp_dynamic_cast<const TMV>(
                          o->getMultiVectorBlock(k), true)
                          ->getConstTpetraMultiVector();
            auto mview = mT->getVector(0)->getData();
            const size_t nloc = mT->getLocalLength();
            for (size_t r = 0; r < nloc; ++r) gview[off + r] = mview[r];
            off += nloc;
        }
    }
}

// Inverse of copyOriginalToGrouped: unpack a grouped product multivector back
// into the original nb-block layout.
inline void copyGroupedToOriginal(
    Teuchos::RCP<const MV>               grouped,
    Teuchos::RCP<MV>                     orig,
    const std::vector<std::vector<int>>& groups)
{
    using TMV = Thyra::TpetraMultiVector<SC, LO, GO, Node>;
    auto gp = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC>>(grouped, true);
    auto o  = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<SC>>(orig, true);
    for (size_t g = 0; g < groups.size(); ++g) {
        auto gT = Teuchos::rcp_dynamic_cast<const TMV>(
                      gp->getMultiVectorBlock(static_cast<int>(g)), true)
                      ->getConstTpetraMultiVector();
        auto gview = gT->getVector(0)->getData();
        size_t off = 0;
        for (int k : groups[g]) {
            auto mT = Teuchos::rcp_dynamic_cast<TMV>(
                          o->getNonconstMultiVectorBlock(k), true)
                          ->getTpetraMultiVector();
            auto mview = mT->getVectorNonConst(0)->getDataNonConst();
            const size_t nloc = mT->getLocalLength();
            for (size_t r = 0; r < nloc; ++r) mview[r] = gview[off + r];
            off += nloc;
        }
    }
}

// Result of buildReconfiguredPrec: a freshly assembled flat nb_new x nb_new
// blocked system in which every merged group is a single monolithic Tpetra
// block — exactly the operator the application would have constructed had it
// been called with this grouping originally — together with the new
// preconditioner and the group membership (needed to repack the RHS into,
// and the solution out of, the grouped layout).
//
// Building the system this way (instead of reordering A_blocked into nested
// sub-blocks with per-apply gather/scatter wrapper inverses) makes the second
// solve structurally identical to a first solve: plain Tpetra blocks, plain
// block inverse ops, flat product vectors. Its factor and iterate times are
// therefore directly comparable to solve 1's — they measure the realistic
// cost of having used the reconfigured grouping from the start.
struct ReconfiguredSystem {
    Teko::LinearOp                 precOp;
    Teko::BlockedLinearOp          flatOp;
    std::vector<std::vector<int>>  groups;
    double                         factor_wall_time_sec;
};

// Build a Stratimikos-backed Ifpack2 inverse factory. This is configuration
// setup only (no factorization happens here, so it is outside the factor-time
// accounting); build it once and reuse it across a sweep rather than
// rebuilding per ordering.
inline Teuchos::RCP<Teko::InverseFactory> makeIfpack2InverseFactory()
{
    Stratimikos::DefaultLinearSolverBuilder builder;
    auto stratParams = Teuchos::rcp(new Teuchos::ParameterList);
    stratParams->set("Linear Solver Type",  "Belos");
    stratParams->set("Preconditioner Type", "Ifpack2");
    builder.setParameterList(stratParams);
    auto invLib = Teko::InverseLibrary::buildFromStratimikos(builder);
    return invLib->getInverseFactory("Ifpack2");
}

// Build the reconfigured system and preconditioner from an ordering vector.
//
// ordering[k] = g means old block k maps to new group g.
// Groups with a single member reuse the original blocks directly. Groups with
// multiple members are assembled into monolithic super-blocks (off-diagonal
// coupling included), so a merged group's diagonal is factored jointly and
// the whole system looks exactly as if the application had been called with
// the merged grouping in the first place.
//
// invFact is the (shared) Ifpack2 inverse factory; only the per-group
// buildInverse calls below are timed into factor_wall_time_sec, so passing a
// prebuilt factory does not affect the reported timing.
inline ReconfiguredSystem buildReconfiguredPrec(
    const std::vector<int>&                            ordering,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>> A_blocked,
    int                                                nb_old,
    const std::string&                                 method,
    Teuchos::RCP<Teko::InverseFactory>                 invFact)
{
    // Validate the externally-supplied ordering before indexing with it:
    // a malformed reconfiguration.json (wrong length, negative entries)
    // would otherwise cause out-of-bounds vector accesses below.
    TEUCHOS_TEST_FOR_EXCEPTION(ordering.size() != static_cast<size_t>(nb_old),
        std::runtime_error,
        "Teko::KrylovSurrogate::buildReconfiguredPrec: reconfiguration.json "
        "ordering has " + std::to_string(ordering.size()) + " entries but "
        "the operator has " + std::to_string(nb_old) + " blocks.");
    for (int k = 0; k < nb_old; ++k)
        TEUCHOS_TEST_FOR_EXCEPTION(ordering[k] < 0, std::runtime_error,
            "Teko::KrylovSurrogate::buildReconfiguredPrec: ordering[" +
            std::to_string(k) + "] = " + std::to_string(ordering[k]) +
            " is negative.");

    // Determine group membership
    const int nb_new = *std::max_element(ordering.begin(), ordering.end()) + 1;
    std::vector<std::vector<int>> groups(nb_new);
    for (int k = 0; k < nb_old; ++k)
        groups[ordering[k]].push_back(k);

    // Block sizes (needed to assemble merged super-blocks with row/col offsets).
    std::vector<int> block_sizes(nb_old);
    for (int k = 0; k < nb_old; ++k)
        block_sizes[k] =
            static_cast<int>(A_blocked->productRange()->getBlock(k)->dim());

    // A reconfiguration.json whose ordering values skip an integer (e.g.
    // [0, 2, 0] — group 1 has no members) would leave groups[1] empty,
    // leading to a 0-sized sub-blocked operator below and a confusing
    // failure deep inside Thyra/Teko. Reject it here with a clear message.
    for (int g = 0; g < nb_new; ++g)
        TEUCHOS_TEST_FOR_EXCEPTION(groups[g].empty(), std::runtime_error,
            "Teko::KrylovSurrogate::buildReconfiguredPrec: reconfiguration.json "
            "ordering has no entries mapping to group " + std::to_string(g) +
            " (group ids must be contiguous starting from 0; max ordering "
            "value implies " + std::to_string(nb_new) + " groups).");

    // invFact (Ifpack2 inverse factory) is supplied by the caller and reused.

    // ── Per-block row maps and per-group maps/spaces ───────────────────────
    std::vector<Teuchos::RCP<const Tpetra::Map<LO, GO, Node>>> rowMapOf(nb_old);
    for (int k = 0; k < nb_old; ++k)
        rowMapOf[k] = blockRowMap(A_blocked, k, nb_old);

    std::vector<std::vector<int>> offsets(nb_new);  // within-group member offsets
    std::vector<Teuchos::RCP<const Tpetra::Map<LO, GO, Node>>>  groupMap(nb_new);
    std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>> groupSpace(nb_new);
    for (int g = 0; g < nb_new; ++g) {
        const auto& mem = groups[g];
        offsets[g].assign(mem.size(), 0);
        for (size_t s = 1; s < mem.size(); ++s)
            offsets[g][s] = offsets[g][s - 1] + block_sizes[mem[s - 1]];
        groupMap[g] = buildGroupMap(mem, rowMapOf, offsets[g]);
        // Singleton groups keep the original Thyra space (so their blocks can
        // be reused as-is); merged groups get a fresh space over the group map.
        groupSpace[g] = (mem.size() == 1)
            ? A_blocked->productRange()->getBlock(mem[0])
            : Teuchos::RCP<const Thyra::VectorSpaceBase<SC>>(
                  Thyra::tpetraVectorSpace<SC, LO, GO, Node>(groupMap[g]));
    }

    // ── Assemble the flat nb_new × nb_new operator ─────────────────────────
    // Singleton×singleton blocks are reused from A_blocked without copying;
    // any block touching a merged group is assembled monolithically.
    auto flat = Thyra::defaultBlockedLinearOp<SC>();
    flat->beginBlockFill(nb_new, nb_new);
    for (int g = 0; g < nb_new; ++g) {
        for (int h = 0; h < nb_new; ++h) {
            if (groups[g].size() == 1 && groups[h].size() == 1) {
                auto blk = A_blocked->getBlock(groups[g][0], groups[h][0]);
                if (!blk.is_null()) {
                    flat->setBlock(g, h, blk);
                    continue;
                }
            } else {
                auto M = assembleGroupBlock(A_blocked, groups[g], groups[h],
                                            offsets[g], offsets[h], rowMapOf,
                                            groupMap[g], groupMap[h]);
                if (!M.is_null()) {
                    flat->setBlock(g, h,
                        Thyra::tpetraLinearOp<SC, LO, GO, Node>(
                            groupSpace[g], groupSpace[h], M));
                    continue;
                }
            }
            flat->setBlock(g, h, Thyra::zero<SC>(groupSpace[g], groupSpace[h]));
        }
    }
    flat->endBlockFill();

    Teko::BlockedLinearOp flatTeko =
        Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<SC>>(flat);

    // ── Factor each group diagonal (this is solve 2's setup cost) ──────────
    double factor_wall_time_sec = 0.0;
    std::vector<Teko::LinearOp> newInvDiag(nb_new);
    for (int g = 0; g < nb_new; ++g) {
        auto diag = flatTeko->getBlock(g, g);
        try {
            const auto factorStart = std::chrono::steady_clock::now();
            newInvDiag[g] = Teko::buildInverse(*invFact, diag);
            factor_wall_time_sec += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - factorStart).count();
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "Teko::KrylovSurrogate::buildReconfiguredPrec: factorization "
                "of group " + std::to_string(g) + " failed (singular?): " +
                e.what());
        }
    }

    // ── Assemble the preconditioner ────────────────────────────────────────
    Teko::LinearOp precOp;
    if (method == "jacobi")
        precOp = Teko::createBlockDiagonalInverseOp(flatTeko, newInvDiag);
    else
        precOp = Teko::createBlockLowerTriInverseOp(flatTeko, newInvDiag);

    return ReconfiguredSystem{precOp, flatTeko, groups, factor_wall_time_sec};
}

// ═══════════════════════════════════════════════════════════════════════════
// Section 5: orchestration — the hook registered with BelosAdaptiveHook
// ═══════════════════════════════════════════════════════════════════════════

// Default location for s<N>_request.json / s<N>_reconfig.json / s<N>_conv.json
// when TEKO_RECONFIG_REQUESTS_DIR is unset. All three file types live together
// in this single directory. Hardcoded for now, relative to the Trilinos
// checkout this package lives in.
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

// Convergence diagnostics for one FGMRES solve, written to
// s<N>_conv.json by writeConvergenceJson().
struct SolveStats {
    int    iterations;
    double initial_residual;
    double final_residual;
    double factor_wall_time_sec;
    double iterate_wall_time_sec;
    double total_wall_time_sec;
};

// Write s<request_id>_conv.json (atomically, via a .tmp file + rename).
// s2 == nullptr means no second solve was attempted (Phase 2 timed out waiting
// for s<request_id>_reconfig.json), and is recorded as "solve2": null.
inline void writeConvergenceJson(
    const std::string& convergence_dir,
    int                 request_id,
    const SolveStats&   s1,
    const SolveStats*   s2)
{
    fs::create_directories(convergence_dir);

    const std::string path     = convergence_dir + "/s" + std::to_string(request_id) + "_conv.json";
    const std::string tmp_path = path + ".tmp";
    std::ofstream ofs(tmp_path);
    TEUCHOS_TEST_FOR_EXCEPTION(!ofs.is_open(), std::runtime_error,
        "Teko::KrylovSurrogate::writeConvergenceJson: cannot open " + tmp_path);

    auto writeStats = [&](const char* name, const SolveStats& s, bool trailing_comma) {
        ofs << "  \"" << name << "\": {\n";
        ofs << "    \"iterations\": " << s.iterations << ",\n";
        ofs << "    \"initial_residual\": " << std::setprecision(17) << s.initial_residual << ",\n";
        ofs << "    \"final_residual\": " << std::setprecision(17) << s.final_residual << ",\n";
        ofs << "    \"factor_wall_time_sec\": " << std::setprecision(17) << s.factor_wall_time_sec << ",\n";
        ofs << "    \"iterate_wall_time_sec\": " << std::setprecision(17) << s.iterate_wall_time_sec << ",\n";
        ofs << "    \"total_wall_time_sec\": " << std::setprecision(17) << s.total_wall_time_sec << ",\n";
        ofs << "    \"wall_time_sec\": " << std::setprecision(17) << s.total_wall_time_sec << "\n";
        ofs << "  }" << (trailing_comma ? ",\n" : "\n");
    };

    ofs << "{\n";
    ofs << "  \"request_id\": " << request_id << ",\n";
    writeStats("solve1", s1, /*trailing_comma=*/true);
    if (s2 != nullptr)
        writeStats("solve2", *s2, /*trailing_comma=*/false);
    else
        ofs << "  \"solve2\": null\n";
    ofs << "}\n";
    ofs.close();

    fs::rename(tmp_path, path);
    std::cout << "[TekoAdaptive] wrote " << path << "\n";
}

// Per-ordering solve outcome (one row of s<N>_solved.json).
struct OrderingResult {
    std::vector<int> ordering;
    int    iterations           = 0;
    bool   converged            = false;
    double initial_residual     = 0.0;
    double final_residual       = 0.0;
    double factor_wall_time_sec = 0.0;
    double iterate_wall_time_sec = 0.0;
    double total_wall_time_sec  = 0.0;
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

// Write s<request_id>_solved.json (atomically): the full test-ordering sweep
// with per-ordering results, the selection_mode, and which ordering was
// selected/returned.
inline void writeSolvedJson(
    const std::string&                 dir,
    int                                request_id,
    const std::string&                 selection_mode,
    const std::vector<OrderingResult>& results,
    int                                selected_index,
    const std::vector<int>&            selected_ordering)
{
    fs::create_directories(dir);
    const std::string path = dir + "/s" + std::to_string(request_id) + "_solved.json";
    const std::string tmp_path = path + ".tmp";
    std::ofstream ofs(tmp_path);
    TEUCHOS_TEST_FOR_EXCEPTION(!ofs.is_open(), std::runtime_error,
        "Teko::KrylovSurrogate::writeSolvedJson: cannot open " + tmp_path);

    auto writeOrdering = [&](const std::vector<int>& o) {
        ofs << "[";
        for (std::size_t k = 0; k < o.size(); ++k) { if (k) ofs << ", "; ofs << o[k]; }
        ofs << "]";
    };

    ofs << "{\n";
    ofs << "  \"selection_mode\": \"" << selection_mode << "\",\n";
    ofs << "  \"request_id\": " << request_id << ",\n";
    ofs << "  \"selected_index\": " << selected_index << ",\n";
    ofs << "  \"selected_ordering\": "; writeOrdering(selected_ordering); ofs << ",\n";
    ofs << "  \"results\": [\n";
    for (std::size_t i = 0; i < results.size(); ++i) {
        const OrderingResult& r = results[i];
        ofs << "    {\n";
        ofs << "      \"ordering\": "; writeOrdering(r.ordering); ofs << ",\n";
        ofs << "      \"iterations\": " << r.iterations << ",\n";
        ofs << "      \"converged\": " << (r.converged ? "true" : "false") << ",\n";
        ofs << "      \"initial_residual\": " << std::setprecision(17) << r.initial_residual << ",\n";
        ofs << "      \"final_residual\": " << std::setprecision(17) << r.final_residual << ",\n";
        ofs << "      \"factor_wall_time_sec\": " << std::setprecision(17) << r.factor_wall_time_sec << ",\n";
        ofs << "      \"iterate_wall_time_sec\": " << std::setprecision(17) << r.iterate_wall_time_sec << ",\n";
        ofs << "      \"total_wall_time_sec\": " << std::setprecision(17) << r.total_wall_time_sec << "\n";
        ofs << "    }" << (i + 1 < results.size() ? "," : "") << "\n";
    }
    ofs << "  ]\n";
    ofs << "}\n";
    ofs.close();

    fs::rename(tmp_path, path);
    std::cout << "[TekoAdaptive] wrote " << path << "\n";
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
