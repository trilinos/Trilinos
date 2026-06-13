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
//   3. Reads the ordering vector, rebuilds the Teko preconditioner.
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
#include "Teko_BlockedReordering.hpp"
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

    const std::string prefix = "s";
    const std::string suffix = "_request.json";
    for (const auto& entry : fs::directory_iterator(requests_dir)) {
        if (!entry.is_regular_file()) continue;
        const std::string name = entry.path().filename().string();
        if (name.size() <= prefix.size() + suffix.size()) continue;
        if (name.compare(0, prefix.size(), prefix) != 0) continue;
        if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0) continue;

        const std::string digits = name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
        if (digits.empty() || !std::all_of(digits.begin(), digits.end(),
                                            [](char c){ return std::isdigit(static_cast<unsigned char>(c)); })) continue;

        const int n = std::stoi(digits);
        if (n + 1 > next) next = n + 1;
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

// Poll for s<request_id>_reconfig.json and read its ordering vector. The file
// is left in place after being read (NOT deleted) — its presence is how
// wait_for_request.py recognizes an already-answered request. Returns empty
// vector on timeout.
inline std::vector<int> waitForOrdering(
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
            // Intentionally not removed: leaving s<N>_reconfig.json in place
            // marks this request as already answered for the watcher.

            // Parse  {"ordering": [0, 2, 1, 1]}
            std::vector<int> ordering;
            const auto lb = content.find('[');
            const auto rb = content.find(']');
            if (lb == std::string::npos || rb == std::string::npos)
                throw std::runtime_error(
                    "[TekoAdaptive] " + path + ": no [...] found");

            std::istringstream ss(content.substr(lb + 1, rb - lb - 1));
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                tok.erase(std::remove_if(tok.begin(), tok.end(),
                                         [](char c){ return std::isspace(c); }),
                          tok.end());
                if (!tok.empty())
                    ordering.push_back(std::stoi(tok));
            }
            std::cout << "[TekoAdaptive] ordering: [";
            for (int k = 0; k < (int)ordering.size(); ++k) {
                if (k) std::cout << ", ";
                std::cout << ordering[k];
            }
            std::cout << "]\n";
            return ordering;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }

    std::cerr << "[TekoAdaptive] timed out waiting for " << path << "\n";
    return {};
}

// ═══════════════════════════════════════════════════════════════════════════
// Section 4: Preconditioner rebuild from ordering
// ═══════════════════════════════════════════════════════════════════════════

// Assemble the coupled super-block of a merged group as a single monolithic
// Tpetra::CrsMatrix, returned wrapped as a Thyra TpetraLinearOp.
//
// For members = [k0, k1, ...] the result is an N x N matrix (N = sum of the
// members' block sizes) whose (si, sj) sub-block is A_blocked->getBlock(
// members[si], members[sj]) — INCLUDING the off-diagonal coupling between
// members — placed at the corresponding row/column offset. This is the matrix
// that Ifpack2 ILUT factors jointly so a merged group is preconditioned as if
// its members were a single block (rather than block-diagonally).
//
// Multi-rank: each member block keeps its own row distribution. The merged
// map's local rows on each rank are the concatenation, in member order, of
// that rank's local rows of every member (with the member's global row index
// shifted by the member's offset). The flat vector layout induced by this map
// is therefore exactly "member-local parts back to back", which is what
// MergedGroupInverseOp::implicitApply's gather/scatter assumes.
inline Teko::LinearOp assembleMergedBlock(
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>> A_blocked,
    const std::vector<int>&                            members,
    const std::vector<int>&                            block_sizes)
{
    using CrsMatrix = Tpetra::CrsMatrix<SC, LO, GO, Node>;
    using TpMap     = Tpetra::Map<LO, GO, Node>;

    const int m = static_cast<int>(members.size());

    // Offsets of each member within the merged block.
    std::vector<int> off(m + 1, 0);
    for (int s = 0; s < m; ++s) off[s + 1] = off[s] + block_sizes[members[s]];

    // Extract every member sub-block (si, sj) as a Tpetra CrsMatrix, mirroring
    // the dynamic-cast chain in computeCHat(). Null blocks (structural zeros)
    // are simply skipped during assembly.
    std::vector<std::vector<Teuchos::RCP<const CrsMatrix>>> sub(
        m, std::vector<Teuchos::RCP<const CrsMatrix>>(m));
    Teuchos::RCP<const Teuchos::Comm<int>> comm;
    for (int si = 0; si < m; ++si) {
        for (int sj = 0; sj < m; ++sj) {
            auto blk = A_blocked->getBlock(members[si], members[sj]);
            if (blk.is_null()) continue;
            auto lo = Teuchos::rcp_dynamic_cast<
                const Thyra::TpetraLinearOp<SC, LO, GO, Node>>(blk);
            if (lo.is_null()) continue;
            auto crs = Teuchos::rcp_dynamic_cast<const CrsMatrix>(
                lo->getConstTpetraOperator());
            if (crs.is_null()) continue;
            sub[si][sj] = crs;
            if (comm.is_null()) comm = crs->getRowMap()->getComm();
        }
    }
    TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), std::runtime_error,
        "Teko::KrylovSurrogate::assembleMergedBlock: no Tpetra blocks found in "
        "the merged group.");

    // Row map of each member = row map of the first non-null block in its row.
    std::vector<Teuchos::RCP<const TpMap>> memberRowMap(m);
    for (int si = 0; si < m; ++si) {
        for (int sj = 0; sj < m; ++sj)
            if (!sub[si][sj].is_null()) { memberRowMap[si] = sub[si][sj]->getRowMap(); break; }
        TEUCHOS_TEST_FOR_EXCEPTION(memberRowMap[si].is_null(), std::runtime_error,
            "Teko::KrylovSurrogate::assembleMergedBlock: member " +
            std::to_string(members[si]) + " has no Tpetra blocks in its row.");
    }

    // Merged map: this rank's GIDs are each member's local rows (as global
    // indices) shifted by the member's offset, concatenated in member order.
    std::vector<GO> myGIDs;
    for (int si = 0; si < m; ++si) {
        const size_t nloc = memberRowMap[si]->getLocalNumElements();
        for (size_t lr = 0; lr < nloc; ++lr)
            myGIDs.push_back(static_cast<GO>(off[si]) +
                             memberRowMap[si]->getGlobalElement(static_cast<LO>(lr)));
    }
    auto map = Teuchos::rcp(new TpMap(
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
        Teuchos::ArrayView<const GO>(myGIDs.data(), myGIDs.size()), 0, comm));

    // Pass 1: count nonzeros per local merged row (same ordering as myGIDs).
    Teuchos::Array<size_t> nnzPerRow(myGIDs.size(), 0);
    {
        size_t lrow = 0;
        for (int si = 0; si < m; ++si) {
            const size_t nloc = memberRowMap[si]->getLocalNumElements();
            for (size_t lr = 0; lr < nloc; ++lr, ++lrow)
                for (int sj = 0; sj < m; ++sj)
                    if (!sub[si][sj].is_null())
                        nnzPerRow[lrow] += sub[si][sj]->getNumEntriesInLocalRow(static_cast<LO>(lr));
        }
    }

    auto M = Teuchos::rcp(new CrsMatrix(map, nnzPerRow()));

    // Pass 2: copy each member block's locally-owned entries in with row/col
    // offsets. Column indices go through the member's column map to recover
    // the member-global column (the local index alone is not the within-block
    // column index once the block is distributed).
    for (int si = 0; si < m; ++si) {
        const size_t nloc = memberRowMap[si]->getLocalNumElements();
        for (int sj = 0; sj < m; ++sj) {
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
                    gcols[k] = static_cast<GO>(off[sj]) + colMap->getGlobalElement(lcols(k));
                    vals[k]  = lvals(k);
                }
                const GO grow = static_cast<GO>(off[si]) +
                                memberRowMap[si]->getGlobalElement(static_cast<LO>(lr));
                M->insertGlobalValues(grow, gcols(), vals());
            }
        }
    }
    M->fillComplete(map, map);

    return Thyra::tpetraLinearOp<SC, LO, GO, Node>(
        Thyra::tpetraVectorSpace<SC, LO, GO, Node>(M->getRangeMap()),
        Thyra::tpetraVectorSpace<SC, LO, GO, Node>(M->getDomainMap()),
        M);
}

// Inverse operator for a merged group that applies a single joint factorization
// of the coupled super-block. The merged block is assembled into one flat
// N x N matrix and factored once (Ifpack2 ILUT) into the flat inverse `Minv`;
// but the surrounding BlockLowerTriInverseOp / BlockDiagonalInverseOp feed this
// op the group's vector as a *nested* product multivector (one sub-block per
// member). implicitApply() therefore gathers the member sub-blocks into a flat
// vector, applies Minv, and scatters the result back — bridging the
// product-of-members layout and the flat monolithic factor.
//
// range()/domain() match the merged group's blocked sub-operator (subTeko),
// exactly as the previous createBlockDiagonalInverseOp(subTeko, ...) did, so no
// other part of buildReconfiguredPrec / adaptiveLoop needs to change.
//
// Multi-rank: assembleMergedBlock builds the merged map so that this rank's
// local rows are the members' local rows concatenated in member order. The
// gather/scatter below therefore copies each member's *local* part (local
// length, not global block size) into consecutive local slots of the flat
// vector — purely local data movement, no communication.
class MergedGroupInverseOp : public Teko::BlockImplicitLinearOp {
public:
    MergedGroupInverseOp(const Teko::BlockedLinearOp& subTeko,
                         const Teko::LinearOp&        Minv,
                         int                          numMembers)
        : Minv_(Minv),
          numMembers_(numMembers),
          productRange_(subTeko->productDomain()),    // inverse: flip range/domain
          productDomain_(subTeko->productRange()) {}

    Teko::VectorSpace range()  const override { return productRange_;  }
    Teko::VectorSpace domain() const override { return productDomain_; }

    using Teko::BlockImplicitLinearOp::implicitApply;
    void implicitApply(const Teko::BlockedMultiVector& src,
                       Teko::BlockedMultiVector&       dst,
                       const double alpha = 1.0, const double beta = 0.0) const override
    {
        using TMV = Thyra::TpetraMultiVector<SC, LO, GO, Node>;
        const int m     = numMembers_;
        const int ncols = src->domain()->dim();

        auto flatSrc = Thyra::createMembers(Minv_->domain(), ncols);
        auto flatDst = Thyra::createMembers(Minv_->range(),  ncols);

        // Gather: member sub-blocks' local parts → flat source.
        {
            auto fT = Teuchos::rcp_dynamic_cast<TMV>(flatSrc, true)
                          ->getTpetraMultiVector();
            size_t o = 0;
            for (int b = 0; b < m; ++b) {
                auto blkT = Teuchos::rcp_dynamic_cast<const TMV>(
                                Teko::getBlock(b, src), true)
                                ->getConstTpetraMultiVector();
                const size_t nloc = blkT->getLocalLength();
                for (int c = 0; c < ncols; ++c) {
                    auto vin  = blkT->getVector(c)->getData();
                    auto vout = fT->getVectorNonConst(c)->getDataNonConst();
                    for (size_t r = 0; r < nloc; ++r)
                        vout[o + r] = vin[r];
                }
                o += nloc;
            }
        }

        // Apply the single joint inverse.
        Teko::applyOp(Minv_, flatSrc, flatDst);

        // Scatter: flat result → member sub-blocks of dst, honoring alpha/beta
        // (dst_block = alpha * (Minv*src)_block + beta * dst_block_old).
        {
            auto fT = Teuchos::rcp_dynamic_cast<const TMV>(flatDst, true)
                          ->getConstTpetraMultiVector();
            size_t o = 0;
            for (int b = 0; b < m; ++b) {
                auto blkT = Teuchos::rcp_dynamic_cast<TMV>(
                                Teko::getBlock(b, dst), true)
                                ->getTpetraMultiVector();
                const size_t nloc = blkT->getLocalLength();
                for (int c = 0; c < ncols; ++c) {
                    auto vin  = fT->getVector(c)->getData();
                    auto vout = blkT->getVectorNonConst(c)->getDataNonConst();
                    for (size_t r = 0; r < nloc; ++r)
                        vout[r] = alpha * vin[o + r] +
                                  (beta != 0.0 ? beta * vout[r] : 0.0);
                }
                o += nloc;
            }
        }
    }

private:
    Teko::LinearOp    Minv_;
    int               numMembers_;
    Teko::VectorSpace productRange_;
    Teko::VectorSpace productDomain_;
};

// Result of buildReconfiguredPrec: the new preconditioner together with the
// reordered operator and the BlockReorderManager that defines its (possibly
// nested) block structure. The second solve's operator/LHS/RHS must be
// reshaped to this same structure (via Teko::buildReorderedMultiVector)
// so that BlockLowerTriInverseOp::implicitApply's block-count assertions
// against reorderedOp hold.
struct ReconfiguredSystem {
    Teko::LinearOp                                 precOp;
    Teko::BlockedLinearOp                          reorderedOp;
    Teuchos::RCP<const Teko::BlockReorderManager>  mgr;
    double                                         factor_wall_time_sec;
};

// Build a new Teko block preconditioner from an ordering vector.
//
// ordering[k] = g means old block k maps to new group g.
// Groups with a single member: singleton (factorize that diagonal block).
// Groups with multiple members: merged — the coupled super-block (off-diagonals
// included) is assembled into one matrix and factored jointly, so the group is
// preconditioned as if its members were a single block.
inline ReconfiguredSystem buildReconfiguredPrec(
    const std::vector<int>&                            ordering,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>> A_blocked,
    int                                                nb_old,
    const std::string&                                 method)
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

    // Build Ifpack2 inverse factory
    Stratimikos::DefaultLinearSolverBuilder builder;
    auto stratParams = Teuchos::rcp(new Teuchos::ParameterList);
    stratParams->set("Linear Solver Type",  "Belos");
    stratParams->set("Preconditioner Type", "Ifpack2");
    builder.setParameterList(stratParams);
    auto invLib  = Teko::InverseLibrary::buildFromStratimikos(builder);
    auto invFact = invLib->getInverseFactory("Ifpack2");

    // ── Build a BlockReorderManager describing the new grouping, then use it
    // to reorder A_blocked into an nb_new x nb_new operator. Merged groups
    // become nested sub-blocked entries. The same manager is reused (via
    // Teko::buildReorderedMultiVector) to reshape the second solve's RHS/LHS
    // so their product-space block counts match tekoNew, which
    // BlockLowerTriInverseOp::implicitApply asserts on.
    auto mgr = Teuchos::rcp(new Teko::BlockReorderManager(nb_new));
    for (int g = 0; g < nb_new; ++g) {
        const auto& mem = groups[g];
        if (mem.size() == 1) {
            mgr->SetBlock(g, mem[0]);
        } else {
            auto sub = Teuchos::rcp(new Teko::BlockReorderManager(static_cast<int>(mem.size())));
            for (int si = 0; si < (int)mem.size(); ++si)
                sub->SetBlock(si, mem[si]);
            mgr->SetBlock(g, sub);
        }
    }

    Teko::LinearOp reordered = Teko::buildReorderedLinearOp(*mgr, A_blocked);
    Teko::BlockedLinearOp tekoNew =
        Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<SC>>(
            Teuchos::rcp_const_cast<Thyra::LinearOpBase<SC>>(reordered), true);

    // ── Build diagonal inverses for each new group ─────────────────────────
    double factor_wall_time_sec = 0.0;
    std::vector<Teko::LinearOp> newInvDiag(nb_new);
    for (int g = 0; g < nb_new; ++g) {
        const auto& mem = groups[g];
        if (mem.size() == 1) {
            // Singleton: factorize A_blocked(b,b) directly.
            const auto factorStart = std::chrono::steady_clock::now();
            newInvDiag[g] =
                Teko::buildInverse(*invFact,
                                   A_blocked->getBlock(mem[0], mem[0]));
            factor_wall_time_sec += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - factorStart).count();
        } else {
            // Merged group: factor the coupled super-block jointly.
            // subTeko is the sub-blocked operator (len(mem) × len(mem),
            // off-diagonals included). It supplies the nested product vector
            // spaces; the actual inverse comes from a single Ifpack2 ILUT
            // factorization of the monolithic assembly of those same blocks.
            auto subOp = Thyra::defaultBlockedLinearOp<SC>();
            subOp->beginBlockFill(static_cast<int>(mem.size()),
                                  static_cast<int>(mem.size()));
            for (int si = 0; si < (int)mem.size(); ++si) {
                for (int sj = 0; sj < (int)mem.size(); ++sj) {
                    auto blk = A_blocked->getBlock(mem[si], mem[sj]);
                    auto ri_sp = A_blocked->productRange()->getBlock(mem[si]);
                    auto di_sp = A_blocked->productDomain()->getBlock(mem[sj]);
                    if (!blk.is_null())
                        subOp->setBlock(si, sj, blk);
                    else
                        subOp->setBlock(si, sj,
                            Thyra::zero<SC>(ri_sp, di_sp));
                }
            }
            subOp->endBlockFill();

            Teko::BlockedLinearOp subTeko =
                Teuchos::rcp_dynamic_cast<Thyra::PhysicallyBlockedLinearOpBase<SC>>(subOp);

            // Assemble the coupled super-block as one matrix and factor it once.
            Teko::LinearOp mergedOp = assembleMergedBlock(A_blocked, mem, block_sizes);
            Teko::LinearOp mergedInv;
            try {
                const auto factorStart = std::chrono::steady_clock::now();
                mergedInv = Teko::buildInverse(*invFact, mergedOp);
                factor_wall_time_sec += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - factorStart).count();
            } catch (const std::exception& e) {
                throw std::runtime_error(
                    "Teko::KrylovSurrogate::buildReconfiguredPrec: joint "
                    "factorization of merged group " + std::to_string(g) +
                    " failed (singular?): " + e.what());
            }

            newInvDiag[g] = Teuchos::rcp(
                new MergedGroupInverseOp(subTeko, mergedInv,
                                         static_cast<int>(mem.size())));
        }
    }

    // ── Assemble the preconditioner ────────────────────────────────────────
    Teko::LinearOp precOp;
    if (method == "jacobi")
        precOp = Teko::createBlockDiagonalInverseOp(tekoNew, newInvDiag);
    else
        precOp = Teko::createBlockLowerTriInverseOp(tekoNew, newInvDiag);

    return ReconfiguredSystem{precOp, tekoNew, mgr, factor_wall_time_sec};
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
// (if it converges) would re-invoke this same hook recursively — and that
// recursive call operates on the *reordered* operator, whose merged-group
// blocks are nested ProductMultiVectors rather than plain TpetraMultiVectors,
// which getBlockTpetraMV() cannot handle. g_inAdaptiveLoop / ScopedGuard
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

// This is the function registered as Belos::AdaptiveHook::HookFn.
// It is called by BelosBlockGmresSolMgr::solve() after a flexible GMRES
// solve converges.
//
// On entry: problem->getLHS() contains the first solve's result.
// On exit:  problem->getLHS() is overwritten with the second solve's result
//           (if one was attempted and converged).
inline void adaptiveLoop(
    const Belos::AdaptiveHook::State&                         state,
    Teuchos::RCP<Belos::AdaptiveHook::Problem>                problem,
    Teuchos::RCP<const Thyra::BlockedLinearOpBase<SC>>        A_blocked,
    Teuchos::RCP<const Teuchos::ParameterList>                orig_params,
    const Belos::AdaptiveHook::SolveMetrics&                  solve1_metrics)
{
    // Bail out immediately if this call is the recursive re-entry from
    // Phase 4's solver2.solve() below (see ScopedAdaptiveLoopGuard).
    if (detail::g_inAdaptiveLoop) return;

    // requests_dir comes from TEKO_RECONFIG_REQUESTS_DIR if set, otherwise
    // kDefaultRequestsDir.
    const char* reconfig_dir_env = std::getenv("TEKO_RECONFIG_REQUESTS_DIR");
    const std::string requests_dir = reconfig_dir_env ? std::string(reconfig_dir_env) : kDefaultRequestsDir;

    if (A_blocked.is_null()) {
        std::cerr << "[TekoAdaptive] operator is not blocked; skipping.\n";
        return;
    }
    if (state.curDim == 0) {
        std::cerr << "[TekoAdaptive] curDim==0; nothing to compute.\n";
        return;
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

    SolveStats s1{state.curDim, b_norm, solve1_metrics.achieved_tol,
                  0.0, solve1_metrics.wall_time_sec, solve1_metrics.wall_time_sec};

    // ── Phase 1: compute C_hat, b_hat and write s<N>_request.json ──────────
    // computeCHat is collective (mvTransMv reductions + distributed applies);
    // its result is replicated identically on every rank, so only rank 0
    // writes the request file.
    CHatData chat = computeCHat(state, A_blocked, problem->getRHS(), block_sizes);
    int request_id = -1;
    if (rank == 0) request_id = writeRequestJson(chat, requests_dir);

    // ── Phase 2: wait for s<N>_reconfig.json ───────────────────────────────
    // The response arrives from an external process via the filesystem, so
    // rank 0 polls and the result is broadcast — the one genuinely
    // rank-0-sourced piece of data in the loop.
    std::vector<int> ordering;
    if (rank == 0) ordering = waitForOrdering(requests_dir, request_id);
    int ordering_size = static_cast<int>(ordering.size());
    Teuchos::broadcast(*comm, 0, 1, &ordering_size);
    ordering.resize(ordering_size);
    if (ordering_size > 0)
        Teuchos::broadcast(*comm, 0, ordering_size, ordering.data());

    if (ordering.empty()) {
        if (rank == 0) {
            std::cerr << "[TekoAdaptive] no ordering received; returning first "
                         "solve result.\n";
            writeConvergenceJson(convergence_dir, request_id, s1, nullptr);
        }
        return;
    }

    // ── Phase 3: rebuild Teko preconditioner ──────────────────────────────
    // Determine method from the existing preconditioner type.
    // Heuristic: if it dynamic-casts to BlockLowerTriInverseOp → "gs",
    //            otherwise "jacobi".
    std::string method = "jacobi";
    {
        auto lti = Teuchos::rcp_dynamic_cast<
            const Teko::BlockLowerTriInverseOp>(problem->getRightPrec());
        if (!lti.is_null()) method = "gs";
    }

    ReconfiguredSystem recon = buildReconfiguredPrec(ordering, A_blocked, nb, method);

    // ── Phase 4: second FGMRES solve ──────────────────────────────────────
    // Prevent the hook from firing recursively during the second solve;
    // reset on scope exit even if solver2.solve() throws.
    ScopedAdaptiveLoopGuard recursion_guard;

    auto newProblem = Teuchos::rcp(new Problem());
    // setOperator/setRightPrec take RCP<const OP>, so a plain dynamic_cast
    // (which can add const) suffices — no const_cast needed.
    newProblem->setOperator(Teuchos::rcp_dynamic_cast<const OP>(recon.reorderedOp));
    newProblem->setRightPrec(Teuchos::rcp_dynamic_cast<const OP>(recon.precOp));

    // Reshape the flat LHS/RHS to match the block structure of
    // recon.reorderedOp (recon.mgr's groups, with merged groups nested as
    // sub-product vectors) — required so blockCount(src) ==
    // blockRowCount(reorderedOp) inside BlockLowerTriInverseOp::implicitApply.
    // x_new_reordered shares storage with x_new (buildReorderedMultiVector
    // wraps the existing leaf blocks without copying), so writes made by the
    // second solve land directly in x_new's buffers.
    auto x_new = Thyra::createMembers(A_blocked->domain(), 1);
    Thyra::assign(x_new.ptr(), SC(0));
    auto x_flat = Teuchos::rcp_dynamic_cast<Thyra::ProductMultiVectorBase<SC>>(x_new, true);
    auto x_new_reordered = Teko::buildReorderedMultiVector(*recon.mgr, x_flat);

    auto b_flat = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC>>(
        problem->getRHS(), true);
    auto b_reordered = Teko::buildReorderedMultiVector(*recon.mgr, b_flat);

    newProblem->setLHS(x_new_reordered);
    newProblem->setRHS(b_reordered);
    newProblem->setProblem();

    // Inherit the first solve's resolved parameters (Maximum Iterations, Num
    // Blocks, Convergence Tolerance, Flexible Gmres, Verbosity, etc.) so the
    // second solve uses the same settings the caller configured for the
    // first, rather than independent hardcoded defaults.
    auto solverParams = Teuchos::rcp(new Teuchos::ParameterList(*orig_params));

    Belos::BlockGmresSolMgr<SC, MV, OP> solver2(newProblem, solverParams);

    // Wall time for solve 2's iteration loop only, mirroring solve1_metrics
    // (which times only solve()'s own loop, excluding preconditioner setup).
    // buildReconfiguredPrec() above is solve 2's "setup" and is excluded here
    // too, so the two wall times are directly comparable.
    const auto solve2Start = std::chrono::steady_clock::now();
    Belos::ReturnType ret2 = solver2.solve();
    const double solve2_iterate_wall_time_sec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solve2Start).count();

    SolveStats s2{solver2.getNumIters(), b_norm, solver2.achievedTol(),
                  recon.factor_wall_time_sec, solve2_iterate_wall_time_sec,
                  recon.factor_wall_time_sec + solve2_iterate_wall_time_sec};
    if (rank == 0) writeConvergenceJson(convergence_dir, request_id, s1, &s2);

    if (ret2 != Belos::Converged) {
        if (rank == 0)
            std::cerr << "[TekoAdaptive] second solve did not converge ("
                      << solver2.getNumIters() << " iters). Using first result.\n";
    } else {
        // Overwrite the original LHS (x_thyra in teko_ext.cpp) with the
        // result of the second solve.
        Thyra::assign(problem->getLHS().ptr(), *x_new);
        if (rank == 0)
            std::cout << "[TekoAdaptive] second solve converged ("
                      << solver2.getNumIters() << " iters).\n";
    }
    // env_guard restores TEKO_RECONFIG_REQUESTS_DIR here (or during unwinding
    // if solver2.solve() threw).
}

} // namespace KrylovSurrogate
} // namespace Teko
