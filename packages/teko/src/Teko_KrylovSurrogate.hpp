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
//   1. Computes C_hat and writes request_<N>.json to requests_dir, where
//      N = nextRequestNumber(requests_dir) (0, 1, 2, ... — never reused).
//   2. Polls for reconfiguration_<N>.json in the same directory.
//   3. Reads the ordering vector, rebuilds the Teko preconditioner.
//   4. Runs a second FGMRES solve and overwrites the LHS with the result.
//   5. Writes convergence_data/conv_<N>.json (sibling of requests_dir) with
//      iteration count, initial/final residual, and wall time for both
//      solves (solve2 is null if step 2 timed out).
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
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

// Filesystem (C++17)
#include <filesystem>
namespace fs = std::filesystem;

// ── Teuchos ───────────────────────────────────────────────────────────────
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
// Section 1: dense helpers
// ═══════════════════════════════════════════════════════════════════════════

// Extract the first nCols columns of block j from a Thyra ProductMultiVector
// into a column-major SerialDenseMatrix (nRows_j × nCols).
inline SDM extractBlockDense(
    Teuchos::RCP<const MV> mv,
    int j,
    int nRows_j,
    int nCols)
{
    auto pmv  = Teuchos::rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<SC>>(mv);
    TEUCHOS_TEST_FOR_EXCEPTION(pmv.is_null(), std::runtime_error,
        "Teko::KrylovSurrogate::extractBlockDense: mv is not a ProductMultiVectorBase");

    auto blk  = pmv->getMultiVectorBlock(j);
    auto tmv  = Teuchos::rcp_dynamic_cast<const Thyra::TpetraMultiVector<SC, LO, GO, Node>>(blk);
    TEUCHOS_TEST_FOR_EXCEPTION(tmv.is_null(), std::runtime_error,
        "Teko::KrylovSurrogate::extractBlockDense: block is not a TpetraMultiVector");

    auto tptr = tmv->getConstTpetraMultiVector();

    // The dense SVD/C_hat computation below indexes this block's data with
    // plain global row indices (0..nRows_j-1). That is only correct if the
    // block's map is entirely local to this rank, i.e. a single-rank
    // (serial) run. On a multi-rank run getData() returns a *local* view
    // whose length need not match nRows_j (the global block dimension),
    // which would silently read/write out of bounds. Fail loudly instead.
    TEUCHOS_TEST_FOR_EXCEPTION(tptr->getMap()->getComm()->getSize() != 1,
        std::runtime_error,
        "Teko::KrylovSurrogate::extractBlockDense: block " + std::to_string(j) +
        " is distributed across " +
        std::to_string(tptr->getMap()->getComm()->getSize()) +
        " MPI ranks; the Krylov-surrogate adaptive loop currently only "
        "supports single-rank (serial) runs.");
    TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(nRows_j) != tptr->getLocalLength(),
        std::runtime_error,
        "Teko::KrylovSurrogate::extractBlockDense: block " + std::to_string(j) +
        " expected " + std::to_string(nRows_j) + " rows but local length is " +
        std::to_string(tptr->getLocalLength()));

    SDM D(nRows_j, nCols);
    for (int col = 0; col < nCols; ++col) {
        auto view = tptr->getVector(col)->getData();
        for (int row = 0; row < nRows_j; ++row)
            D(row, col) = view[row];
    }
    return D;
}

// ═══════════════════════════════════════════════════════════════════════════
// Section 2: CHat data structure and computation
// ═══════════════════════════════════════════════════════════════════════════

struct CHatData {
    int                nb;          // number of blocks
    int                krylov_dim;  // curDim used
    std::vector<int>   block_sizes;
    std::vector<int>   ranks;       // r_j per block

    // Q[j] is (block_sizes[j] × ranks[j]) — left singular vectors of V_j.
    // Internal use only (needed to form C_hat and b_hat below); NOT written
    // to request.json since its size is O(block_sizes[j]) and would make the
    // request scale with the global problem size n.
    std::vector<SDM>               Q;
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
    chat.Q.resize(nb);
    chat.b_hat.resize(nb);
    chat.blocks.assign(nb, std::vector<SDM>(nb));

    Teuchos::LAPACK<int, SC> lapack;

    // ── Per-block SVD of V_j ──────────────────────────────────────────────
    // V[:,0:curDim] are the curDim Krylov vectors.
    // V_j = rows of block j of V → (block_sizes[j] × curDim).
    // SVD: V_j = U Sigma W^T  (thin, JOBU='S', JOBVT='S').
    //   U:   block_sizes[j] × curDim  (left singular vectors)
    //   S:   curDim                    (singular values)
    //   VT:  curDim × curDim           (right singular vectors, rows)

    // Also store Z_j for later use.
    std::vector<SDM> Z_dense(nb);          // Z_dense[j] = n_j × curDim
    std::vector<SDM> WTrunc(nb);           // WTrunc[j]  = curDim × r_j (right sv)
    std::vector<std::vector<SC>> Sigma(nb); // Sigma[j]   = first r_j values

    for (int j = 0; j < nb; ++j) {
        const int nj = block_sizes[j];

        // Extract V_j and Z_j
        SDM Vj = extractBlockDense(state.V, j, nj, curDim);
        Z_dense[j] = extractBlockDense(state.Z, j, nj, curDim);

        // GESVD: overwrites Vj.  Need a mutable copy of the data.
        std::vector<SC> A_buf(nj * curDim);
        for (int c = 0; c < curDim; ++c)
            for (int r = 0; r < nj; ++r)
                A_buf[c * nj + r] = Vj(r, c);  // column-major

        const int k   = std::min(nj, curDim);
        std::vector<SC> S(k);
        std::vector<SC> U_buf(nj * k);           // nj × k, column-major
        std::vector<SC> VT_buf(k * curDim);      // k × curDim, column-major

        // Workspace query
        int lwork = -1, info = 0;
        SC work_query;
        lapack.GESVD('S', 'S',
            nj, curDim,
            A_buf.data(), nj,
            S.data(),
            U_buf.data(), nj,
            VT_buf.data(), k,
            &work_query, lwork,
            nullptr,  // RWORK ignored for double
            &info);
        lwork = static_cast<int>(work_query);
        std::vector<SC> work(lwork);
        lapack.GESVD('S', 'S',
            nj, curDim,
            A_buf.data(), nj,
            S.data(),
            U_buf.data(), nj,
            VT_buf.data(), k,
            work.data(), lwork,
            nullptr,
            &info);
        TEUCHOS_TEST_FOR_EXCEPTION(info != 0, std::runtime_error,
            "Teko::KrylovSurrogate::computeCHat: GESVD failed (info=" +
            std::to_string(info) + ") for block " + std::to_string(j));

        // Truncate
        const SC threshold = (S[0] > 0.0 ? svd_tol * S[0] : 0.0);
        int rj = 0;
        while (rj < k && S[rj] >= threshold) ++rj;
        if (rj == 0) rj = 1;  // keep at least one direction

        chat.ranks[j] = rj;
        Sigma[j].assign(S.begin(), S.begin() + rj);

        // Q[j] = U[:,0:rj]  (nj × rj), column-major
        chat.Q[j].shape(nj, rj);
        for (int c = 0; c < rj; ++c)
            for (int r = 0; r < nj; ++r)
                chat.Q[j](r, c) = U_buf[c * nj + r];

        // b_hat[j] = Q[j]^T * b_j   (rj × 1)
        SDM bj = extractBlockDense(b, j, nj, 1);
        chat.b_hat[j].shape(rj, 1);
        chat.b_hat[j].multiply(Teuchos::TRANS, Teuchos::NO_TRANS,
                                1.0, chat.Q[j], bj, 0.0);

        // WTrunc[j] = VT[0:rj,:]^T  stored as (curDim × rj).
        // GESVD returns VT as (k × curDim), column-major with ldvt=k, i.e.
        // VT(row, col) = VT_buf[col * k + row]. We want
        // WTrunc[j](kd, r) = VT(r, kd) = VT_buf[kd * k + r].
        WTrunc[j].shape(curDim, rj);
        for (int kd = 0; kd < curDim; ++kd)
            for (int r = 0; r < rj; ++r)
                WTrunc[j](kd, r) = VT_buf[kd * k + r];
    }

    // ── Per (i,j) block of C_hat ──────────────────────────────────────────
    // C_hat[i][j] = Q_i^T  (A_ij Z_j W_j Sigma_j^{-1})
    // where A_ij Z_j is n_i × curDim, then × W_j (curDim × rj) → n_i × rj,
    // scaled by Sigma_j^{-1}, then left-multiplied by Q_i^T (ri × ni).

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

            const int nj = block_sizes[j];
            const int ni = block_sizes[i];

            // Build Tpetra MV for Z_j (n_j × curDim)
            auto Zj_tmv = Teuchos::rcp(new TpetraMV(A_crs->getDomainMap(), curDim));
            // As in extractBlockDense: the indexing below assumes block j's
            // domain map is entirely local (single-rank) and matches nj.
            TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(nj) != Zj_tmv->getLocalLength(),
                std::runtime_error,
                "Teko::KrylovSurrogate::computeCHat: A_blocked->getBlock(" +
                std::to_string(i) + "," + std::to_string(j) + ")'s domain map "
                "has local length " + std::to_string(Zj_tmv->getLocalLength()) +
                " but block_sizes[" + std::to_string(j) + "]=" + std::to_string(nj));
            for (int col = 0; col < curDim; ++col) {
                auto view = Zj_tmv->getVectorNonConst(col)->getDataNonConst();
                for (int row = 0; row < nj; ++row)
                    view[row] = Z_dense[j](row, col);
            }

            // AijZj = A_ij * Z_j  (n_i × curDim)
            auto AijZj_tmv = Teuchos::rcp(new TpetraMV(A_crs->getRangeMap(), curDim));
            TEUCHOS_TEST_FOR_EXCEPTION(static_cast<size_t>(ni) != AijZj_tmv->getLocalLength(),
                std::runtime_error,
                "Teko::KrylovSurrogate::computeCHat: A_blocked->getBlock(" +
                std::to_string(i) + "," + std::to_string(j) + ")'s range map "
                "has local length " + std::to_string(AijZj_tmv->getLocalLength()) +
                " but block_sizes[" + std::to_string(i) + "]=" + std::to_string(ni));
            A_crs->apply(*Zj_tmv, *AijZj_tmv, Teuchos::NO_TRANS, 1.0, 0.0);

            SDM AijZj(ni, curDim);
            for (int col = 0; col < curDim; ++col) {
                auto view = AijZj_tmv->getVector(col)->getData();
                for (int row = 0; row < ni; ++row)
                    AijZj(row, col) = view[row];
            }

            // T = AijZj * WTrunc[j]  (n_i × rj)   via DGEMM
            // T(ni, rj) = AijZj(ni, curDim) * WTrunc[j](curDim, rj)
            SDM T(ni, rj);
            T.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS,
                       1.0, AijZj, WTrunc[j], 0.0);

            // Scale columns of T by 1/Sigma[j][c]
            for (int c = 0; c < rj; ++c) {
                const SC inv_s = (std::abs(Sigma[j][c]) > 1e-14)
                                 ? 1.0 / Sigma[j][c] : 0.0;
                for (int r = 0; r < ni; ++r)
                    T(r, c) *= inv_s;
            }

            // C_hat[i][j] = Q_i^T * T   (ri × rj)
            // Q_i is (ni × ri), so Q_i^T is (ri × ni)
            chat.blocks[i][j].multiply(Teuchos::TRANS, Teuchos::NO_TRANS,
                                       1.0, chat.Q[i], T, 0.0);
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

// Scan requests_dir for existing request_<N>.json files and return
// max(N) + 1, or 0 if none exist. This is how the C++ side picks the id for
// the next request, so request ids are unique and monotonically increasing
// across the lifetime of requests_dir (they are never reused, even if old
// request_<N>.json files are deleted).
inline int nextRequestNumber(const std::string& requests_dir)
{
    int next = 0;
    if (!fs::exists(requests_dir)) return next;

    const std::string prefix = "request_";
    const std::string suffix = ".json";
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

// Write request_<N>.json to requests_dir, where N = nextRequestNumber(requests_dir).
// Returns N so the caller can wait for the matching reconfiguration_<N>.json
// and tag convergence_data/conv_<N>.json with the same id.
//
// Format:
// {
//   "n_blocks": nb,
//   "block_sizes": [...],
//   "krylov_dim": m,
//   "ranks": [...],
//   "C_hat": { "0_0": [...], "0_1": [...], ... },
//   "b_hat": { "0": [...], "1": [...], ... }
// }
//
// Every entry scales with nb (number of blocks) and ranks[j] (Krylov-derived
// rank, <= krylov_dim) only — nothing here scales with the global problem
// size n, so this stays small even for very large systems.
inline int writeRequestJson(const CHatData& chat, const std::string& requests_dir)
{
    const int request_id = nextRequestNumber(requests_dir);

    // Write to a temp file and atomically rename into place. Without this, a
    // reader's fs::exists(path) check (in wait_for_request.py) can observe
    // the file mid-write and parse incomplete JSON. A same-directory
    // rename() is atomic on POSIX, so readers only ever see the file fully
    // written or not at all.
    const std::string path     = requests_dir + "/request_" + std::to_string(request_id) + ".json";
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

    // C_hat blocks
    ofs << "  \"C_hat\": {\n";
    {
        bool first = true;
        for (int i = 0; i < chat.nb; ++i) {
            for (int j = 0; j < chat.nb; ++j) {
                if (!first) ofs << ",\n";
                first = false;
                ofs << "    \"" << i << "_" << j << "\":\n";
                writeDenseJson(ofs, chat.blocks[i][j], 4);
            }
        }
    }
    ofs << "\n  },\n";

    // b_hat vectors
    ofs << "  \"b_hat\": {\n";
    for (int j = 0; j < chat.nb; ++j) {
        ofs << "    \"" << j << "\": ";
        writeVectorJson(ofs, chat.b_hat[j]);
        if (j + 1 < chat.nb) ofs << ",";
        ofs << "\n";
    }
    ofs << "  }\n";
    ofs << "}\n";
    ofs.close();

    fs::rename(tmp_path, path);

    std::cout << "[TekoAdaptive] wrote " << path << "\n";
    return request_id;
}

// Poll for reconfiguration_<request_id>.json, read ordering vector, delete
// file. Returns empty vector on timeout.
inline std::vector<int> waitForOrdering(
    const std::string& requests_dir,
    int                request_id,
    int                timeout_s = 600)
{
    const std::string path = requests_dir + "/reconfiguration_" + std::to_string(request_id) + ".json";
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
            std::remove(path.c_str());

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
};

// Build a new Teko block preconditioner from an ordering vector.
//
// ordering[k] = g means old block k maps to new group g.
// Groups with a single member: singleton.
// Groups with multiple members: merged (block-diagonal sub-inverse).
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
    std::vector<Teko::LinearOp> newInvDiag(nb_new);
    for (int g = 0; g < nb_new; ++g) {
        const auto& mem = groups[g];
        if (mem.size() == 1) {
            // Singleton: factorize A_blocked(b,b) directly.
            newInvDiag[g] =
                Teko::buildInverse(*invFact,
                                   A_blocked->getBlock(mem[0], mem[0]));
        } else {
            // Merged group: block-diagonal inverse of sub-blocks.
            // Sub-operator = sub-blocked diagonal (len(mem) × len(mem)).
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

            std::vector<Teko::LinearOp> subInvDiag(mem.size());
            for (int si = 0; si < (int)mem.size(); ++si)
                subInvDiag[si] = Teko::buildInverse(*invFact,
                                                     A_blocked->getBlock(mem[si], mem[si]));

            newInvDiag[g] = Teko::createBlockDiagonalInverseOp(subTeko, subInvDiag);
        }
    }

    // ── Assemble the preconditioner ────────────────────────────────────────
    Teko::LinearOp precOp;
    if (method == "jacobi")
        precOp = Teko::createBlockDiagonalInverseOp(tekoNew, newInvDiag);
    else
        precOp = Teko::createBlockLowerTriInverseOp(tekoNew, newInvDiag);

    return ReconfiguredSystem{precOp, tekoNew, mgr};
}

// ═══════════════════════════════════════════════════════════════════════════
// Section 5: orchestration — the hook registered with BelosAdaptiveHook
// ═══════════════════════════════════════════════════════════════════════════

// Default location for request_<N>.json / reconfiguration_<N>.json /
// convergence_data/ when TEKO_RECONFIG_REQUESTS_DIR is unset. Hardcoded for
// now, relative to the Trilinos checkout this package lives in.
constexpr const char* kDefaultRequestsDir = "/workspace/Trilinos/teko-reconfig";

// adaptiveLoop()'s Phase 4 below runs a second flexible FGMRES solve, which
// (if it converges) would re-invoke this same hook recursively — and that
// recursive call operates on the *reordered* operator, whose merged-group
// blocks are nested ProductMultiVectors rather than plain TpetraMultiVectors,
// which extractBlockDense() cannot handle. g_inAdaptiveLoop / ScopedGuard
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
// convergence_data/conv_<N>.json by writeConvergenceJson().
struct SolveStats {
    int    iterations;
    double initial_residual;
    double final_residual;
    double wall_time_sec;
};

// Write convergence_data/conv_<request_id>.json (atomically, via a .tmp file
// + rename). s2 == nullptr means no second solve was attempted (Phase 2 timed
// out waiting for reconfiguration_<request_id>.json), and is recorded as
// "solve2": null.
inline void writeConvergenceJson(
    const std::string& convergence_dir,
    int                 request_id,
    const SolveStats&   s1,
    const SolveStats*   s2)
{
    fs::create_directories(convergence_dir);

    const std::string path     = convergence_dir + "/conv_" + std::to_string(request_id) + ".json";
    const std::string tmp_path = path + ".tmp";
    std::ofstream ofs(tmp_path);
    TEUCHOS_TEST_FOR_EXCEPTION(!ofs.is_open(), std::runtime_error,
        "Teko::KrylovSurrogate::writeConvergenceJson: cannot open " + tmp_path);

    auto writeStats = [&](const char* name, const SolveStats& s, bool trailing_comma) {
        ofs << "  \"" << name << "\": {\n";
        ofs << "    \"iterations\": " << s.iterations << ",\n";
        ofs << "    \"initial_residual\": " << std::setprecision(17) << s.initial_residual << ",\n";
        ofs << "    \"final_residual\": " << std::setprecision(17) << s.final_residual << ",\n";
        ofs << "    \"wall_time_sec\": " << std::setprecision(17) << s.wall_time_sec << "\n";
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

    fs::create_directories(requests_dir);
    const std::string convergence_dir =
        (fs::path(requests_dir).parent_path() / "convergence_data").string();

    // Determine block structure from A_blocked.
    const int nb = A_blocked->productRange()->numBlocks();
    std::vector<int> block_sizes(nb);
    for (int j = 0; j < nb; ++j)
        block_sizes[j] =
            static_cast<int>(A_blocked->productRange()->getBlock(j)->dim());

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

    SolveStats s1{state.curDim, b_norm, solve1_metrics.achieved_tol, solve1_metrics.wall_time_sec};

    // ── Phase 1: compute C_hat, b_hat and write request_<N>.json ──────────
    CHatData chat = computeCHat(state, A_blocked, problem->getRHS(), block_sizes);
    const int request_id = writeRequestJson(chat, requests_dir);

    // ── Phase 2: wait for reconfiguration_<N>.json ─────────────────────────
    std::vector<int> ordering = waitForOrdering(requests_dir, request_id);
    if (ordering.empty()) {
        std::cerr << "[TekoAdaptive] no ordering received; returning first "
                     "solve result.\n";
        writeConvergenceJson(convergence_dir, request_id, s1, nullptr);
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
    const double solve2_wall_time_sec = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solve2Start).count();

    SolveStats s2{solver2.getNumIters(), b_norm, solver2.achievedTol(), solve2_wall_time_sec};
    writeConvergenceJson(convergence_dir, request_id, s1, &s2);

    if (ret2 != Belos::Converged)
        std::cerr << "[TekoAdaptive] second solve did not converge ("
                  << solver2.getNumIters() << " iters). Using first result.\n";
    else {
        // Overwrite the original LHS (x_thyra in teko_ext.cpp) with the
        // result of the second solve.
        Thyra::assign(problem->getLHS().ptr(), *x_new);
        std::cout << "[TekoAdaptive] second solve converged ("
                  << solver2.getNumIters() << " iters).\n";
    }
    // env_guard restores TEKO_RECONFIG_REQUESTS_DIR here (or during unwinding
    // if solver2.solve() threw).
}

} // namespace KrylovSurrogate
} // namespace Teko
