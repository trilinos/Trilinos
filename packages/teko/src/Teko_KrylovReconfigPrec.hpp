// Teko_KrylovReconfigPrec.hpp — build the reconfigured ("as-if-original")
// blocked system and its preconditioner from an ordering vector.
//
// An ordering[k] = g maps original block k into new group g. Singleton groups
// reuse the original Tpetra blocks; multi-member groups are assembled into one
// monolithic Tpetra::CrsMatrix (off-diagonal coupling included) and factored
// jointly — exactly the operator the application would have built had it been
// called with this grouping from the start, so the re-solve's factor/iterate
// times are directly comparable to the first solve's. The result
// (ReconfiguredSystem) carries the flat operator, its preconditioner, the group
// membership (to repack vectors), and the factorization wall time.
//
// This is multi-rank aware: each member block keeps its own row distribution
// and only locally-owned rows are inserted. The type aliases (SC/LO/GO/Node,
// MV) come from Teko_KrylovReducedModel.hpp, which this header includes.
#pragma once

// ── Standard library ──────────────────────────────────────────────────────
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <string>
#include <vector>

// ── Teuchos ───────────────────────────────────────────────────────────────
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// ── Thyra ─────────────────────────────────────────────────────────────────
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
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
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_KrylovReducedModel.hpp"  // type aliases (SC/LO/GO/Node, MV, OP)
#include "Teko_Utilities.hpp"

// ── Stratimikos ───────────────────────────────────────────────────────────
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

namespace Teko {
namespace KrylovSurrogate {

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

} // namespace KrylovSurrogate
} // namespace Teko
