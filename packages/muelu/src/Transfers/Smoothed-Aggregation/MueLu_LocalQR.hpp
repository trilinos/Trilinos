#ifndef MUELU_LOCALQR_HPP
#define MUELU_LOCALQR_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "Xpetra_ConfigDefs.hpp"

#include "KokkosBlas1_set.hpp"
#include "KokkosBatched_QR_FormQ_TeamVector_Internal.hpp"
#include "KokkosBatched_ApplyQ_Decl.hpp"
#include "KokkosBatched_SetIdentity_Decl.hpp"
#include "KokkosBatched_SetIdentity_Impl.hpp"
#include "Kokkos_DualView.hpp"
#include "Kokkos_Pair.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "KokkosBatched_QR_Decl.hpp"
#include "KokkosBatched_QR_Serial_Impl.hpp"
#include "KokkosBatched_QR_TeamVector_Impl.hpp"
#include "KokkosBatched_QR_FormQ_TeamVector_Internal.hpp"

// #define MUELU_DEBUG_QR 1

namespace MueLu::LocalQR {

template <class LocalOrdinal, class View>
class ReduceMaxFunctor {
 public:
  ReduceMaxFunctor(View view)
    : view_(view) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const LocalOrdinal& i, LocalOrdinal& vmax) const {
    if (vmax < view_(i))
      vmax = view_(i);
  }

  KOKKOS_INLINE_FUNCTION
  void join(LocalOrdinal& dst, const LocalOrdinal& src) const {
    if (dst < src) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(LocalOrdinal& dst) const {
    dst = 0;
  }

 private:
  View view_;
};

// local QR decomposition
// local QR decomposition
template <class LOType, class GOType, class SCType, class DeviceType, class NspType, class aggRowsType, class maxAggDofSizeType, class agg2RowMapLOType, class statusType, class rowsType, class rowsAuxType, class colsAuxType, class valsAuxType>
class LocalQRDecompFunctor {
 private:
  using LO = LOType;
  using GO = GOType;
  using SC = SCType;

  using execution_space = typename DeviceType::execution_space;
  using impl_SC         = typename Kokkos::ArithTraits<SC>::val_type;
  using impl_ATS        = Kokkos::ArithTraits<impl_SC>;
  using Magnitude       = typename impl_ATS::magnitudeType;

  using shared_matrix = Kokkos::View<impl_SC**, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;
  using shared_vector = Kokkos::View<impl_SC*, typename execution_space::scratch_memory_space, Kokkos::MemoryUnmanaged>;

  using temp_view_type1 = Kokkos::View<impl_SC*, execution_space>;
  using temp_view_type  = Kokkos::View<impl_SC**, execution_space>;

  NspType fineNS;
  NspType coarseNS;
  aggRowsType aggRows;
  maxAggDofSizeType maxAggDofSize;  //< maximum number of dofs in aggregate (max size of aggregate * numDofsPerNode)
  agg2RowMapLOType agg2RowMapLO;
  statusType statusAtomic;
  rowsType rows;
  rowsAuxType rowsAux;
  colsAuxType colsAux;
  valsAuxType valsAux;
  bool doQRStep;
  temp_view_type1 tau;
  temp_view_type work;

 public:
  LocalQRDecompFunctor(NspType fineNS_, NspType coarseNS_, aggRowsType aggRows_, maxAggDofSizeType maxAggDofSize_, agg2RowMapLOType agg2RowMapLO_, statusType statusAtomic_, rowsType rows_, rowsAuxType rowsAux_, colsAuxType colsAux_, valsAuxType valsAux_, bool doQRStep_)
    : fineNS(fineNS_)
    , coarseNS(coarseNS_)
    , aggRows(aggRows_)
    , maxAggDofSize(maxAggDofSize_)
    , agg2RowMapLO(agg2RowMapLO_)
    , statusAtomic(statusAtomic_)
    , rows(rows_)
    , rowsAux(rowsAux_)
    , colsAux(colsAux_)
    , valsAux(valsAux_)
    , doQRStep(doQRStep_) {
    work = temp_view_type("work", aggRows_.extent(0), maxAggDofSize_);
    tau  = temp_view_type1("tau", fineNS.extent(0));
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<execution_space>::member_type& thread, size_t& nnz) const {
    using member_type = typename Kokkos::TeamPolicy<execution_space>::member_type;

    auto agg = thread.league_rank();

    const auto aggOffset = aggRows(agg);
    // size of aggregate: number of DOFs in aggregate
    const auto aggSize = aggRows(agg + 1) - aggOffset;

    const impl_SC one  = impl_ATS::one();
    const impl_SC zero = impl_ATS::zero();

    const int m = aggSize;
    const int n = fineNS.extent(1);

    // calculate row offset for coarse nullspace
    Xpetra::global_size_t offset = agg * n;

    if (doQRStep) {
      // A is m x n
      // Q is m x m
      // R is m x n

      // A (initially) gets overwritten with R in QR
      shared_matrix r(thread.team_shmem(), m, n);
      // Q
      shared_matrix q(thread.team_shmem(), m, m);

      // Extract the piece of the nullspace corresponding to the aggregate
      for (int j = 0; j < n; j++)
        for (int k = 0; k < m; k++)
          r(k, j) = fineNS(agg2RowMapLO(aggOffset + k), j);

#ifdef MUELU_DEBUG_QR
      Kokkos::printf("\nA %d x %d\n", m, n);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
          Kokkos::printf(" %5.3lf ", r(i, j));
        Kokkos::printf("\n");
      }
#endif
      thread.team_barrier();

      if (m >= n) {
        // team_tau has size n
        auto team_tau = Kokkos::subview(tau, Kokkos::make_pair(aggOffset, aggOffset + n));
        KokkosBlas::TeamVectorSet<member_type>::invoke(thread, zero, team_tau);

        // team_work has size m
        auto team_work = Kokkos::subview(work, agg, Kokkos::make_pair(0, m));
        KokkosBlas::TeamVectorSet<member_type>::invoke(thread, zero, team_work);

        // Calculate QR
        KokkosBatched::TeamVectorQR<member_type, KokkosBlas::Algo::QR::Unblocked>::invoke(thread, r, team_tau, team_work);
        thread.team_barrier();

#ifdef MUELU_DEBUG_QR
        Kokkos::printf("\nQR %d x %d\n", m, n);
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < n; j++)
            Kokkos::printf(" %5.3lf ", r(i, j));
          Kokkos::printf("\n");
        }
#endif

        // Initialize Q to identity
        KokkosBatched::TeamSetIdentity<member_type>::invoke(thread, q);
        thread.team_barrier();

        // Form Q
        KokkosBatched::TeamVectorApplyQ<member_type, KokkosBatched::Side::Left, KokkosBlas::Trans::NoTranspose, KokkosBlas::Algo::ApplyQ::Unblocked>::invoke(thread, r, team_tau, q, team_work);
        thread.team_barrier();

#ifdef MUELU_DEBUG_QR
        Kokkos::printf("\nQ1 %d x %d\n", m, m);
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++)
            Kokkos::printf(" %5.3lf ", q(i, j));
          Kokkos::printf("\n");
        }

        // Initialize Q to identity
        KokkosBatched::TeamSetIdentity<member_type>::invoke(thread, q);
        thread.team_barrier();

        for (int k = 0; k < m; ++k) {
          auto q_col = Kokkos::subview(q, Kokkos::ALL(), k);
          KokkosBatched::TeamVectorApplyQ<member_type, KokkosBatched::Side::Left, KokkosBlas::Trans::NoTranspose, KokkosBlas::Algo::ApplyQ::Unblocked>::invoke(thread, r, team_tau, q_col, team_work);
        }

        Kokkos::printf("\nQ2 %d x %d\n", m, m);
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++)
            Kokkos::printf(" %5.3lf ", q(i, j));
          Kokkos::printf("\n");
        }

        Kokkos::printf("\ntau %d \n", n);
        for (int i = 0; i < n; i++) {
          Kokkos::printf(" %5.3lf ", team_tau(i));
          Kokkos::printf("\n");
        }

        Kokkos::printf("\work %d \n", team_work.extent(0));
        for (int i = 0; i < team_work.extent(0); i++) {
          Kokkos::printf(" %5.3lf ", team_work(i));
          Kokkos::printf("\n");
        }
#endif

        // Build coarse nullspace using the upper triangular part of R
        for (int j = 0; j < n; j++) {
          for (int k = 0; k < n; k++)
            coarseNS(offset + k, j) = (k <= j) ? r(k, j) : zero;
        }

      } else {
        // Special handling for m < n (i.e. single node aggregates in structural mechanics)

        // The local QR decomposition is not possible in the "overconstrained"
        // case (i.e. number of columns in qr > number of rowsAux), which
        // corresponds to #DOFs in Aggregate < n. For usual problems this
        // is only possible for single node aggregates in structural mechanics.
        // (Similar problems may arise in discontinuous Galerkin problems...)
        // We bypass the QR decomposition and use an identity block in the
        // tentative prolongator for the single node aggregate and transfer the
        // corresponding fine level null space information 1-to-1 to the coarse
        // level null space part.

        // NOTE: The resulting tentative prolongation operator has
        // (m*DofsPerNode-n) zero columns leading to a singular
        // coarse level operator A.  To deal with that one has the following
        // options:
        // - Use the "RepairMainDiagonal" flag in the RAPFactory (default:
        //   false) to add some identity block to the diagonal of the zero rowsAux
        //   in the coarse level operator A, such that standard level smoothers
        //   can be used again.
        // - Use special (projection-based) level smoothers, which can deal
        //   with singular matrices (very application specific)
        // - Adapt the code below to avoid zero columns. However, we do not
        //   support a variable number of DOFs per node in MueLu/Xpetra which
        //   makes the implementation really hard.
        //
        // FIXME: do we need to check for singularity here somehow? Zero
        // columns would be easy but linear dependency would require proper QR.

        // R = extended (by adding identity rowsAux) qr
        for (int j = 0; j < n; j++)
          for (int k = 0; k < n; k++)
            if (k < m)
              coarseNS(offset + k, j) = r(k, j);
            else
              coarseNS(offset + k, j) = (k == j ? one : zero);

        // Q = I (rectangular)
        for (int i = 0; i < m; i++)
          for (int j = 0; j < n; j++)
            q(i, j) = (j == i ? one : zero);
      }

#ifdef MUELU_DEBUG_QR
      // Kokkos::printf("\nQ3 %d x %d\n", m, m);
      // for (int i = 0; i < m; i++) {
      //   for (int j = 0; j < m; j++)
      //     Kokkos::printf(" %5.3lf ", q(i, j));
      //   Kokkos::printf("\n");
      // }

      Kokkos::printf("\nR %d x %d\n", m, n);
      for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
          Kokkos::printf(" %5.3lf ", i <= j ? r(i, j) : zero);
        Kokkos::printf("\n");
      }

#endif

      // Process each row in the local Q factor and fill helper arrays to assemble P
      for (int j = 0; j < m; j++) {
        LO localRow     = agg2RowMapLO(aggRows(agg) + j);
        size_t rowStart = rowsAux(localRow);
        size_t lnnz     = 0;
        for (int k = 0; k < n; k++) {
          // skip zeros
          if (q(j, k) != zero) {
            colsAux(rowStart + lnnz) = offset + k;
            valsAux(rowStart + lnnz) = q(j, k);
            lnnz++;
          }
        }
        rows(localRow + 1) = lnnz;
        nnz += lnnz;
      }
    } else {
      /////////////////////////////
      //      "no-QR" option     //
      /////////////////////////////
      // Local Q factor is just the fine nullspace support over the current aggregate.
      // Local R factor is the identity.
      // TODO I have not implemented any special handling for aggregates that are too
      // TODO small to locally support the nullspace, as is done in the standard QR
      // TODO case above.

      for (int j = 0; j < m; j++) {
        LO localRow     = agg2RowMapLO(aggRows(agg) + j);
        size_t rowStart = rowsAux(localRow);
        size_t lnnz     = 0;
        for (int k = 0; k < n; k++) {
          const impl_SC qr_jk = fineNS(localRow, k);
          // skip zeros
          if (qr_jk != zero) {
            colsAux(rowStart + lnnz) = offset + k;
            valsAux(rowStart + lnnz) = qr_jk;
            lnnz++;
          }
        }
        rows(localRow + 1) = lnnz;
        nnz += lnnz;
      }

      for (int j = 0; j < n; j++)
        coarseNS(offset + j, j) = one;
    }
  }

  // amount of shared memory
  size_t team_shmem_size(int /* team_size */) const {
    if (doQRStep) {
      int m = maxAggDofSize;
      int n = fineNS.extent(1);
      return shared_matrix::shmem_size(m, n) +  // r
             shared_matrix::shmem_size(m, m);   // q
    } else
      return 0;
  }
};

}  // namespace MueLu::LocalQR

#endif
