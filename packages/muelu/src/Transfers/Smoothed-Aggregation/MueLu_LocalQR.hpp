#ifndef MUELU_LOCALQR_HPP
#define MUELU_LOCALQR_HPP

#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>
#include "Xpetra_ConfigDefs.hpp"

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
    , doQRStep(doQRStep_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<execution_space>::member_type& thread, size_t& nnz) const {
    auto agg = thread.league_rank();

    // size of aggregate: number of DOFs in aggregate
    auto aggSize = aggRows(agg + 1) - aggRows(agg);

    const impl_SC one  = impl_ATS::one();
    const impl_SC two  = one + one;
    const impl_SC zero = impl_ATS::zero();
    const auto zeroM   = impl_ATS::magnitude(zero);

    int m = aggSize;
    int n = fineNS.extent(1);

    // calculate row offset for coarse nullspace
    Xpetra::global_size_t offset = agg * n;

    if (doQRStep) {
      // Extract the piece of the nullspace corresponding to the aggregate
      shared_matrix r(thread.team_shmem(), m, n);  // A (initially), R (at the end)
      for (int j = 0; j < n; j++)
        for (int k = 0; k < m; k++)
          r(k, j) = fineNS(agg2RowMapLO(aggRows(agg) + k), j);
#if 0
          printf("A\n");
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
              printf(" %5.3lf ", r(i,j));
            printf("\n");
          }
#endif

      // Calculate QR decomposition (standard)
      shared_matrix q(thread.team_shmem(), m, m);  // Q
      if (m >= n) {
        bool isSingular = false;

        // Initialize Q^T
        auto qt = q;
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++)
            qt(i, j) = zero;
          qt(i, i) = one;
        }

        for (int k = 0; k < n; k++) {  // we ignore "n" instead of "n-1" to normalize
          // FIXME_KOKKOS: use team
          Magnitude s = zeroM;
          Magnitude norm;
          Magnitude norm_x;
          for (int i = k + 1; i < m; i++)
            s += pow(impl_ATS::magnitude(r(i, k)), 2);
          norm = sqrt(pow(impl_ATS::magnitude(r(k, k)), 2) + s);

          if (norm == zero) {
            isSingular = true;
            break;
          }

          r(k, k) -= norm * one;

          norm_x = sqrt(pow(impl_ATS::magnitude(r(k, k)), 2) + s);
          if (norm_x == zeroM) {
            // We have a single diagonal element in the column.
            // No reflections required. Just need to restor r(k,k).
            r(k, k) = norm * one;
            continue;
          }

          // FIXME_KOKKOS: use team
          for (int i = k; i < m; i++)
            r(i, k) /= norm_x;

          // Update R(k:m,k+1:n)
          for (int j = k + 1; j < n; j++) {
            // FIXME_KOKKOS: use team in the loops
            impl_SC si = zero;
            for (int i = k; i < m; i++)
              si += r(i, k) * r(i, j);
            for (int i = k; i < m; i++)
              r(i, j) -= two * si * r(i, k);
          }

          // Update Q^T (k:m,k:m)
          for (int j = k; j < m; j++) {
            // FIXME_KOKKOS: use team in the loops
            impl_SC si = zero;
            for (int i = k; i < m; i++)
              si += r(i, k) * qt(i, j);
            for (int i = k; i < m; i++)
              qt(i, j) -= two * si * r(i, k);
          }

          // Fix R(k:m,k)
          r(k, k) = norm * one;
          for (int i = k + 1; i < m; i++)
            r(i, k) = zero;
        }

#if 0
            // Q = (Q^T)^T
            for (int i = 0; i < m; i++)
              for (int j = 0; j < i; j++) {
                impl_SC tmp  = qt(i,j);
                qt(i,j) = qt(j,i);
                qt(j,i) = tmp;
              }
#endif

        // Build coarse nullspace using the upper triangular part of R
        for (int j = 0; j < n; j++)
          for (int k = 0; k <= j; k++)
            coarseNS(offset + k, j) = r(k, j);

        if (isSingular) {
          statusAtomic(1) = true;
          return;
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

#if 0
          printf("R\n");
          for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
              printf(" %5.3lf ", coarseNS(i,j));
            printf("\n");
          }

          printf("Q\n");
          for (int i = 0; i < aggSize; i++) {
            for (int j = 0; j < aggSize; j++)
              printf(" %5.3lf ", q(i,j));
            printf("\n");
          }
#endif
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
