//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOSGRAPH_EXPLICIT_COARSEN_IMPL_HPP
#define KOKKOSGRAPH_EXPLICIT_COARSEN_IMPL_HPP

namespace KokkosGraph {
namespace Impl {

template <typename lno_t, typename size_type, typename device_t, typename fine_rowmap_t, typename fine_entries_t,
          typename labels_t, typename coarse_rowmap_t, typename coarse_entries_t, typename ordinal_view_t>
struct ExplicitGraphCoarsening {
  using exec_space     = typename device_t::execution_space;
  using range_pol      = Kokkos::RangePolicy<exec_space>;
  using team_pol       = Kokkos::TeamPolicy<exec_space>;
  using team_member_t  = typename team_pol::member_type;
  using bitset_t       = Kokkos::Bitset<device_t>;
  using const_bitset_t = Kokkos::ConstBitset<device_t>;

  struct ClusterSizeFunctor {
    ClusterSizeFunctor(const ordinal_view_t& counts_, const labels_t& vertClusters_)
        : counts(counts_), vertClusters(vertClusters_) {}
    KOKKOS_INLINE_FUNCTION void operator()(const lno_t i) const { Kokkos::atomic_increment(&counts(vertClusters(i))); }
    ordinal_view_t counts;
    labels_t vertClusters;
  };

  struct FillClusterVertsFunctor {
    FillClusterVertsFunctor(const ordinal_view_t& clusterOffsets_, const ordinal_view_t& clusterVerts_,
                            const labels_t& vertClusters_, const ordinal_view_t& insertCounts_)
        : clusterOffsets(clusterOffsets_),
          clusterVerts(clusterVerts_),
          vertClusters(vertClusters_),
          insertCounts(insertCounts_) {}
    KOKKOS_INLINE_FUNCTION void operator()(const lno_t i) const {
      lno_t cluster        = vertClusters(i);
      lno_t offset         = clusterOffsets(cluster) + Kokkos::atomic_fetch_add(&insertCounts(cluster), 1);
      clusterVerts(offset) = i;
    }
    ordinal_view_t clusterOffsets;
    ordinal_view_t clusterVerts;
    labels_t vertClusters;
    ordinal_view_t insertCounts;
  };

  struct BuildCrossClusterMaskFunctor {
    BuildCrossClusterMaskFunctor(const fine_rowmap_t& rowmap_, const fine_entries_t& colinds_,
                                 const ordinal_view_t& clusterOffsets_, const ordinal_view_t& clusterVerts_,
                                 const labels_t& vertClusters_, const bitset_t& mask_)
        : numRows(rowmap_.extent(0) - 1),
          rowmap(rowmap_),
          colinds(colinds_),
          clusterOffsets(clusterOffsets_),
          clusterVerts(clusterVerts_),
          vertClusters(vertClusters_),
          mask(mask_) {}

    // Used a fixed-size hash set in shared memory
    KOKKOS_INLINE_FUNCTION constexpr int tableSize() const {
      // Should always be a power-of-two, so that X % tableSize() reduces to a
      // bitwise and.
      return 512;
    }

    // Given a cluster index, get the hash table index.
    // This is the 32-bit xorshift RNG, but it works as a hash function.
    KOKKOS_INLINE_FUNCTION unsigned xorshiftHash(lno_t cluster) const {
      unsigned x = cluster;
      x ^= x << 13;
      x ^= x >> 17;
      x ^= x << 5;
      return x;
    }

    KOKKOS_INLINE_FUNCTION bool lookup(lno_t cluster, int* table) const {
      unsigned h = xorshiftHash(cluster);
      for (unsigned i = h; i < h + 2; i++) {
        if (table[i % tableSize()] == cluster) return true;
      }
      return false;
    }

    // Try to insert the edge between cluster (team's cluster) and neighbor
    // (neighboring cluster) by inserting nei into the table.
    KOKKOS_INLINE_FUNCTION bool insert(lno_t cluster, lno_t nei, int* table) const {
      unsigned h = xorshiftHash(nei);
      for (unsigned i = h; i < h + 2; i++) {
        if (Kokkos::atomic_compare_exchange_strong<int>(&table[i % tableSize()], cluster, nei)) return true;
      }
      return false;
    }

    KOKKOS_INLINE_FUNCTION void operator()(const team_member_t t) const {
      lno_t cluster     = t.league_rank();
      lno_t clusterSize = clusterOffsets(cluster + 1) - clusterOffsets(cluster);
      // Use a fixed-size hash table per thread to accumulate neighbor of the
      // cluster. If it fills up (very unlikely) then just count every remaining
      // edge going to another cluster not already in the table; this provides a
      // reasonable upper bound for overallocating the cluster graph. each
      // thread handles a cluster
      int* table = (int*)t.team_shmem().get_shmem(tableSize() * sizeof(int));
      // mark every entry as cluster (self-loop) to represent free/empty
      Kokkos::parallel_for(Kokkos::TeamVectorRange(t, tableSize()), [&](const lno_t i) { table[i] = cluster; });
      t.team_barrier();
      // now, for each row belonging to the cluster, iterate through the
      // neighbors
      Kokkos::parallel_for(Kokkos::TeamThreadRange(t, clusterSize), [&](const lno_t i) {
        lno_t row    = clusterVerts(clusterOffsets(cluster) + i);
        lno_t rowDeg = rowmap(row + 1) - rowmap(row);
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(t, rowDeg), [&](const lno_t j) {
          lno_t nei = colinds(rowmap(row) + j);
          // Remote neighbors are not included
          if (nei >= numRows) return;
          lno_t neiCluster = vertClusters(nei);
          if (neiCluster != cluster) {
            // Have a neighbor. Try to find it in the
            // table.
            if (!lookup(neiCluster, table)) {
              // Not in the table. Try to insert it.
              insert(cluster, neiCluster, table);
              // Whether or not insertion succeeded,
              // this is a cross-cluster edge possibly
              // not seen before
              mask.set(rowmap(row) + j);
            }
          }
        });
      });
    }

    size_t team_shmem_size(int /*teamSize*/) const { return tableSize() * sizeof(int); }

    lno_t numRows;
    fine_rowmap_t rowmap;
    fine_entries_t colinds;
    ordinal_view_t clusterOffsets;
    ordinal_view_t clusterVerts;
    labels_t vertClusters;
    bitset_t mask;
  };

  struct FillClusterEntriesFunctor {
    FillClusterEntriesFunctor(const fine_rowmap_t& rowmap_, const fine_entries_t& colinds_,
                              const coarse_rowmap_t& clusterRowmap_, const coarse_entries_t& clusterEntries_,
                              const ordinal_view_t& clusterOffsets_, const ordinal_view_t& clusterVerts_,
                              const labels_t& vertClusters_, const bitset_t& edgeMask_)
        : rowmap(rowmap_),
          colinds(colinds_),
          clusterRowmap(clusterRowmap_),
          clusterEntries(clusterEntries_),
          clusterOffsets(clusterOffsets_),
          clusterVerts(clusterVerts_),
          vertClusters(vertClusters_),
          edgeMask(edgeMask_) {}
    // Run this scan over entries in clusterVerts (reordered point rows)
    KOKKOS_INLINE_FUNCTION void operator()(const lno_t i, lno_t& lcount, const bool& finalPass) const {
      lno_t numRows      = rowmap.extent(0) - 1;
      lno_t row          = clusterVerts(i);
      size_type rowStart = rowmap(row);
      size_type rowEnd   = rowmap(row + 1);
      lno_t cluster      = vertClusters(row);
      lno_t clusterStart = clusterOffsets(cluster);
      // Count the number of entries in this row.
      // This is how much lcount will be increased by,
      // yielding the offset corresponding to
      // these point entries in the cluster entries.
      lno_t rowEntries = 0;
      for (size_type j = rowStart; j < rowEnd; j++) {
        if (edgeMask.test(j)) rowEntries++;
      }
      if (finalPass) {
        // if this is the last row in the cluster, update the upper bound in
        // clusterRowmap
        if (i == clusterStart) {
          clusterRowmap(cluster) = lcount;
        }
        lno_t clusterEdge = lcount;
        // populate clusterEntries for these edges
        for (size_type j = rowStart; j < rowEnd; j++) {
          if (edgeMask.test(j)) {
            clusterEntries(clusterEdge++) = vertClusters(colinds(j));
          }
        }
      }
      // update the scan result at the end (exclusive)
      lcount += rowEntries;
      if (i == numRows - 1 && finalPass) {
        // on the very last row, set the last entry of the cluster rowmap
        clusterRowmap(clusterRowmap.extent(0) - 1) = lcount;
      }
    }
    fine_rowmap_t rowmap;
    fine_entries_t colinds;
    coarse_rowmap_t clusterRowmap;
    coarse_entries_t clusterEntries;
    ordinal_view_t clusterOffsets;
    ordinal_view_t clusterVerts;
    labels_t vertClusters;
    const_bitset_t edgeMask;
  };

  // Constructor just does the computation and outputs to coarseRowmap,
  // coarseEntries.
  ExplicitGraphCoarsening(const fine_rowmap_t& fineRowmap, const fine_entries_t& fineEntries, const labels_t& labels,
                          lno_t numCoarseVerts) {
    lno_t numFineVerts = fineRowmap.extent(0);
    if (numFineVerts <= 1) {
      coarseRowmap  = coarse_rowmap_t();
      coarseEntries = coarse_entries_t();
      return;
    }
    numFineVerts--;
    clusterOffsets = ordinal_view_t("Cluster offsets", numCoarseVerts + 1);
    clusterVerts   = ordinal_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Cluster verts"), numFineVerts);
    Kokkos::parallel_for(range_pol(0, numFineVerts), ClusterSizeFunctor(clusterOffsets, labels));
    KokkosKernels::Impl::exclusive_parallel_prefix_sum<ordinal_view_t, exec_space>(numCoarseVerts + 1, clusterOffsets);
    {
      ordinal_view_t tempInsertCounts("Temporary cluster insert counts", numCoarseVerts);
      Kokkos::parallel_for(range_pol(0, numFineVerts),
                           FillClusterVertsFunctor(clusterOffsets, clusterVerts, labels, tempInsertCounts));
    }
    // Determine the set of edges (in the point graph) that cross between two
    // distinct clusters
    int vectorSize = KokkosKernels::Impl::kk_get_suggested_vector_size(
        numFineVerts, fineEntries.extent(0), KokkosKernels::Impl::kk_get_exec_space_type<exec_space>());
    bitset_t crossClusterEdgeMask(fineEntries.extent(0));
    size_type numClusterEdges;
    {
      BuildCrossClusterMaskFunctor buildEdgeMask(fineRowmap, fineEntries, clusterOffsets, clusterVerts, labels,
                                                 crossClusterEdgeMask);
      int sharedPerTeam =
          buildEdgeMask.team_shmem_size(0);  // using team-size = 0 for since no per-thread shared is used.
      int teamSize =
          KokkosKernels::Impl::get_suggested_team_size<team_pol>(buildEdgeMask, vectorSize, sharedPerTeam, 0);
      Kokkos::parallel_for(
          team_pol(numCoarseVerts, teamSize, vectorSize).set_scratch_size(0, Kokkos::PerTeam(sharedPerTeam)),
          buildEdgeMask);
      numClusterEdges = crossClusterEdgeMask.count();
    }
    coarseRowmap =
        coarse_rowmap_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Cluster graph rowmap"), numCoarseVerts + 1);
    coarseEntries =
        coarse_entries_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "Cluster graph colinds"), numClusterEdges);
    Kokkos::parallel_scan(range_pol(0, numFineVerts),
                          FillClusterEntriesFunctor(fineRowmap, fineEntries, coarseRowmap, coarseEntries,
                                                    clusterOffsets, clusterVerts, labels, crossClusterEdgeMask));
  }

  coarse_rowmap_t coarseRowmap;
  coarse_entries_t coarseEntries;
  ordinal_view_t clusterOffsets;
  ordinal_view_t clusterVerts;
};

}  // namespace Impl
}  // namespace KokkosGraph

#endif
