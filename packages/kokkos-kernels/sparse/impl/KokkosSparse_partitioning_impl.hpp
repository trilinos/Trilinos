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

#include "KokkosKernels_Utils.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_Random.hpp>
#include "KokkosBlas1_fill.hpp"
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosKernels_Uniform_Initialized_MemoryPool.hpp"

#ifndef _KOKKOS_PARTITIONING_IMP_HPP
#define _KOKKOS_PARTITIONING_IMP_HPP

namespace KokkosSparse {

namespace Impl {

// Fill a view such that v(i) = i
// Does the same thing as std::iota(begin, end)
template <typename View, typename Ordinal>
struct IotaFunctor {
  IotaFunctor(View& v_) : v(v_) {}
  KOKKOS_INLINE_FUNCTION void operator()(const Ordinal i) const { v(i) = i; }
  View v;
};

template <typename HandleType, typename lno_row_view_t, typename lno_nnz_view_t>
struct BalloonClustering {
  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename HandleType::size_type size_type;
  typedef typename HandleType::nnz_lno_t nnz_lno_t;

  typedef typename lno_row_view_t::value_type offset_t;

  typedef typename HandleType::row_lno_temp_work_view_t row_lno_temp_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;
  typedef
      typename HandleType::row_lno_persistent_work_host_view_t row_lno_persistent_work_host_view_t;  // Host view type

  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;
  typedef
      typename HandleType::nnz_lno_persistent_work_host_view_t nnz_lno_persistent_work_host_view_t;  // Host view type

  typedef nnz_lno_persistent_work_view_t nnz_view_t;
  typedef Kokkos::View<float*, MyPersistentMemorySpace> float_view_t;
  // typedef Kokkos::View<nnz_lno_t, MyTempMemorySpace, Kokkos::MemoryTraits<0>>
  // single_view_t; typedef Kokkos::View<nnz_lno_t, Kokkos::HostSpace,
  // Kokkos::MemoryTraits<0>> single_view_host_t;

  typedef Kokkos::RangePolicy<MyExecSpace> my_exec_space;
  typedef Kokkos::Bitset<MyExecSpace> bitset_t;

  typedef Kokkos::RangePolicy<MyExecSpace> range_policy_t;
  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t;
  typedef typename team_policy_t::member_type team_member_t;

  BalloonClustering(size_type numRows_, const lno_row_view_t& rowmap_, const lno_nnz_view_t& colinds_)
      : numRows(numRows_), rowmap(rowmap_), colinds(colinds_), randPool(0xDEADBEEF) {}

  nnz_lno_t numRows;
  lno_row_view_t rowmap;
  lno_nnz_view_t colinds;

  typedef Kokkos::Random_XorShift64_Pool<MyExecSpace> RandPool;
  RandPool randPool;

  struct InitRootsTag {};   // select roots; set their distances to 0
  struct RandomFillTag {};  // assign non-roots to random clusters, and assign large random distances
  struct UpdatePressureTag {};
  struct BalloonTag {};  // run the "balloon" procedure, where each cluster
                         // tries to inflate up to clusterSize

  struct BalloonFunctor {
    BalloonFunctor(const nnz_view_t& vertClusters_, const nnz_view_t& clusterCounts_, const nnz_view_t& distances_,
                   const lno_row_view_t& row_map_, const lno_nnz_view_t& col_inds_, const float_view_t& pressure_,
                   nnz_lno_t clusterSize_, RandPool& randPool_)
        : vertClusters(vertClusters_),
          clusterCounts(clusterCounts_),
          distances(distances_),
          row_map(row_map_),
          col_inds(col_inds_),
          pressure(pressure_),
          clusterSize(clusterSize_),
          numRows(row_map.extent(0) - 1),
          vertLocks(numRows),
          randPool(randPool_) {
      numClusters    = (numRows + clusterSize - 1) / clusterSize;
      avgClusterSize = (double)numRows / numClusters;
      iter           = 0;
    }

    // Run init version over the number of clusters.
    KOKKOS_INLINE_FUNCTION void operator()(const InitRootsTag, const nnz_lno_t i) const {
      nnz_lno_t root;
      auto state = randPool.get_state();
      do {
        root = state.rand(numRows);
      } while (!Kokkos::atomic_compare_exchange_strong(&vertClusters(root), numClusters, i));
      randPool.free_state(state);
      distances(root) = 0;
      pressure(root)  = 1;
    }

    KOKKOS_INLINE_FUNCTION void operator()(const RandomFillTag, const nnz_lno_t i) const {
      if (vertClusters(i) == numClusters) {
        auto state        = randPool.get_state();
        nnz_lno_t cluster = state.rand(numClusters);
        randPool.free_state(state);
        vertClusters(i) = cluster;
        Kokkos::atomic_increment(&clusterCounts(cluster));
        distances(i) = numRows;
        pressure(i)  = 0.1;
      }
    };

    KOKKOS_INLINE_FUNCTION void operator()(const UpdatePressureTag, const nnz_lno_t i) const {
      nnz_lno_t cluster = vertClusters(i);
      if (cluster == numClusters) {
        // unassigned vertices have 0 pressure
        return;
      }
      // count the number of neighbors in the same cluster
      nnz_lno_t sameClusterNeighbors = 0;
      for (size_type j = row_map(i); j < row_map(i + 1); j++) {
        nnz_lno_t nei = col_inds(j);
        if (nei < numRows && nei != i && vertClusters(nei) == cluster) {
          // while we're at it, minimize distance to root as in Djikstra's
          if (distances(nei) + 1 < distances(i)) distances(i) = distances(nei) + 1;
          sameClusterNeighbors++;
        }
      }
      nnz_lno_t curSize = clusterCounts(cluster);
      // update pressure, if cluster is undersized
      nnz_lno_t shortage   = clusterSize - curSize;
      float pressureChange = (shortage * shortage * (1.0f + 0.2f * sameClusterNeighbors)) / (1.0f + distances(i));
      if (shortage > 0) pressure(i) += pressureChange;
    }

    KOKKOS_INLINE_FUNCTION void operator()(const BalloonTag, const nnz_lno_t i, double& sizeDeviation) const {
      nnz_lno_t cluster = vertClusters(i);
      if (cluster == numClusters) return;
      // find the weakest affinity neighbor
      nnz_lno_t weakNei        = numRows;
      float weakestPressure    = pressure(i);
      nnz_lno_t weakNeiCluster = numClusters;
      for (size_type j = row_map(i); j < row_map(i + 1); j++) {
        nnz_lno_t nei = col_inds(j);
        // to annex another vertex, it must be a non-root in a different cluster
        if (nei < numRows && nei != i && vertClusters(nei) != cluster && pressure(nei) < weakestPressure &&
            distances(nei) != 0) {
          weakNei         = nei;
          weakestPressure = pressure(nei);
          weakNeiCluster  = vertClusters(nei);
        }
      }
      if (weakNei != numRows && clusterCounts(cluster) < clusterSize) {
        // this cluster will take over weakNei
        if (vertLocks.set(i)) {
          if (vertLocks.set(weakNei)) {
            Kokkos::atomic_increment(&clusterCounts(cluster));
            if (weakNeiCluster != numClusters) Kokkos::atomic_decrement(&clusterCounts(weakNeiCluster));
            vertClusters(weakNei) = cluster;
            pressure(i) -= pressure(weakNei);
            pressure(weakNei)  = pressure(i);
            distances(weakNei) = distances(i) + 1;
            vertLocks.reset(weakNei);
          }
          vertLocks.reset(i);
        }
      }
      if (distances(i) == 0) {
        // roots update sizeDeviation on behalf of the cluster
        double deviation = clusterCounts(cluster) - avgClusterSize;
        sizeDeviation += deviation * deviation;
      }
    }

    nnz_view_t vertClusters;
    nnz_view_t clusterCounts;
    nnz_view_t distances;
    // row_map/col_inds of input graph (read-only)
    lno_row_view_t row_map;
    lno_nnz_view_t col_inds;
    float_view_t pressure;
    // constants
    nnz_lno_t clusterSize;
    nnz_lno_t numClusters;
    nnz_lno_t numRows;
    nnz_lno_t iter;
    Kokkos::Bitset<MyExecSpace> vertLocks;
    RandPool randPool;
    double avgClusterSize;
  };

  nnz_view_t run(nnz_lno_t clusterSize) {
    nnz_view_t vertClusters("Vertex cluster labels", numRows);
    // For the sake of completeness, handle the clusterSize = 1 case by
    // generating a trivial (identity) clustering.
    if (clusterSize == 1) {
      Kokkos::parallel_for(Kokkos::RangePolicy<MyExecSpace>(0, numRows),
                           IotaFunctor<nnz_view_t, nnz_lno_t>(vertClusters));
      return vertClusters;
    }
    nnz_lno_t numClusters = (numRows + clusterSize - 1) / clusterSize;
    nnz_view_t distances("Root distances", numRows);
    nnz_view_t clusterCounts("Vertices per cluster", numClusters);
    float_view_t pressure("Cluster pressure", numRows);
    Kokkos::deep_copy(clusterCounts, 1);
    Kokkos::deep_copy(vertClusters, (nnz_lno_t)numClusters);
    Kokkos::deep_copy(distances, numRows);
    BalloonFunctor funct(vertClusters, clusterCounts, distances, rowmap, colinds, pressure, clusterSize, randPool);
    Kokkos::Timer globalTimer;
    Kokkos::Timer timer;
    timer.reset();
    Kokkos::parallel_for(Kokkos::RangePolicy<MyExecSpace, InitRootsTag>(0, numClusters), funct);
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    MyExecSpace().fence();
    std::cout << "Creating roots: " << timer.seconds() << '\n';
    timer.reset();
#endif
    double stoppingRMS = sqrt(numClusters * (0.02 * clusterSize) * (0.02 * clusterSize));
    double deviation   = (double)numClusters * (clusterSize - 1) * (clusterSize - 1);
    int regressions    = 0;
    while (true) {
      Kokkos::parallel_for(Kokkos::RangePolicy<MyExecSpace, UpdatePressureTag>(0, numRows), funct);
      double iterDeviation = 0;
      Kokkos::parallel_reduce(Kokkos::RangePolicy<MyExecSpace, BalloonTag>(0, numRows), funct,
                              Kokkos::Sum<double>(iterDeviation));
      if (iterDeviation <= stoppingRMS || iterDeviation == deviation) {
        // got within 2% RMS of optimal, or stagnated
        deviation = iterDeviation;
        break;
      } else if (iterDeviation >= deviation) {
        regressions++;
        if (regressions == 3) {
          deviation = iterDeviation;
          break;
        }
      }
      deviation = iterDeviation;
      funct.iter++;
    }
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    MyExecSpace().fence();
    std::cout << "Expanding clusters for " << funct.iter << " iterations: " << timer.seconds() << '\n';
    timer.reset();
#endif
    Kokkos::parallel_for(Kokkos::RangePolicy<MyExecSpace, RandomFillTag>(0, numRows), funct);
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    MyExecSpace().fence();
    std::cout << "Randomly assigning clusters to remaining: " << timer.seconds() << '\n';
    std::cout << "Clustering total: " << globalTimer.seconds() << "\n\n";
#endif
    return vertClusters;
  }
};

}  // namespace Impl
}  // namespace KokkosSparse

#endif
