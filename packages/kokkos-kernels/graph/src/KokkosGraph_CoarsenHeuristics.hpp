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
#pragma once
// exclude from Cuda builds without lambdas enabled
#if !defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#include <limits>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_Functional.hpp>
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosGraph_MIS2.hpp"

namespace KokkosGraph {

namespace Experimental {

template <class crsMat>
class coarsen_heuristics {
 public:
  // define internal types
  using matrix_t       = crsMat;
  using exec_space     = typename matrix_t::execution_space;
  using mem_space      = typename matrix_t::memory_space;
  using Device         = typename matrix_t::device_type;
  using ordinal_t      = typename matrix_t::ordinal_type;
  using edge_offset_t  = typename matrix_t::size_type;
  using scalar_t       = typename matrix_t::value_type;
  using vtx_view_t     = typename Kokkos::View<ordinal_t*, Device>;
  using wgt_view_t     = typename Kokkos::View<scalar_t*, Device>;
  using edge_view_t    = typename Kokkos::View<edge_offset_t*, Device>;
  using edge_subview_t = typename Kokkos::View<edge_offset_t, Device>;
  using rand_view_t    = typename Kokkos::View<uint64_t*, Device>;
  using graph_type     = typename matrix_t::staticcrsgraph_type;
  using policy_t       = typename Kokkos::RangePolicy<exec_space>;
  using team_policy_t  = typename Kokkos::TeamPolicy<exec_space>;
  using member         = typename team_policy_t::member_type;
  using part_view_t    = typename Kokkos::View<int*, Device>;
  using pool_t         = Kokkos::Random_XorShift64_Pool<Device>;
  using gen_t          = typename pool_t::generator_type;
  using hasher_t       = Kokkos::pod_hash<ordinal_t>;
  static constexpr ordinal_t get_null_val() {
    if (std::is_signed<ordinal_t>::value) {
      return -1;
    } else {
      return std::numeric_limits<ordinal_t>::max();
    }
  }
  static constexpr ordinal_t ORD_MAX = get_null_val();

  static vtx_view_t generate_permutation(ordinal_t n, pool_t rand_pool) {
    rand_view_t randoms("randoms", n);

    Kokkos::parallel_for(
        "create random entries", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          gen_t generator = rand_pool.get_state();
          randoms(i)      = generator.urand64();
          rand_pool.free_state(generator);
        });

    int t_buckets = 2 * n;
    vtx_view_t buckets("buckets", t_buckets);
    Kokkos::parallel_for(
        "init buckets", policy_t(0, t_buckets), KOKKOS_LAMBDA(ordinal_t i) { buckets(i) = ORD_MAX; });

    uint64_t max         = std::numeric_limits<uint64_t>::max();
    uint64_t bucket_size = max / t_buckets;
    Kokkos::parallel_for(
        "insert buckets", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          ordinal_t bucket = randoms(i) / bucket_size;
          // find the next open bucket
          for (;; bucket++) {
            if (bucket >= t_buckets) bucket -= t_buckets;
            if (buckets(bucket) == ORD_MAX) {
              // attempt to insert into bucket
              if (Kokkos::atomic_compare_exchange_strong(&buckets(bucket), ORD_MAX, i)) {
                break;
              }
            }
          }
        });

    vtx_view_t permute("permutation", n);
    Kokkos::parallel_scan(
        "extract permutation", policy_t(0, t_buckets),
        KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
          if (buckets(i) != ORD_MAX) {
            if (final) {
              permute(update) = buckets(i);
            }
            update++;
          }
        });

    return permute;
  }

  // create a mapping when some vertices are already mapped
  // hn is a list of vertices such that vertex i wants to aggregate with vertex
  // hn(i)
  static ordinal_t parallel_map_construct_prefilled(vtx_view_t vcmap, const ordinal_t n, const vtx_view_t vperm,
                                                    const vtx_view_t hn,
                                                    Kokkos::View<ordinal_t, Device> nvertices_coarse) {
    vtx_view_t match("match", n);
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          if (vcmap(i) == ORD_MAX) {
            match(i) = ORD_MAX;
          } else {
            match(i) = n + 1;
          }
        });
    ordinal_t perm_length = vperm.extent(0);

    // construct mapping using heaviest edges
    int swap             = 1;
    vtx_view_t curr_perm = vperm;
    while (perm_length > 0) {
      vtx_view_t next_perm("next perm", perm_length);
      Kokkos::View<ordinal_t, Device> next_length("next_length");

      Kokkos::parallel_for(
          policy_t(0, perm_length), KOKKOS_LAMBDA(ordinal_t i) {
            ordinal_t u   = curr_perm(i);
            ordinal_t v   = hn(u);
            int condition = u < v;
            // need to enforce an ordering condition to allow hard-stall
            // conditions to be broken
            if (condition ^ swap) {
              if (Kokkos::atomic_compare_exchange_strong(&match(u), ORD_MAX, v)) {
                if (u == v || Kokkos::atomic_compare_exchange_strong(&match(v), ORD_MAX, u)) {
                  ordinal_t cv = Kokkos::atomic_fetch_add(&nvertices_coarse(), 1);
                  vcmap(u)     = cv;
                  vcmap(v)     = cv;
                } else {
                  if (vcmap(v) != ORD_MAX) {
                    vcmap(u) = vcmap(v);
                  } else {
                    match(u) = ORD_MAX;
                  }
                }
              }
            }
          });
      Kokkos::fence();
      // add the ones that failed to be reprocessed next round
      // maybe count these then create next_perm to save memory?
      Kokkos::parallel_for(
          policy_t(0, perm_length), KOKKOS_LAMBDA(ordinal_t i) {
            ordinal_t u = curr_perm(i);
            if (vcmap(u) == ORD_MAX) {
              ordinal_t add_next  = Kokkos::atomic_fetch_add(&next_length(), 1);
              next_perm(add_next) = u;
            }
          });
      Kokkos::fence();
      swap = swap ^ 1;
      Kokkos::deep_copy(perm_length, next_length);
      curr_perm = next_perm;
    }
    ordinal_t nc = 0;
    Kokkos::deep_copy(nc, nvertices_coarse);
    return nc;
  }

  // hn is a list of vertices such that vertex i wants to aggregate with vertex
  // hn(i)
  static ordinal_t parallel_map_construct(vtx_view_t vcmap, const ordinal_t n, const vtx_view_t vperm,
                                          const vtx_view_t hn, const vtx_view_t ordering) {
    vtx_view_t match("match", n);
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) { match(i) = ORD_MAX; });
    ordinal_t perm_length = n;
    Kokkos::View<ordinal_t, Device> nvertices_coarse("nvertices");

    // construct mapping using heaviest edges
    int swap             = 1;
    vtx_view_t curr_perm = vperm;
    while (perm_length > 0) {
      vtx_view_t next_perm("next perm", perm_length);
      Kokkos::View<ordinal_t, Device> next_length("next_length");

      Kokkos::parallel_for(
          policy_t(0, perm_length), KOKKOS_LAMBDA(ordinal_t i) {
            ordinal_t u   = curr_perm(i);
            ordinal_t v   = hn(u);
            int condition = ordering(u) < ordering(v);
            // need to enforce an ordering condition to allow hard-stall
            // conditions to be broken
            if (condition ^ swap) {
              if (Kokkos::atomic_compare_exchange_strong(&match(u), ORD_MAX, v)) {
                if (u == v || Kokkos::atomic_compare_exchange_strong(&match(v), ORD_MAX, u)) {
                  ordinal_t cv = u;
                  if (v < u) {
                    cv = v;
                  }
                  vcmap(u) = cv;
                  vcmap(v) = cv;
                } else {
                  if (vcmap(v) != ORD_MAX) {
                    vcmap(u) = vcmap(v);
                  } else {
                    match(u) = ORD_MAX;
                  }
                }
              }
            }
          });
      Kokkos::fence();
      // add the ones that failed to be reprocessed next round
      // maybe count these then create next_perm to save memory?
      Kokkos::parallel_scan(
          policy_t(0, perm_length), KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
            ordinal_t u = curr_perm(i);
            if (vcmap(u) == ORD_MAX) {
              if (final) {
                next_perm(update) = u;
              }
              update++;
            }
            if (final && (i + 1) == perm_length) {
              next_length() = update;
            }
          });
      Kokkos::fence();
      swap = swap ^ 1;
      Kokkos::deep_copy(perm_length, next_length);
      curr_perm = next_perm;
    }
    Kokkos::parallel_scan(
        "assign aggregates", policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t u, ordinal_t& update, const bool final) {
          if (vcmap(u) == u) {
            if (final) {
              vcmap(u) = update;
            }
            update++;
          } else if (final) {
            vcmap(u) = vcmap(u) + n;
          }
          if (final && (u + 1) == n) {
            nvertices_coarse() = update;
          }
        });
    Kokkos::parallel_for(
        "propagate aggregates", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t u) {
          if (vcmap(u) >= n) {
            ordinal_t c_id = vcmap(u) - n;
            vcmap(u)       = vcmap(c_id);
          }
        });
    ordinal_t nc = 0;
    Kokkos::deep_copy(nc, nvertices_coarse);
    return nc;
  }

  static part_view_t GOSH_clusters(const matrix_t& g) {
    // finds the central vertices for GOSH clusters
    // approximately this is a maximal independent set (if you pretend edges
    // whose endpoints both exceed degree thresholds don't exist) IS vertices
    // are preferred to be vertices with high degree, so it should be small

    ordinal_t n = g.numRows();

    // 0: unassigned
    // 1: in IS
    //-1: adjacent to an IS vertex
    part_view_t state("psuedo is membership", n);

    ordinal_t unassigned_total = n;

    // gonna keep this as an edge view in case we wanna do weighted degree
    edge_view_t degrees("degrees", n);
    vtx_view_t unassigned("unassigned vertices", n);
    Kokkos::parallel_for(
        "populate degrees", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          degrees(i)    = g.graph.row_map(i + 1) - g.graph.row_map(i);
          unassigned(i) = i;
        });
    edge_offset_t threshold = g.nnz() / g.numRows();

    while (unassigned_total > 0) {
      part_view_t tuple_state("tuple state", n);
      edge_view_t tuple_degree("tuple rand", n);
      vtx_view_t tuple_idx("tuple index", n);

      part_view_t tuple_state_update("tuple state", n);
      edge_view_t tuple_degree_update("tuple rand", n);
      vtx_view_t tuple_idx_update("tuple index", n);
      Kokkos::parallel_for(
          policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t i) {
            tuple_state(i)  = state(i);
            tuple_degree(i) = degrees(i);
            tuple_idx(i)    = i;
          });

      Kokkos::parallel_for(
          policy_t(0, unassigned_total), KOKKOS_LAMBDA(const ordinal_t i) {
            ordinal_t u              = unassigned(i);
            int max_state            = tuple_state(u);
            edge_offset_t max_degree = tuple_degree(u);
            ordinal_t max_idx        = tuple_idx(u);

            for (edge_offset_t j = g.graph.row_map(u); j < g.graph.row_map(u + 1); j++) {
              ordinal_t v = g.graph.entries(j);
              bool is_max = false;
              if (tuple_state(v) > max_state) {
                is_max = true;
              } else if (tuple_state(v) == max_state) {
                if (tuple_degree(v) > max_degree) {
                  is_max = true;
                } else if (tuple_degree(v) == max_degree) {
                  if (tuple_idx(v) > max_idx) {
                    is_max = true;
                  }
                }
              }
              // pretend edges between two vertices exceeding threshold do not
              // exist
              if (degrees(u) > threshold && degrees(v) > threshold) {
                is_max = false;
              }
              if (is_max) {
                max_state  = tuple_state(v);
                max_degree = tuple_degree(v);
                max_idx    = tuple_idx(v);
              }
            }
            tuple_state_update(u)  = max_state;
            tuple_degree_update(u) = max_degree;
            tuple_idx_update(u)    = max_idx;
          });

      Kokkos::parallel_for(
          policy_t(0, unassigned_total), KOKKOS_LAMBDA(const ordinal_t i) {
            ordinal_t u     = unassigned(i);
            tuple_state(u)  = tuple_state_update(u);
            tuple_degree(u) = tuple_degree_update(u);
            tuple_idx(u)    = tuple_idx_update(u);
          });

      ordinal_t next_unassigned_total = 0;
      Kokkos::parallel_reduce(
          policy_t(0, unassigned_total),
          KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& thread_sum) {
            ordinal_t u = unassigned(i);
            if (state(u) == 0) {
              if (tuple_idx(u) == u) {
                state(u) = 1;
              }
              // check if at least one of neighbors are in the IS or will be
              // placed into the IS
              else if (tuple_state(u) == 1 || tuple_idx(tuple_idx(u)) == tuple_idx(u)) {
                state(u) = -1;
              }
            }
            if (state(u) == 0) {
              thread_sum++;
            }
          },
          next_unassigned_total);

      vtx_view_t next_unassigned("next unassigned", next_unassigned_total);
      Kokkos::parallel_scan(
          "create next unassigned", policy_t(0, unassigned_total),
          KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
            ordinal_t u = unassigned(i);
            if (state(u) == 0) {
              if (final) {
                next_unassigned(update) = u;
              }
              update++;
            }
          });
      unassigned_total = next_unassigned_total;
      unassigned       = next_unassigned;
    }
    return state;
  }

  static matrix_t coarsen_mis_2(const matrix_t& g) {
    ordinal_t n = g.numRows();

    typename matrix_t::staticcrsgraph_type::entries_type::non_const_value_type nc = 0;
    vtx_view_t vcmap =
        KokkosGraph::graph_mis2_aggregate<Device, typename matrix_t::staticcrsgraph_type::row_map_type,
                                          typename matrix_t::staticcrsgraph_type::entries_type, vtx_view_t>(
            g.graph.row_map, g.graph.entries, nc);

    edge_view_t row_map("interpolate row map", n + 1);

    Kokkos::parallel_for(
        policy_t(0, n + 1), KOKKOS_LAMBDA(ordinal_t u) { row_map(u) = u; });

    vtx_view_t entries("interpolate entries", n);
    wgt_view_t values("interpolate values", n);
    // compute the interpolation weights
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t u) {
          entries(u) = vcmap(u);
          values(u)  = 1.0;
        });

    graph_type graph(entries, row_map);
    matrix_t interp("interpolate", nc, values, graph);

    return interp;
  }

  static matrix_t coarsen_GOSH(const matrix_t& g) {
    ordinal_t n = g.numRows();

    part_view_t colors = GOSH_clusters(g);

    Kokkos::View<ordinal_t, Device> nvc("nvertices_coarse");
    vtx_view_t vcmap("vcmap", n);

    int first_color = 1;

    // create aggregates for color 1
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          if (colors(i) == first_color) {
            vcmap(i) = Kokkos::atomic_fetch_add(&nvc(), 1);
          } else {
            vcmap(i) = ORD_MAX;
          }
        });

    // add unaggregated vertices to aggregate of highest degree neighbor
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          if (colors(i) != first_color) {
            // could use a thread team here
            edge_offset_t max_degree = 0;
            for (edge_offset_t j = g.graph.row_map(i); j < g.graph.row_map(i + 1); j++) {
              ordinal_t v          = g.graph.entries(j);
              edge_offset_t degree = g.graph.row_map(v + 1) - g.graph.row_map(v);
              if (colors(v) == first_color && degree > max_degree) {
                max_degree = degree;
                vcmap(i)   = vcmap(v);
              }
            }
          }
        });

    ordinal_t nc = 0;
    Kokkos::deep_copy(nc, nvc);

    edge_view_t row_map("interpolate row map", n + 1);

    Kokkos::parallel_for(
        policy_t(0, n + 1), KOKKOS_LAMBDA(ordinal_t u) { row_map(u) = u; });

    vtx_view_t entries("interpolate entries", n);
    wgt_view_t values("interpolate values", n);
    // compute the interpolation weights
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t u) {
          entries(u) = vcmap(u);
          values(u)  = 1.0;
        });

    graph_type graph(entries, row_map);
    matrix_t interp("interpolate", nc, values, graph);

    return interp;
  }

  static matrix_t coarsen_GOSH_v2(const matrix_t& g) {
    ordinal_t n = g.numRows();

    Kokkos::View<ordinal_t, Device> nvc("nvertices_coarse");
    vtx_view_t vcmap("vcmap", n);

    edge_offset_t threshold_d = g.nnz() / n;
    if (threshold_d < 50) {
      threshold_d = 50;
    }
    // create aggregates for large degree vtx
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          if (g.graph.row_map(i + 1) - g.graph.row_map(i) > threshold_d) {
            ordinal_t cv = Kokkos::atomic_fetch_add(&nvc(), 1);
            vcmap(i)     = cv;
          } else {
            vcmap(i) = ORD_MAX;
          }
        });

    // add vertex to max wgt neighbor's aggregate
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          if (vcmap(i) == ORD_MAX) {
            ordinal_t argmax = ORD_MAX;
            scalar_t max_w   = 0;
            for (edge_offset_t j = g.graph.row_map(i); j < g.graph.row_map(i + 1); j++) {
              ordinal_t v   = g.graph.entries(j);
              ordinal_t wgt = g.values(j);
              if (vcmap(v) != ORD_MAX) {
                if (wgt >= max_w) {
                  max_w  = wgt;
                  argmax = v;
                }
              }
            }
            if (argmax != ORD_MAX) {
              vcmap(i) = vcmap(argmax);
            }
          }
        });

    // add vertex to max degree neighbor's aggregate
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          if (vcmap(i) == ORD_MAX) {
            ordinal_t argmax    = ORD_MAX;
            edge_offset_t max_d = 0;
            for (edge_offset_t j = g.graph.row_map(i); j < g.graph.row_map(i + 1); j++) {
              ordinal_t v          = g.graph.entries(j);
              edge_offset_t degree = g.graph.row_map(v + 1) - g.graph.row_map(v);
              if (vcmap(v) != ORD_MAX) {
                if (degree >= max_d) {
                  max_d  = degree;
                  argmax = v;
                }
              }
            }
            if (argmax != ORD_MAX) {
              vcmap(i) = vcmap(argmax);
            }
          }
        });

    // add neighbors of each aggregated vertex to aggregate
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          if (vcmap(i) != ORD_MAX) {
            for (edge_offset_t j = g.graph.row_map(i); j < g.graph.row_map(i + 1); j++) {
              ordinal_t v = g.graph.entries(j);
              if (vcmap(v) == ORD_MAX) {
                vcmap(v) = vcmap(i);
              }
            }
          }
        });

    ordinal_t remaining_total = 0;

    Kokkos::parallel_reduce(
        "count remaining", policy_t(0, n),
        KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& sum) {
          if (vcmap(i) == ORD_MAX) {
            sum++;
          }
        },
        remaining_total);

    vtx_view_t remaining("remaining vtx", remaining_total);

    Kokkos::parallel_scan(
        "count remaining", policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
          if (vcmap(i) == ORD_MAX) {
            if (final) {
              remaining(update) = i;
            }
            update++;
          }
        });

    vtx_view_t hn("heaviest neighbors", n);

    pool_t rand_pool(std::time(nullptr));

    Kokkos::parallel_for(
        "fill hn", policy_t(0, remaining_total), KOKKOS_LAMBDA(ordinal_t r_idx) {
          // select heaviest neighbor with ties randomly broken
          ordinal_t i       = remaining(r_idx);
          ordinal_t hn_i    = ORD_MAX;
          uint64_t max_rand = 0;
          scalar_t max_ewt  = 0;

          edge_offset_t end_offset = g.graph.row_map(i + 1);
          for (edge_offset_t j = g.graph.row_map(i); j < end_offset; j++) {
            scalar_t wgt    = g.values(j);
            ordinal_t v     = g.graph.entries(j);
            gen_t generator = rand_pool.get_state();
            uint64_t rand   = generator.urand64();
            rand_pool.free_state(generator);
            bool choose = false;
            if (max_ewt < wgt) {
              choose = true;
            } else if (max_ewt == wgt && max_rand <= rand) {
              choose = true;
            }

            if (choose) {
              max_ewt  = wgt;
              max_rand = rand;
              hn_i     = v;
            }
          }
          hn(i) = hn_i;
        });

    ordinal_t nc = parallel_map_construct_prefilled(vcmap, n, remaining, hn, nvc);
    Kokkos::deep_copy(nc, nvc);

    edge_view_t row_map("interpolate row map", n + 1);

    Kokkos::parallel_for(
        policy_t(0, n + 1), KOKKOS_LAMBDA(ordinal_t u) { row_map(u) = u; });

    vtx_view_t entries("interpolate entries", n);
    wgt_view_t values("interpolate values", n);
    // compute the interpolation weights
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t u) {
          entries(u) = vcmap(u);
          values(u)  = 1.0;
        });

    graph_type graph(entries, row_map);
    matrix_t interp("interpolate", nc, values, graph);

    return interp;
  }

  static matrix_t coarsen_HEC(const matrix_t& g, bool uniform_weights) {
    ordinal_t n = g.numRows();

    vtx_view_t hn("heavies", n);

    vtx_view_t vcmap("vcmap", n);

    Kokkos::parallel_for(
        "initialize vcmap", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) { vcmap(i) = ORD_MAX; });

    pool_t rand_pool(std::time(nullptr));

    vtx_view_t vperm = generate_permutation(n, rand_pool);

    vtx_view_t reverse_map("reversed", n);
    Kokkos::parallel_for(
        "construct reverse map", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) { reverse_map(vperm(i)) = i; });

    if (uniform_weights) {
      // all weights equal at this level so choose heaviest edge randomly
      Kokkos::parallel_for(
          "Random HN", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
            gen_t generator    = rand_pool.get_state();
            ordinal_t adj_size = g.graph.row_map(i + 1) - g.graph.row_map(i);
            if (adj_size > 0) {
              ordinal_t offset = g.graph.row_map(i) + (generator.urand64() % adj_size);
              hn(i)            = g.graph.entries(offset);
            } else {
              hn(i) = generator.urand64() % n;
            }
            rand_pool.free_state(generator);
          });
    } else {
      Kokkos::parallel_for(
          "Heaviest HN", team_policy_t(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member& thread) {
            ordinal_t i        = thread.league_rank();
            ordinal_t adj_size = g.graph.row_map(i + 1) - g.graph.row_map(i);
            if (adj_size > 0) {
              edge_offset_t end = g.graph.row_map(i + 1);
              typename Kokkos::MaxLoc<scalar_t, edge_offset_t, Device>::value_type argmax{};
              Kokkos::parallel_reduce(
                  Kokkos::TeamThreadRange(thread, g.graph.row_map(i), end),
                  [=](const edge_offset_t idx, Kokkos::ValLocScalar<scalar_t, edge_offset_t>& local) {
                    scalar_t wgt = g.values(idx);
                    if (wgt >= local.val) {
                      local.val = wgt;
                      local.loc = idx;
                    }
                  },
                  Kokkos::MaxLoc<scalar_t, edge_offset_t, Device>(argmax));
              Kokkos::single(Kokkos::PerTeam(thread), [=]() {
                ordinal_t h = g.graph.entries(argmax.loc);
                hn(i)       = h;
              });
            } else {
              gen_t generator = rand_pool.get_state();
              hn(i)           = generator.urand64() % n;
              rand_pool.free_state(generator);
            }
          });
    }
    ordinal_t nc = 0;
    nc           = parallel_map_construct(vcmap, n, vperm, hn, reverse_map);

    edge_view_t row_map("interpolate row map", n + 1);

    Kokkos::parallel_for(
        policy_t(0, n + 1), KOKKOS_LAMBDA(ordinal_t u) { row_map(u) = u; });

    vtx_view_t entries("interpolate entries", n);
    wgt_view_t values("interpolate values", n);
    // compute the interpolation weights
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t u) {
          entries(u) = vcmap(u);
          values(u)  = 1.0;
        });

    graph_type graph(entries, row_map);
    matrix_t interp("interpolate", nc, values, graph);

    return interp;
  }

  static ordinal_t countInf(vtx_view_t target) {
    ordinal_t totalInf = 0;

    Kokkos::parallel_reduce(
        policy_t(0, target.extent(0)),
        KOKKOS_LAMBDA(ordinal_t i, ordinal_t & thread_sum) {
          if (target(i) == ORD_MAX) {
            thread_sum++;
          }
        },
        totalInf);

    return totalInf;
  }

  struct MatchByHashSorted {
    vtx_view_t vcmap, unmapped;
    Kokkos::View<uint32_t*, Device> hashes;
    ordinal_t unmapped_total;
    Kokkos::View<ordinal_t, Device> nvertices_coarse;
    MatchByHashSorted(vtx_view_t _vcmap, vtx_view_t _unmapped, Kokkos::View<uint32_t*, Device> _hashes,
                      ordinal_t _unmapped_total, Kokkos::View<ordinal_t, Device> _nvertices_coarse)
        : vcmap(_vcmap),
          unmapped(_unmapped),
          hashes(_hashes),
          unmapped_total(_unmapped_total),
          nvertices_coarse(_nvertices_coarse) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const ordinal_t i, ordinal_t& update, const bool final) const {
      ordinal_t u         = unmapped(i);
      ordinal_t tentative = 0;
      if (i == 0) {
        tentative = i;
      } else if (hashes(i - 1) != hashes(i)) {
        tentative = i;
      }

      if (tentative > update) {
        update = tentative;
      }

      if (final) {
        // update should contain the index of the first hash that equals
        // hash(i), could be i we want to determine if i is an odd offset from
        // update
        ordinal_t isOddOffset = (i - update) & 1;
        // if even (0 counts as even) we match unmapped(i) with unmapped(i+1) if
        // hash(i) == hash(i+1) if odd do nothing
        if (isOddOffset == 0) {
          if (i + 1 < unmapped_total) {
            if (hashes(i) == hashes(i + 1)) {
              ordinal_t v = unmapped(i + 1);
              vcmap(u)    = Kokkos::atomic_fetch_add(&nvertices_coarse(), 1);
              vcmap(v)    = vcmap(u);
            }
          }
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void join(ordinal_t& update, const ordinal_t& input) const {
      if (input > update) update = input;
    }
  };

  static matrix_t coarsen_match(const matrix_t& g, bool uniform_weights, int match_choice) {
    ordinal_t n = g.numRows();

    vtx_view_t hn("heavies", n);

    vtx_view_t vcmap("vcmap", n);

    Kokkos::parallel_for(
        "initialize vcmap", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) { vcmap(i) = ORD_MAX; });

    rand_view_t randoms("randoms", n);

    pool_t rand_pool(std::time(nullptr));

    vtx_view_t vperm = generate_permutation(n, rand_pool);

    vtx_view_t reverse_map("reversed", n);
    Kokkos::parallel_for(
        "construct reverse map", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) { reverse_map(vperm(i)) = i; });

    if (uniform_weights) {
      // all weights equal at this level so choose heaviest edge randomly
      Kokkos::parallel_for(
          "Random HN", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
            gen_t generator    = rand_pool.get_state();
            ordinal_t adj_size = g.graph.row_map(i + 1) - g.graph.row_map(i);
            ordinal_t offset   = g.graph.row_map(i) + (generator.urand64() % adj_size);
            hn(i)              = g.graph.entries(offset);
            rand_pool.free_state(generator);
          });
    } else {
      Kokkos::parallel_for(
          "Heaviest HN", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
            ordinal_t hn_i   = g.graph.entries(g.graph.row_map(i));
            scalar_t max_ewt = g.values(g.graph.row_map(i));

            edge_offset_t end_offset = g.graph.row_map(i + 1);  // +g.edges_per_source[i];

            for (edge_offset_t j = g.graph.row_map(i) + 1; j < end_offset; j++) {
              if (max_ewt < g.values(j)) {
                max_ewt = g.values(j);
                hn_i    = g.graph.entries(j);
              }
            }
            hn(i) = hn_i;
          });
    }
    vtx_view_t match("match", n);
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) { match(i) = ORD_MAX; });

    ordinal_t perm_length = n;

    Kokkos::View<ordinal_t, Device> nvertices_coarse("nvertices");

    // construct mapping using heaviest edges
    int swap = 1;
    while (perm_length > 0) {
      vtx_view_t next_perm("next perm", perm_length);
      Kokkos::View<ordinal_t, Device> next_length("next_length");

      // match vertices with heaviest unmatched edge
      Kokkos::parallel_for(
          policy_t(0, perm_length), KOKKOS_LAMBDA(ordinal_t i) {
            ordinal_t u   = vperm(i);
            ordinal_t v   = hn(u);
            int condition = reverse_map(u) < reverse_map(v);
            // need to enforce an ordering condition to allow hard-stall
            // conditions to be broken
            if (condition ^ swap) {
              if (Kokkos::atomic_compare_exchange_strong(&match(u), ORD_MAX, v)) {
                if (u == v || Kokkos::atomic_compare_exchange_strong(&match(v), ORD_MAX, u)) {
                  // u == v avoids problems if there is a self-loop edge
                  ordinal_t cv = Kokkos::atomic_fetch_add(&nvertices_coarse(), 1);
                  vcmap(u)     = cv;
                  vcmap(v)     = cv;
                } else {
                  match(u) = ORD_MAX;
                }
              }
            }
          });
      Kokkos::fence();

      // add the ones that failed to be reprocessed next round
      // maybe count these then create next_perm to save memory?
      Kokkos::parallel_for(
          policy_t(0, perm_length), KOKKOS_LAMBDA(ordinal_t i) {
            ordinal_t u = vperm(i);
            if (vcmap(u) == ORD_MAX) {
              ordinal_t h = ORD_MAX;

              if (uniform_weights) {
                ordinal_t max_ewt = 0;
                // we have to iterate over the edges anyways because we need to
                // check if any are unmatched! so instead of randomly choosing a
                // heaviest edge, we instead use the reverse permutation order
                // as the weight
                for (edge_offset_t j = g.graph.row_map(u); j < g.graph.row_map(u + 1); j++) {
                  ordinal_t v = g.graph.entries(j);
                  // v must be unmatched to be considered
                  if (vcmap(v) == ORD_MAX) {
                    // using <= so that zero weight edges may still be chosen
                    if (max_ewt <= reverse_map(v)) {
                      max_ewt = reverse_map(v);
                      h       = v;
                    }
                  }
                }
              } else {
                scalar_t max_ewt = 0;
                for (edge_offset_t j = g.graph.row_map(u); j < g.graph.row_map(u + 1); j++) {
                  ordinal_t v = g.graph.entries(j);
                  // v must be unmatched to be considered
                  if (vcmap(v) == ORD_MAX) {
                    // using <= so that zero weight edges may still be chosen
                    if (max_ewt <= g.values(j)) {
                      max_ewt = g.values(j);
                      h       = v;
                    }
                  }
                }
              }

              if (h != ORD_MAX) {
                ordinal_t add_next  = Kokkos::atomic_fetch_add(&next_length(), 1);
                next_perm(add_next) = u;
                hn(u)               = h;
              }
            }
          });
      Kokkos::fence();
      swap = swap ^ 1;
      Kokkos::deep_copy(perm_length, next_length);
      vperm = next_perm;
    }

    if (match_choice == 1) {
      ordinal_t unmapped   = countInf(vcmap);
      double unmappedRatio = static_cast<double>(unmapped) / static_cast<double>(n);

      // leaf matches
      if (unmappedRatio > 0.25) {
        Kokkos::parallel_for(
            policy_t(0, n), KOKKOS_LAMBDA(ordinal_t u) {
              if (vcmap(u) != ORD_MAX) {
                ordinal_t lastLeaf = ORD_MAX;
                for (edge_offset_t j = g.graph.row_map(u); j < g.graph.row_map(u + 1); j++) {
                  ordinal_t v = g.graph.entries(j);
                  // v must be unmatched to be considered
                  if (vcmap(v) == ORD_MAX) {
                    // must be degree 1 to be a leaf
                    if (g.graph.row_map(v + 1) - g.graph.row_map(v) == 1) {
                      if (lastLeaf == ORD_MAX) {
                        lastLeaf = v;
                      } else {
                        vcmap(lastLeaf) = Kokkos::atomic_fetch_add(&nvertices_coarse(), 1);
                        vcmap(v)        = vcmap(lastLeaf);
                        lastLeaf        = ORD_MAX;
                      }
                    }
                  }
                }
              }
            });
      }

      unmapped      = countInf(vcmap);
      unmappedRatio = static_cast<double>(unmapped) / static_cast<double>(n);

      // twin matches
      if (unmappedRatio > 0.25) {
        vtx_view_t unmappedVtx("unmapped vertices", unmapped);
        Kokkos::View<uint32_t*, Device> hashes("hashes", unmapped);

        Kokkos::View<ordinal_t, Device> unmappedIdx("unmapped index");
        hasher_t hasher;
        // compute digests of adjacency lists
        Kokkos::parallel_for(
            "create digests", team_policy_t(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member& thread) {
              ordinal_t u = thread.league_rank();
              if (vcmap(u) == ORD_MAX) {
                uint32_t hash = 0;
                Kokkos::parallel_reduce(
                    Kokkos::TeamThreadRange(thread, g.graph.row_map(u), g.graph.row_map(u + 1)),
                    [=](const edge_offset_t j, uint32_t& thread_sum) { thread_sum += hasher(g.graph.entries(j)); },
                    hash);
                Kokkos::single(Kokkos::PerTeam(thread), [=]() {
                  ordinal_t idx    = Kokkos::atomic_fetch_add(&unmappedIdx(), 1);
                  unmappedVtx(idx) = u;
                  hashes(idx)      = hash;
                });
              }
            });
        uint32_t max = std::numeric_limits<uint32_t>::max();
        typedef Kokkos::BinOp1D<Kokkos::View<uint32_t*, Device> > BinOp;
        BinOp bin_op(unmapped, 0, max);
        // VERY important that final parameter is true
        Kokkos::BinSort<Kokkos::View<uint32_t*, Device>, BinOp, exec_space, ordinal_t> sorter(hashes, bin_op, true);
        sorter.create_permute_vector();
        sorter.template sort<Kokkos::View<uint32_t*, Device> >(hashes);
        sorter.template sort<vtx_view_t>(unmappedVtx);

        MatchByHashSorted matchTwinFunctor(vcmap, unmappedVtx, hashes, unmapped, nvertices_coarse);
        Kokkos::parallel_scan("match twins", policy_t(0, unmapped), matchTwinFunctor);
      }

      unmapped      = countInf(vcmap);
      unmappedRatio = static_cast<double>(unmapped) / static_cast<double>(n);

      // relative matches
      if (unmappedRatio > 0.25) {
        // get possibly mappable vertices of unmapped
        vtx_view_t mappableVtx("mappable vertices", unmapped);
        Kokkos::parallel_scan(
            "get unmapped", policy_t(0, n), KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
              if (vcmap(i) == ORD_MAX) {
                if (final) {
                  mappableVtx(update) = i;
                }

                update++;
              }
            });

        ordinal_t mappable_count = unmapped;
        do {
          Kokkos::parallel_for(
              "reset hn", policy_t(0, mappable_count), KOKKOS_LAMBDA(ordinal_t i) {
                ordinal_t u = mappableVtx(i);
                hn(u)       = ORD_MAX;
              });

          // choose relatives for unmapped vertices
          Kokkos::parallel_for(
              "assign relatives", policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
                if (vcmap(i) != ORD_MAX) {
                  ordinal_t last_free = ORD_MAX;
                  for (edge_offset_t j = g.graph.row_map(i); j < g.graph.row_map(i + 1); j++) {
                    ordinal_t v = g.graph.entries(j);
                    if (vcmap(v) == ORD_MAX) {
                      if (last_free != ORD_MAX) {
                        // there can be multiple threads updating this but it
                        // doesn't matter as long as they have some value
                        hn(last_free) = v;
                        hn(v)         = last_free;
                        last_free     = ORD_MAX;
                      } else {
                        last_free = v;
                      }
                    }
                  }
                }
              });

          // create a list of all mappable vertices according to set entries of
          // hn
          ordinal_t old_mappable = mappable_count;
          mappable_count         = 0;
          Kokkos::parallel_reduce(
              "count mappable", policy_t(0, old_mappable),
              KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& thread_sum) {
                ordinal_t u = mappableVtx(i);
                if (hn(u) != ORD_MAX) {
                  thread_sum++;
                }
              },
              mappable_count);

          vtx_view_t nextMappable("next mappable vertices", mappable_count);

          Kokkos::parallel_scan(
              "get next mappable", policy_t(0, old_mappable),
              KOKKOS_LAMBDA(const ordinal_t i, ordinal_t& update, const bool final) {
                ordinal_t u = mappableVtx(i);
                if (hn(u) != ORD_MAX) {
                  if (final) {
                    nextMappable(update) = u;
                  }

                  update++;
                }
              });
          mappableVtx = nextMappable;

          // match vertices with chosen relative
          if (mappable_count > 0) {
            Kokkos::parallel_for(
                policy_t(0, mappable_count), KOKKOS_LAMBDA(ordinal_t i) {
                  ordinal_t u   = mappableVtx(i);
                  ordinal_t v   = hn(u);
                  int condition = reverse_map(u) < reverse_map(v);
                  // need to enforce an ordering condition to allow hard-stall
                  // conditions to be broken
                  if (condition ^ swap) {
                    if (Kokkos::atomic_compare_exchange_strong(&match(u), ORD_MAX, v)) {
                      if (Kokkos::atomic_compare_exchange_strong(&match(v), ORD_MAX, u)) {
                        ordinal_t cv = Kokkos::atomic_fetch_add(&nvertices_coarse(), 1);
                        vcmap(u)     = cv;
                        vcmap(v)     = cv;
                      } else {
                        match(u) = ORD_MAX;
                      }
                    }
                  }
                });
          }
          Kokkos::fence();
          swap = swap ^ 1;
        } while (mappable_count > 0);
      }
    }

    // create singleton aggregates of remaining unmatched vertices
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t i) {
          if (vcmap(i) == ORD_MAX) {
            vcmap(i) = Kokkos::atomic_fetch_add(&nvertices_coarse(), 1);
          }
        });

    ordinal_t nc = 0;
    Kokkos::deep_copy(nc, nvertices_coarse);

    edge_view_t row_map("interpolate row map", n + 1);

    Kokkos::parallel_for(
        policy_t(0, n + 1), KOKKOS_LAMBDA(ordinal_t u) { row_map(u) = u; });

    vtx_view_t entries("interpolate entries", n);
    wgt_view_t values("interpolate values", n);
    // compute the interpolation weights
    Kokkos::parallel_for(
        policy_t(0, n), KOKKOS_LAMBDA(ordinal_t u) {
          entries(u) = vcmap(u);
          values(u)  = 1.0;
        });

    graph_type graph(entries, row_map);
    matrix_t interp("interpolate", nc, values, graph);

    return interp;
  }
};

}  // end namespace Experimental
}  // end namespace KokkosGraph
// exclude from Cuda builds without lambdas enabled
#endif
