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

#if defined(KOKKOS_ENABLE_OPENMP)
#ifdef KOKKOSKERNELS_HAVE_OUTER
#include <parallel/multiseq_selection.h>
#include <parallel/multiway_merge.h>
#include <parallel/merge.h>
#include <parallel/multiway_mergesort.h>
#endif
#endif

namespace KokkosSparse {

namespace Impl {

#if defined(KOKKOS_ENABLE_OPENMP)
#ifdef KOKKOSKERNELS_HAVE_OUTER
template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::Triplet {
  nnz_lno_t src, dst;  //, block;
  scalar_t val;
  bool operator<(const Triplet &a) const { return (this->src < a.src) || (this->src == a.src && this->dst < a.dst); }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_col_view_t, typename a_nnz_view_t, typename a_scalar_view_t, typename b_row_view_t,
          typename b_nnz_view_t, typename b_scalar_view_t, typename flop_row_view_t>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::OuterProduct {
  nnz_lno_t begin, end, team_work_size;
  a_col_view_t a_col_xadj;
  a_nnz_view_t a_adj;
  a_scalar_view_t a_vals;

  b_row_view_t b_row_xadj;
  b_nnz_view_t b_adj;
  b_scalar_view_t b_vals;

  flop_row_view_t flop_per_row;

  Kokkos::View<Triplet *, MyTempMemorySpace> triplets;

  OuterProduct(nnz_lno_t begin_, nnz_lno_t end_, nnz_lno_t team_work_size_, a_col_view_t a_col_xadj_,
               a_nnz_view_t a_adj_, a_scalar_view_t a_vals_, b_row_view_t b_row_xadj_, b_nnz_view_t b_adj_,
               b_scalar_view_t b_vals_, flop_row_view_t flop_per_row_,
               Kokkos::View<Triplet *, MyTempMemorySpace> triplets_)
      : begin(begin_),
        end(end_),
        team_work_size(team_work_size_),
        a_col_xadj(a_col_xadj_),
        a_adj(a_adj_),
        a_vals(a_vals_),
        b_row_xadj(b_row_xadj_),
        b_adj(b_adj_),
        b_vals(b_vals_),
        flop_per_row(flop_per_row_),
        triplets(triplets_) {}

  // assumes that the vector lane is 1, as in cpus
  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member_t &teamMember) const {
    nnz_lno_t team_row_begin     = teamMember.league_rank() * team_work_size + begin;
    const nnz_lno_t team_row_end = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, end);

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t &row_col_ind) {
          const size_type a_col_begin = a_col_xadj[row_col_ind];
          const size_type a_col_end   = a_col_xadj[row_col_ind + 1];
          size_t write_begin_index    = flop_per_row[row_col_ind] - flop_per_row[begin];

          const nnz_lno_t a_col_size = a_col_end - a_col_begin;

          const size_type b_row_begin = b_row_xadj[row_col_ind];
          const size_type b_row_end   = b_row_xadj[row_col_ind + 1];
          const nnz_lno_t b_row_size  = b_row_end - b_row_begin;

          for (nnz_lno_t i = 0; i < a_col_size; ++i) {
            size_type a_nnz_ind = i + a_col_begin;

            nnz_lno_t a_row = a_adj[a_nnz_ind];
            scalar_t a_val  = a_vals[a_nnz_ind];

            Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, b_row_size), [&](const nnz_lno_t b_ind) {
              size_type b_nnz_ind       = b_ind + b_row_begin;
              nnz_lno_t b_col           = b_adj[b_nnz_ind];
              scalar_t b_val            = b_vals[b_nnz_ind];
              scalar_t c_val            = b_val * a_val;
              size_t write_index        = write_begin_index + b_ind;
              triplets[write_index].src = a_row;
              triplets[write_index].dst = b_col;
              triplets[write_index].val = c_val;
            });
            write_begin_index += b_row_size;
          }
        });
  }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename triplet_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::sort_triplets(triplet_view_t triplets,
                                                                          size_t num_triplets) {
  // std::sort (triplets.data(), triplets.data() + num_triplets);
  typedef typename triplet_view_t::value_type element_t;
  __gnu_parallel::parallel_sort_mwms<false, true, element_t *>(triplets.data(), triplets.data() + num_triplets,
                                                               std::less<element_t>(), this->concurrency);
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename host_triplet_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_,
                  b_scalar_nnz_view_t_>::merge_triplets_on_slow_memory(host_triplet_view_t *triplets, size_t num_blocks,
                                                                       size_t overall_size,
                                                                       host_triplet_view_t output_triplets) {
  typedef typename host_triplet_view_t::value_type element_t;
  typedef typename std::pair<element_t *, element_t *> _RAIterTriple;
  std::vector<_RAIterTriple> seqs(num_blocks);
  typedef typename std::vector<_RAIterTriple>::iterator _RAIterTripleIterator;

  for (size_t i = 0; i < num_blocks; ++i) {
    seqs[i].first  = triplets[i].data();
    seqs[i].second = triplets[i].data() + triplets[i].extent(0);
  }

  __gnu_parallel::multiway_merge
      //<_RAIterTripleIterator,
      // element_t*,
      // uint64_t,
      // std::less<element_t> >
      (seqs.begin(), seqs.end(), output_triplets.data(), overall_size, std::less<element_t>());
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename triplet_view_t>
size_t KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::collapse_triplets(triplet_view_t triplets,
                                                                                size_t num_triplets) {
  if (0) return this->collapse_triplets_omp(triplets, num_triplets);
  // if(0)
  else {
    nnz_lno_t src = 0, dst = 0;
    scalar_t val = 0;
    if (num_triplets > 0) {
      src = triplets[0].src;
      dst = triplets[0].dst;
      val = triplets[0].val;
    }

    size_t write_index = 0;
    for (size_t i = 1; i < num_triplets; ++i) {
      if (src == triplets[i].src && dst == triplets[i].dst) {
        val += triplets[i].val;
      } else {
        triplets[write_index].src   = src;
        triplets[write_index].dst   = dst;
        triplets[write_index++].val = val;
        src                         = triplets[i].src;
        dst                         = triplets[i].dst;
        val                         = triplets[i].val;
      }
    }

    triplets[write_index].src   = src;
    triplets[write_index].dst   = dst;
    triplets[write_index++].val = val;
    std::cout << "write_index:" << write_index << std::endl;
    return write_index;
  }
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename triplet_view_t, typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
size_t
KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
             b_scalar_nnz_view_t_>::final_collapse_triplets_omp(triplet_view_t triplets, size_t num_triplets,
                                                                c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_,
                                                                c_scalar_nnz_view_t &valuesC_) {
  int tnum = 1;
#pragma omp parallel
  { tnum = omp_get_num_threads(); }
  std::vector<size_t> num_triplets_begin_index(tnum + 1, 0);
  size_t chunksize = num_triplets / tnum + 1;

  for (int i = 1; i < tnum; ++i) {
    size_t begin = chunksize * i;

    nnz_lno_t src = 0, dst = 0;

    src = triplets[begin].src;
    dst = triplets[begin].dst;
    while (begin > 0) {
      nnz_lno_t src2 = 0, dst2 = 0;

      src2 = triplets[begin - 1].src;
      dst2 = triplets[begin - 1].dst;

      if (src == src2 && dst2 == dst) {
        --begin;
      } else {
        break;
      }
    }
    num_triplets_begin_index[i] = begin;
  }
  num_triplets_begin_index[tnum] = num_triplets;

  std::vector<size_t> num_collapsed_triplets_per_thread(tnum + 1);
  std::vector<size_t> write_triplets_pos_per_thread(tnum + 1);

#pragma omp parallel
  {
    int tid = 0;
    tid     = omp_get_thread_num();

    size_t begin               = num_triplets_begin_index[tid];
    size_t end                 = num_triplets_begin_index[tid + 1];
    size_t triplet_write_index = begin;

    nnz_lno_t src = 0, dst = 0;
    scalar_t val = 0;
    if (begin < end) {
      src = triplets[begin].src;
      dst = triplets[begin].dst;
      val = triplets[begin].val;
    }
    for (size_t i = begin + 1; i < end; ++i) {
      if (src == triplets[i].src && dst == triplets[i].dst) {
        val += triplets[i].val;
      } else {
        triplets[triplet_write_index].src   = src;
        triplets[triplet_write_index].dst   = dst;
        triplets[triplet_write_index++].val = val;

        src = triplets[i].src;
        dst = triplets[i].dst;
        val = triplets[i].val;
      }
    }
    if (begin < end) {
      triplets[triplet_write_index].src   = src;
      triplets[triplet_write_index].dst   = dst;
      triplets[triplet_write_index++].val = val;
    }
    num_collapsed_triplets_per_thread[tid] = triplet_write_index;
  }

  size_t overall_size = 0;
  for (int i = 0; i < tnum; ++i) {
    write_triplets_pos_per_thread[i] = overall_size;
    overall_size += num_collapsed_triplets_per_thread[i] - num_triplets_begin_index[i];
  }
  write_triplets_pos_per_thread[tnum] = overall_size;
  // std::cout << "overall_size:" << overall_size << std::endl;

  // size_t write_index = num_collapsed_triplets_per_thread[0];

  entriesC_   = c_lno_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "entriesC_"), overall_size);
  valuesC_    = c_scalar_nnz_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "valuesC_"), overall_size);
  rowmapC_(0) = 0;
  rowmapC_(this->a_row_cnt) = overall_size;
#pragma omp parallel
  {
    int tid                     = 0;
    tid                         = omp_get_thread_num();
    size_t begin                = num_triplets_begin_index[tid];
    size_t end                  = num_collapsed_triplets_per_thread[tid];
    size_t write_index          = write_triplets_pos_per_thread[tid];
    nnz_lno_t current_row_index = triplets[begin].src;
    rowmapC_(current_row_index) = write_index;

    for (size_t j = begin; j < end; ++j) {
      while (triplets[j].src != current_row_index) {
        rowmapC_(++current_row_index) = write_index + j - begin;
      }
      entriesC_[write_index + j - begin] = triplets[j].dst;
      valuesC_[write_index + j - begin]  = triplets[j].val;
    }

    // std::cout << "current_row_index:" << current_row_index << " rowend:" <<
    // write_index + end - begin << std::endl; write_index += end - begin;
  }

  return overall_size;
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename triplet_view_t>
size_t KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::collapse_triplets_omp(triplet_view_t triplets,
                                                                                    size_t num_triplets,
                                                                                    triplet_view_t out_triplets) {
  int tnum = 1;
#pragma omp parallel
  { tnum = omp_get_num_threads(); }
  std::vector<size_t> num_triplets_begin_index(tnum + 1, 0);
  size_t chunksize = num_triplets / tnum + 1;

  for (int i = 1; i < tnum; ++i) {
    size_t begin = chunksize * i;

    nnz_lno_t src = 0, dst = 0;

    src = triplets[begin].src;
    dst = triplets[begin].dst;
    while (begin > 0) {
      nnz_lno_t src2 = 0, dst2 = 0;

      src2 = triplets[begin - 1].src;
      dst2 = triplets[begin - 1].dst;

      if (src == src2 && dst2 == dst) {
        --begin;
      } else {
        break;
      }
    }
    num_triplets_begin_index[i] = begin;
  }
  num_triplets_begin_index[tnum] = num_triplets;

  std::vector<size_t> num_collapsed_triplets_per_thread(tnum + 1);
  std::vector<size_t> write_triplets_pos_per_thread(tnum + 1);

#pragma omp parallel
  {
    int tid = 0;
    tid     = omp_get_thread_num();

    size_t begin               = num_triplets_begin_index[tid];
    size_t end                 = num_triplets_begin_index[tid + 1];
    size_t triplet_write_index = begin;

    nnz_lno_t src = 0, dst = 0;
    scalar_t val = 0;
    if (begin < end) {
      src = triplets[begin].src;
      dst = triplets[begin].dst;
      val = triplets[begin].val;
    }
    for (size_t i = begin + 1; i < end; ++i) {
      if (src == triplets[i].src && dst == triplets[i].dst) {
        val += triplets[i].val;
      } else {
        triplets[triplet_write_index].src   = src;
        triplets[triplet_write_index].dst   = dst;
        triplets[triplet_write_index++].val = val;

        src = triplets[i].src;
        dst = triplets[i].dst;
        val = triplets[i].val;
      }
    }
    if (begin < end) {
      triplets[triplet_write_index].src   = src;
      triplets[triplet_write_index].dst   = dst;
      triplets[triplet_write_index++].val = val;
    }
    num_collapsed_triplets_per_thread[tid] = triplet_write_index;
  }

  size_t overall_size = 0;
  for (int i = 0; i < tnum; ++i) {
    write_triplets_pos_per_thread[i] = overall_size;
    overall_size += num_collapsed_triplets_per_thread[i] - num_triplets_begin_index[i];
  }
  write_triplets_pos_per_thread[tnum] = overall_size;
  // std::cout << "overall_size:" << overall_size << std::endl;

  // size_t write_index = num_collapsed_triplets_per_thread[0];

#pragma omp parallel
  {
    int tid            = 0;
    tid                = omp_get_thread_num();
    size_t begin       = num_triplets_begin_index[tid];
    size_t end         = num_collapsed_triplets_per_thread[tid];
    size_t write_index = write_triplets_pos_per_thread[tid];

    for (size_t j = begin; j < end; ++j) {
      out_triplets[write_index + j - begin].src = triplets[j].src;
      out_triplets[write_index + j - begin].dst = triplets[j].dst;
      out_triplets[write_index + j - begin].val = triplets[j].val;
    }
    // write_index += end - begin;
  }
  return overall_size;
}

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename a_row_view_t, typename b_row_view_t, typename flop_row_view_t>
struct KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                    b_lno_nnz_view_t_, b_scalar_nnz_view_t_>::FlopsPerRowOuter {
  nnz_lno_t n;            // num rows
  a_row_view_t row_mapA;  // row pointers of a
  b_row_view_t row_mapB;
  flop_row_view_t flop_per_row;
  FlopsPerRowOuter(nnz_lno_t n_, a_row_view_t row_mapA_, b_row_view_t row_mapB_, flop_row_view_t flop_per_row_)
      : n(n_), row_mapA(row_mapA_), row_mapB(row_mapB_), flop_per_row(flop_per_row_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t &i) const {
    flop_per_row(i) = (row_mapB(i + 1) - row_mapB(i)) * (row_mapA(i + 1) - row_mapA(i));
  }

  KOKKOS_INLINE_FUNCTION
  void join(size_type &dst, const size_type &src) const {
    if (dst < src) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init(size_type &dst) const { dst = min_val; }
};

template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void KokkosSPGEMM<HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_,
                  b_lno_nnz_view_t_,
                  b_scalar_nnz_view_t_>::KokkosSPGEMM_numeric_outer(c_row_view_t &rowmapC_, c_lno_nnz_view_t &entriesC_,
                                                                    c_scalar_nnz_view_t &valuesC_,
                                                                    KokkosKernels::Impl::ExecSpaceType my_exec_space) {
  // const size_t block_size = 300000000;
  size_t block_size = 200000000;
  char *env_p;
  if (env_p = std::getenv("BS")) {
    std::cout << "resetting blocksize:" << block_size << " to:" << env_p << std::endl;
    block_size = atoi(env_p);
  }

  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace> mySlowMemory;
  typedef Kokkos::Device<Kokkos::OpenMP, Kokkos::OpenMP::memory_space> myFastMemory;

  ////////////////////
  // TRANSPOSE MATRIX
  ////////////////////
  // get the suggested vectorlane size based on the execution space, and average
  // number of nnzs per row.
  int suggested_vector_size = this->handle->get_suggested_vector_size(b_row_cnt, entriesB.extent(0));
  // get the suggested team size.
  int suggested_team_size = this->handle->get_suggested_team_size(suggested_vector_size);
  // get the chunk size suggested by the handle.
  nnz_lno_t team_row_chunk_size = this->handle->get_team_work_size(suggested_team_size, this->concurrency, b_row_cnt);

  // step-1 tranpose the first matrix.
  Kokkos::Timer timer1, timer_all;
  row_lno_temp_work_view_t transpose_col_xadj("transpose_col_xadj", b_row_cnt + 1);
  nnz_lno_temp_work_view_t transpose_col_adj(Kokkos::view_alloc(Kokkos::WithoutInitializing, "transpose_col_adj"),
                                             entriesA.extent(0));
  scalar_temp_work_view_t tranpose_vals(Kokkos::view_alloc(Kokkos::WithoutInitializing, "transpose_col_values"),
                                        entriesA.extent(0));

  KokkosSparse::Impl::transpose_matrix<const_a_lno_row_view_t, const_a_lno_nnz_view_t, const_a_scalar_nnz_view_t,
                                       row_lno_temp_work_view_t, nnz_lno_temp_work_view_t, scalar_temp_work_view_t,
                                       row_lno_temp_work_view_t, MyExecSpace>(
      a_row_cnt, b_row_cnt, row_mapA, entriesA, valsA, transpose_col_xadj, transpose_col_adj, tranpose_vals);

  MyExecSpace().fence();
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << "\t\tTranspose FlopsPerRowOuterCal BlockPartition FastAllocation";
    std::cout << " Outer Sort Collapse CopyToSLOW MultiWayMerge FinalCollapse Overall" << std::endl;
  }
  if (KOKKOSKERNELS_VERBOSE) {
    // std::cout << "\t\tTranspose TIME:" << timer1.seconds() << std::endl;
    std::cout << "\t\t" << timer1.seconds() << " ";
  }
  ////////////////////
  // TRANSPOSE MATRIX OVER
  ////////////////////

  ////////////////////////////////////////
  // CALCULATE FLOPS PER OUTER PRODUCT
  ////////////////////////////////////////
  timer1.reset();
  typedef Kokkos::View<size_t *, mySlowMemory> size_t_view_t;
  size_t_view_t flop_per_row(Kokkos::view_alloc(Kokkos::WithoutInitializing, "flops per row"), b_row_cnt);
  FlopsPerRowOuter<row_lno_temp_work_view_t, const_b_lno_row_view_t, size_t_view_t> fpr(b_row_cnt, transpose_col_xadj,
                                                                                        row_mapB, flop_per_row);
  Kokkos::parallel_for(Kokkos::RangePolicy<MyExecSpace>(0, b_row_cnt), fpr);
  MyExecSpace().fence();

  KokkosKernels::Impl::exclusive_parallel_prefix_sum<size_t_view_t, MyExecSpace>(b_row_cnt + 1, flop_per_row);

  auto num_flops   = Kokkos::subview(flop_per_row, b_row_cnt);
  auto h_num_flops = Kokkos::create_mirror_view(num_flops);
  Kokkos::deep_copy(h_num_flops, num_flops);
  MyExecSpace().fence();
  size_t num_required_flops = h_num_flops();

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << timer1.seconds() << " ";
    // std::cout << "\t\tnum_required_flops:" << num_required_flops <<
    // std::endl; std::cout << "\t\tFlopsPerRowOuter TIME:" << timer1.seconds()
    // << std::endl;
  }

  timer1.reset();

  ////////////////////////////////////////
  // DETERMINE BEGINNING AND END of EACH BLOCK
  ////////////////////////////////////////
  // Now create the blocks
  size_t num_blocks = num_required_flops / block_size;
  if (num_required_flops % block_size > 0) ++num_blocks;

  std::vector<size_t> block_xadj((num_blocks + 1) * 2);
  nnz_lno_t total_num_blocks = 0;
  size_t current_block_begin = 0;
  size_t current_block_end   = current_block_begin + block_size;
  block_xadj[0]              = 0;
  // identify blocks
  for (nnz_lno_t i = 1; i < this->b_row_cnt; ++i) {
    if (flop_per_row(i) > current_block_end) {
      block_xadj[++total_num_blocks] = i - 1;

      current_block_begin = flop_per_row(i - 1);
      current_block_end   = current_block_begin + block_size;
    }
  }
  block_xadj[++total_num_blocks] = this->b_row_cnt;
  if (KOKKOSKERNELS_VERBOSE) {
    // std::cout << "\t\tLower Bound num_blocks:" << num_blocks << "
    // total_num_blocks:" << total_num_blocks<< std::endl; std::cout <<
    // "\t\tCALCULATE BLOCKS:" << timer1.seconds() << std::endl;
    std::cout << timer1.seconds() << " ";
  }

  ////////////////////////////////////////
  // ALLOCATE FAST MEMORY TRIPLETS
  ////////////////////////////////////////
  timer1.reset();
  typedef Kokkos::View<Triplet *, myFastMemory> fast_triplet_view_t;
  fast_triplet_view_t fast_memory_triplets(
      // Kokkos::view_alloc(Kokkos::WithoutInitializing,
      "triplets"
      //    )
      ,
      block_size);
  fast_triplet_view_t collapsed_fast_memory_triplets(
      // Kokkos::view_alloc(Kokkos::WithoutInitializing,
      "triplets"
      //    )
      ,
      block_size);

  MyExecSpace().fence();
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << timer1.seconds() << " ";
    // std::cout << "\t\tAllocation WITH FIRST TOUCH TIME:" << timer1.seconds()
    // << std::endl<< std::endl;
  }

  typedef Kokkos::View<Triplet *, mySlowMemory> slow_triplet_view_t;
  std::vector<slow_triplet_view_t> host_triplet_arrays(total_num_blocks);
  size_t overall_size = 0;

  double outerproducttime = 0, sorttime = 0, collapse_time = 0, copy_to_slow_mem_time = 0;
  for (nnz_lno_t bi = 0; bi < total_num_blocks; ++bi) {
    timer1.reset();
    nnz_lno_t begin        = block_xadj[bi];
    nnz_lno_t end          = block_xadj[bi + 1];
    size_t num_block_flops = flop_per_row[end] - flop_per_row[begin];
    // std::cout << "\t\tBLOCK:" << bi << " Begin:" << begin << " End:" << end
    // << std::endl;
    ////////////////////////////////////////
    // OUTER PRODUCT FOR BLOCK BI
    ////////////////////////////////////////
    OuterProduct<row_lno_temp_work_view_t, nnz_lno_temp_work_view_t, scalar_temp_work_view_t, const_b_lno_row_view_t,
                 const_b_lno_nnz_view_t, const_b_scalar_nnz_view_t, size_t_view_t>
        outer_product(begin, end, team_row_chunk_size, transpose_col_xadj, transpose_col_adj, tranpose_vals, row_mapB,
                      entriesB, valsB, flop_per_row, fast_memory_triplets);

    if (this->use_dynamic_schedule)
      Kokkos::parallel_for(
          dynamic_team_policy_t(b_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
          outer_product);
    else
      Kokkos::parallel_for(
          team_policy_t(b_row_cnt / team_row_chunk_size + 1, suggested_team_size, suggested_vector_size),
          outer_product);

    MyExecSpace().fence();
    if (KOKKOSKERNELS_VERBOSE) {
      outerproducttime += timer1.seconds();
      // std::cout << "\t\tOuter Product TIME:" << timer1.seconds() <<
      // std::endl;
    }

    ////////////////////////////////////////
    // SORT TRIPLETS IN FAST MEMORY
    ////////////////////////////////////////
    //
    // std::string filename = char (bi  + '0') + ".triplets";

    // BELOW CODE writes the blocks to a file.
    /*
    std::fstream ff;
    std::stringstream ss;
    ss << bi;
    std::string filename = ss.str();
    filename += ".triplets";
    ff.open(filename.c_str(), std::fstream::out);
    ff << num_block_flops << std::endl;
    for (size_t i = 0; i < num_block_flops; ++i){
      ff << fast_memory_triplets[i].src << " " << fast_memory_triplets[i].dst <<
    " " << fast_memory_triplets[i].val << std::endl;
    }
    ff.close();
     */

    timer1.reset();
    this->sort_triplets(fast_memory_triplets, num_block_flops);
    MyExecSpace().fence();
    if (KOKKOSKERNELS_VERBOSE) {
      sorttime += timer1.seconds();
      // std::cout << "\t\tTriplet Sort Time:" << timer1.seconds() << std::endl;
    }

    ////////////////////////////////////////
    // COLLAPSE TRIPLETS, MERGE SIMILAR ONES
    ////////////////////////////////////////
    timer1.reset();
    size_type outsize =
        this->collapse_triplets_omp(fast_memory_triplets, num_block_flops, collapsed_fast_memory_triplets);
    overall_size += outsize;
    MyExecSpace().fence();
    if (KOKKOSKERNELS_VERBOSE) {
      // std::cout << "\t\toutsize:" << outsize << std::endl;
      collapse_time += timer1.seconds();
      // std::cout << "\t\tTriplet Collapse Time:" << timer1.seconds() <<
      // std::endl;
    }

    ////////////////////////////////////////
    // COPY TO SLOW MEMORY FOR COLLAPSED TRIPLETS
    ////////////////////////////////////////
    timer1.reset();

    host_triplet_arrays[bi] = slow_triplet_view_t(Kokkos::view_alloc(Kokkos::WithoutInitializing, "thv"), outsize);
    KokkosKernels::Impl::kk_copy_vector<fast_triplet_view_t, slow_triplet_view_t, MyExecSpace>(
        outsize, collapsed_fast_memory_triplets, host_triplet_arrays[bi]);
    MyExecSpace().fence();
    if (KOKKOSKERNELS_VERBOSE) {
      copy_to_slow_mem_time += timer1.seconds();
      // std::cout << "\t\tTriplet Copy To Slow Memory Time:" <<
      // timer1.seconds()  << std::endl << std::endl;
    }
  }
  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << outerproducttime << " ";
    std::cout << sorttime << " ";
    std::cout << collapse_time << " ";
    std::cout << copy_to_slow_mem_time << " ";

    // std::cout << "\t\tOuter Product TIME:" << outerproducttime << std::endl;
    // std::cout << "\t\tTriplet Sort Time:" << sorttime << std::endl;
    // std::cout << "\t\tTriplet Collapse Time:" << collapse_time << std::endl;
    // std::cout << "\t\tTriplet Copy To Slow Memory Time:" <<
    // copy_to_slow_mem_time  << std::endl << std::endl;
  }

  ////////////////////////////////////////
  // ALLOCATE MEMORY FOR ALL TRIPLETS IN SLOW
  // AND MULTI-WAY-MERGE THEM
  ////////////////////////////////////////
  timer1.reset();
  slow_triplet_view_t output_triplets("output_triplets", overall_size);

  this->merge_triplets_on_slow_memory(&(host_triplet_arrays[0]), total_num_blocks, overall_size, output_triplets);
  if (KOKKOSKERNELS_VERBOSE) {
    // std::cout << "\t\tMultiway Merge Time:" << timer1.seconds() << std::endl;
    std::cout << timer1.seconds() << " ";
  }
  timer1.reset();

  ////////////////////////////////////////
  // COLLAPSE TRIPLETS, MERGE SAME ONES
  ////////////////////////////////////////
  // triplet_host_view_t
  // collapsed_output_triplets("collapsed_output_triplets",overall_size);

  size_type outsize = this->final_collapse_triplets_omp(output_triplets, overall_size, rowmapC_, entriesC_, valuesC_);
  if (KOKKOSKERNELS_VERBOSE) {
    // std::cout << "\t\tOutsize:" << outsize << " overall_size:" <<
    // overall_size<<  std::endl; std::cout << "\t\tFinal Collapse:" <<
    // timer1.seconds() << std::endl;
    std::cout << timer1.seconds() << " ";
  }

  if (KOKKOSKERNELS_VERBOSE) {
    std::cout << timer_all.seconds() << std::endl;
    // std::cout << "\t\tNumeric TIME:" << timer_all.seconds() << std::endl;
  }

  std::cout << "\t\tLower Bound num_blocks:" << num_blocks << " total_num_blocks:" << total_num_blocks << std::endl;
}
#else
template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void KokkosSPGEMM<
    HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
    b_scalar_nnz_view_t_>::KokkosSPGEMM_numeric_outer(c_row_view_t& /*rowmapC_*/, c_lno_nnz_view_t& /*entriesC_*/,
                                                      c_scalar_nnz_view_t& /*valuesC_*/,
                                                      KokkosKernels::Impl::ExecSpaceType /*my_exec_space_*/) {
  throw std::runtime_error("Cannot run outer product. ENABLE openmp and outer product to run\n");
}
#endif
#else
template <typename HandleType, typename a_row_view_t_, typename a_lno_nnz_view_t_, typename a_scalar_nnz_view_t_,
          typename b_lno_row_view_t_, typename b_lno_nnz_view_t_, typename b_scalar_nnz_view_t_>
template <typename c_row_view_t, typename c_lno_nnz_view_t, typename c_scalar_nnz_view_t>
void KokkosSPGEMM<
    HandleType, a_row_view_t_, a_lno_nnz_view_t_, a_scalar_nnz_view_t_, b_lno_row_view_t_, b_lno_nnz_view_t_,
    b_scalar_nnz_view_t_>::KokkosSPGEMM_numeric_outer(c_row_view_t& /*rowmapC_*/, c_lno_nnz_view_t& /*entriesC_*/,
                                                      c_scalar_nnz_view_t& /*valuesC_*/,
                                                      KokkosKernels::Impl::ExecSpaceType /*my_exec_space_*/) {
  throw std::runtime_error("Cannot run outer product. ENABLE openmp and outer product to run\n");
}
#endif

}  // namespace Impl

}  // namespace KokkosSparse
