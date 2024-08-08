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

#ifndef _KOKKOSGSIMP_HPP
#define _KOKKOSGSIMP_HPP

#include "KokkosKernels_Utils.hpp"
#include "KokkosSparse_Utils.hpp"
#include <Kokkos_Core.hpp>
#include <Kokkos_Bitset.hpp>
#include "KokkosGraph_Distance1Color.hpp"
#include "KokkosKernels_Uniform_Initialized_MemoryPool.hpp"
#include "KokkosKernels_BitUtils.hpp"
#include "KokkosKernels_SimpleUtils.hpp"
#include "KokkosSparse_SortCrs.hpp"

// FOR DEBUGGING
#include "KokkosBlas1_nrm2.hpp"

namespace KokkosSparse {
namespace Impl {

template <typename HandleType, typename lno_row_view_t_, typename lno_nnz_view_t_, typename scalar_nnz_view_t_,
          KokkosSparse::SparseMatrixFormat format = KokkosSparse::SparseMatrixFormat::CRS>
class PointGaussSeidel {
 public:
  typedef lno_row_view_t_ in_lno_row_view_t;
  typedef lno_nnz_view_t_ in_lno_nnz_view_t;
  typedef scalar_nnz_view_t_ in_scalar_nnz_view_t;

  typedef typename HandleType::HandleExecSpace MyExecSpace;
  typedef typename HandleType::HandleTempMemorySpace MyTempMemorySpace;
  typedef typename HandleType::HandlePersistentMemorySpace MyPersistentMemorySpace;

  typedef typename in_lno_row_view_t::non_const_value_type row_lno_t;

  typedef typename HandleType::size_type size_type;
  typedef typename HandleType::nnz_lno_t nnz_lno_t;
  typedef typename HandleType::nnz_scalar_t nnz_scalar_t;

  typedef typename in_lno_row_view_t::const_type const_lno_row_view_t;
  typedef typename in_lno_row_view_t::non_const_type non_const_lno_row_view_t;

  typedef typename lno_nnz_view_t_::const_type const_lno_nnz_view_t;
  typedef typename lno_nnz_view_t_::non_const_type non_const_lno_nnz_view_t;

  typedef typename scalar_nnz_view_t_::const_type const_scalar_nnz_view_t;
  typedef typename scalar_nnz_view_t_::non_const_type non_const_scalar_nnz_view_t;

  typedef typename HandleType::row_lno_temp_work_view_t row_lno_temp_work_view_t;
  typedef typename HandleType::row_lno_persistent_work_view_t row_lno_persistent_work_view_t;
  typedef
      typename HandleType::row_lno_persistent_work_host_view_t row_lno_persistent_work_host_view_t;  // Host view type

  typedef typename HandleType::nnz_lno_temp_work_view_t nnz_lno_temp_work_view_t;
  typedef typename HandleType::nnz_lno_persistent_work_view_t nnz_lno_persistent_work_view_t;
  typedef
      typename HandleType::nnz_lno_persistent_work_host_view_t nnz_lno_persistent_work_host_view_t;  // Host view type

  typedef typename HandleType::scalar_temp_work_view_t scalar_temp_work_view_t;
  typedef typename HandleType::scalar_persistent_work_view2d_t scalar_persistent_work_view2d_t;
  typedef typename HandleType::scalar_persistent_work_view_t scalar_persistent_work_view_t;

  typedef Kokkos::RangePolicy<MyExecSpace> range_policy_t;
  typedef typename HandleType::GraphColoringHandleType::color_view_t color_view_t;
  typedef typename HandleType::GraphColoringHandleType::color_t color_t;
  typedef Kokkos::Bitset<MyExecSpace> bitset_t;
  typedef Kokkos::ConstBitset<MyExecSpace> const_bitset_t;

  typedef Kokkos::TeamPolicy<MyExecSpace> team_policy_t;
  typedef typename team_policy_t::member_type team_member_t;

  struct BlockTag {};
  struct BigBlockTag {};
  struct LongRowTag {};

  typedef Kokkos::TeamPolicy<BlockTag, MyExecSpace> block_apply_team_policy_t;
  typedef Kokkos::TeamPolicy<BigBlockTag, MyExecSpace> bigblock_apply_team_policy_t;
  typedef Kokkos::RangePolicy<LongRowTag, MyExecSpace> longrow_apply_range_policy_t;
  typedef Kokkos::TeamPolicy<LongRowTag, MyExecSpace> longrow_apply_team_policy_t;
  typedef KokkosKernels::Impl::UniformMemoryPool<MyTempMemorySpace, nnz_scalar_t> pool_memory_space;

  typedef typename KokkosSparse::Impl::MatrixRowIndex<format, nnz_lno_t, size_type> RowIndex;

 private:
  HandleType* handle;

  // Get the specialized PointGaussSeidel handle from the main handle
  typename HandleType::PointGaussSeidelHandleType* get_gs_handle() {
    auto gsHandle = dynamic_cast<typename HandleType::PointGaussSeidelHandleType*>(this->handle->get_gs_handle());
    if (!gsHandle) {
      throw std::runtime_error(
          "PointGaussSeidel: GS handle has not been created, or is set up for "
          "Cluster GS.");
    }
    return gsHandle;
  }

  nnz_lno_t num_rows, num_cols;

  const_lno_row_view_t row_map;
  const_lno_nnz_view_t entries;
  const_scalar_nnz_view_t values;

  const_scalar_nnz_view_t given_inverse_diagonal;

  bool have_diagonal_given;
  bool is_symmetric;

  // Batch size for column applies. Used as a stack array size, so must be a
  // compile-time constant.
  static constexpr nnz_lno_t apply_batch_size = 8;

 public:
  struct PSGS {
    row_lno_persistent_work_view_t _xadj;
    nnz_lno_persistent_work_view_t _adj;      // CSR storage of the graph.
    scalar_persistent_work_view_t _adj_vals;  // CSR storage of the graph.

    scalar_persistent_work_view2d_t _Xvector /*output*/;
    scalar_persistent_work_view2d_t _Yvector;

    scalar_persistent_work_view_t _permuted_inverse_diagonal;

    nnz_scalar_t omega;

    scalar_persistent_work_view_t _long_row_x;  // Results of simple Ax matvec over long rows.
    nnz_lno_t _long_row_col;                    // Which X/Y column is now being processed for
                                                // long rows.
    nnz_lno_t _color_set_begin;                 //(only used for long rows): where the current
                                                // set of rows begins
    nnz_lno_t _long_row_par;

    PSGS(row_lno_persistent_work_view_t xadj_, nnz_lno_persistent_work_view_t adj_,
         scalar_persistent_work_view_t adj_vals_, scalar_persistent_work_view2d_t Xvector_,
         scalar_persistent_work_view2d_t Yvector_, nnz_lno_persistent_work_view_t /* color_adj_ */, nnz_scalar_t omega_,
         scalar_persistent_work_view_t permuted_inverse_diagonal_)
        : _xadj(xadj_),
          _adj(adj_),
          _adj_vals(adj_vals_),
          _Xvector(Xvector_),
          _Yvector(Yvector_),
          _permuted_inverse_diagonal(permuted_inverse_diagonal_),
          omega(omega_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t ii) const {
      size_type row_begin                = _xadj(ii);
      size_type row_end                  = _xadj(ii + 1);
      nnz_scalar_t sum[apply_batch_size] = {0};
      nnz_lno_t num_vecs                 = _Xvector.extent(1);
      for (nnz_lno_t batch_start = 0; batch_start < num_vecs; batch_start += apply_batch_size) {
        nnz_lno_t this_batch_size = apply_batch_size;
        if (batch_start + this_batch_size >= num_vecs) this_batch_size = num_vecs - batch_start;
        // the current batch of columns given by: batch_start, this_batch_size
        for (nnz_lno_t i = 0; i < this_batch_size; i++) sum[i] = _Yvector(ii, batch_start + i);
        for (size_type adjind = row_begin; adjind < row_end; ++adjind) {
          nnz_lno_t colIndex = _adj(adjind);
          nnz_scalar_t val   = _adj_vals(adjind);
          for (nnz_lno_t i = 0; i < this_batch_size; i++) sum[i] -= val * _Xvector(colIndex, batch_start + i);
        }
        nnz_scalar_t invDiagonalVal = _permuted_inverse_diagonal(ii);
        for (nnz_lno_t i = 0; i < this_batch_size; i++)
          _Xvector(ii, batch_start + i) += omega * sum[i] * invDiagonalVal;
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const LongRowTag&, const nnz_lno_t i) const {
      nnz_lno_t row         = _color_set_begin + i / _long_row_par;
      nnz_lno_t chunk       = i % _long_row_par;
      size_type row_begin   = _xadj(row);
      size_type row_end     = _xadj(row + 1);
      size_type chunk_begin = row_begin + (row_end - row_begin) * chunk / _long_row_par;
      size_type chunk_end   = row_begin + (row_end - row_begin) * (chunk + 1) / _long_row_par;
      if (chunk_end > row_end) chunk_end = row_end;
      nnz_scalar_t localSum{};
      for (size_type j = chunk_begin; j < chunk_end; j++) {
        nnz_lno_t colIndex = _adj(j);
        localSum += _adj_vals(j) * _Xvector(colIndex, _long_row_col);
      }
      Kokkos::atomic_add(&_long_row_x(row - _color_set_begin), localSum);
    }
  };

  struct Team_PSGS {
    row_lno_persistent_work_view_t _xadj;
    nnz_lno_persistent_work_view_t _adj;      // CSR storage of the graph.
    scalar_persistent_work_view_t _adj_vals;  // CSR storage of the graph.

    scalar_persistent_work_view2d_t _Xvector /*output*/;
    scalar_persistent_work_view2d_t _Yvector;
    nnz_lno_t _color_set_begin;
    nnz_lno_t _color_set_end;

    scalar_persistent_work_view_t _permuted_inverse_diagonal;
    nnz_lno_t block_size;
    nnz_lno_t team_work_size;
    const size_t shared_memory_size;

    int suggested_team_size;
    const size_t thread_shared_memory_scalar_size;
    int vector_size;
    const pool_memory_space pool;
    const nnz_lno_t num_max_vals_in_l1, num_max_vals_in_l2;
    bool is_backward;

    nnz_scalar_t omega;

    typedef typename KokkosKernels::Impl::array_sum_reduce<nnz_scalar_t, apply_batch_size> batch_sum;

    nnz_lno_persistent_work_view_t _long_rows;
    scalar_persistent_work_view_t _long_row_x;
    nnz_lno_t _long_row_col;  // Which X/Y column is now being processed for
                              // long rows.
    nnz_lno_t _long_row_par;

    Team_PSGS(row_lno_persistent_work_view_t xadj_, nnz_lno_persistent_work_view_t adj_,
              scalar_persistent_work_view_t adj_vals_, scalar_persistent_work_view2d_t Xvector_,
              scalar_persistent_work_view2d_t Yvector_, nnz_lno_t color_set_begin, nnz_lno_t color_set_end,
              scalar_persistent_work_view_t permuted_inverse_diagonal_, pool_memory_space pms,
              nnz_lno_t _num_max_vals_in_l1 = 0, nnz_lno_t _num_max_vals_in_l2 = 0,
              nnz_scalar_t omega_ = Kokkos::ArithTraits<nnz_scalar_t>::one(),

              nnz_lno_t block_size_ = 1, nnz_lno_t team_work_size_ = 1, size_t shared_memory_size_ = 16,
              int suggested_team_size_ = 1, int vector_size_ = 1)
        : _xadj(xadj_),
          _adj(adj_),
          _adj_vals(adj_vals_),
          _Xvector(Xvector_),
          _Yvector(Yvector_),
          _color_set_begin(color_set_begin),
          _color_set_end(color_set_end),
          _permuted_inverse_diagonal(permuted_inverse_diagonal_),
          block_size(block_size_),
          team_work_size(team_work_size_),
          shared_memory_size(shared_memory_size_),
          suggested_team_size(suggested_team_size_),
          thread_shared_memory_scalar_size(((shared_memory_size / suggested_team_size / 8) * 8) / sizeof(nnz_scalar_t)),
          vector_size(vector_size_),
          pool(pms),
          num_max_vals_in_l1(_num_max_vals_in_l1),
          num_max_vals_in_l2(_num_max_vals_in_l2),
          is_backward(false),
          omega(omega_) {}

    // Do a Gauss-Seidel step on a single row, for X/Y columns
    // colStart:colStart+N-1 (inclusive) Specializing this on the batch size
    // allows the best reuse of matrix accesses, while also using the correct
    // width array_sum_reduce.
    template <int N>
    KOKKOS_INLINE_FUNCTION void runColBatch(const team_member_t& teamMember, nnz_lno_t row, nnz_lno_t colStart) const {
      typedef KokkosKernels::Impl::array_sum_reduce<nnz_scalar_t, N> reducer;
      size_type row_begin = _xadj(row);
      size_type row_end   = _xadj(row + 1);
      reducer sum;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(teamMember, row_end - row_begin),
          [&](size_type i, reducer& lsum) {
            size_type adjind   = row_begin + i;
            nnz_lno_t colIndex = _adj(adjind);
            nnz_scalar_t val   = _adj_vals(adjind);
            for (int j = 0; j < N; j++) lsum.data[j] += val * _Xvector(colIndex, colStart + j);
          },
          sum);
      Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
        nnz_scalar_t invDiagonalVal = _permuted_inverse_diagonal(row);
        for (int i = 0; i < N; i++) {
          _Xvector(row, colStart + i) += omega * (_Yvector(row, colStart + i) - sum.data[i]) * invDiagonalVal;
        }
      });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t& teamMember) const {
      nnz_lno_t row = teamMember.league_rank() * teamMember.team_size() + teamMember.team_rank() + _color_set_begin;
      if (row >= _color_set_end) return;
      nnz_lno_t num_vecs = _Xvector.extent(1);
      for (nnz_lno_t batch_start = 0; batch_start < num_vecs;) {
        switch (num_vecs - batch_start) {
#define COL_BATCH_CASE(n)                         \
  case n:                                         \
    runColBatch<n>(teamMember, row, batch_start); \
    batch_start += n;                             \
    break;
          COL_BATCH_CASE(1)
          COL_BATCH_CASE(2)
          COL_BATCH_CASE(3)
          COL_BATCH_CASE(4)
          COL_BATCH_CASE(5)
          COL_BATCH_CASE(6)
          COL_BATCH_CASE(7)
#undef COL_BATCH_CASE
          default: runColBatch<8>(teamMember, row, batch_start); batch_start += 8;
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const BigBlockTag&, const team_member_t& teamMember) const {
      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + _color_set_begin;
      const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, _color_set_end);
      // get the shared memory and shift it based on the thread index so that
      // each thread has private memory.
      nnz_scalar_t* all_shared_memory = (nnz_scalar_t*)(teamMember.team_shmem().get_shmem(shared_memory_size));

      all_shared_memory += thread_shared_memory_scalar_size * teamMember.team_rank();

      // store the diagonal positions, because we want to update them on shared
      // memory if we update them on global memory.
      nnz_lno_t* diagonal_positions = (nnz_lno_t*)all_shared_memory;
      all_shared_memory             = (nnz_scalar_t*)(((nnz_lno_t*)all_shared_memory) + ((block_size / 8) + 1) * 8);

      nnz_scalar_t* all_global_memory = NULL;

      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t ii) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, block_size),
                             [&](nnz_lno_t i) { diagonal_positions[i] = -1; });

        RowIndex row(block_size, _xadj(ii), _xadj(ii + 1));
        nnz_lno_t l1_val_size = row.size() * block_size, l2_val_size = 0;
        // if the current row size is larger than shared memory size,
        // than allocate l2 vector.
        if (row.size() > num_max_vals_in_l1) {
          volatile nnz_scalar_t* tmp = NULL;
          while (tmp == NULL) {
            Kokkos::single(
                Kokkos::PerThread(teamMember),
                [&](volatile nnz_scalar_t*& memptr) { memptr = (volatile nnz_scalar_t*)(pool.allocate_chunk(ii)); },
                tmp);
          }
          all_global_memory = (nnz_scalar_t*)tmp;
          l1_val_size       = num_max_vals_in_l1 * block_size;
          l2_val_size       = (row.size() * block_size - l1_val_size);
        }
        for (nnz_lno_t vec = 0; vec < (nnz_lno_t)_Xvector.extent(1); vec++) {
          // bring values to l1 vector
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, l1_val_size), [&](nnz_lno_t i) {
            size_type adjind   = i / block_size + row.begin();
            nnz_lno_t colIndex = _adj(adjind);

            if (colIndex == ii) {
              diagonal_positions[i % block_size] = i;
            }
            all_shared_memory[i] = _Xvector(colIndex * block_size + i % block_size, vec);
          });
          // bring values to l2 vector.
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, l2_val_size), [&](nnz_lno_t k) {
            nnz_lno_t i = l1_val_size + k;

            size_type adjind   = i / block_size + row.begin();
            nnz_lno_t colIndex = _adj(adjind);

            if (colIndex == ii) {
              diagonal_positions[i % block_size] = i;
            }
            all_global_memory[k] = _Xvector(colIndex * block_size + i % block_size, vec);
          });

          // sequentially solve in the block.
          // this respects backward and forward sweeps.
          for (int m = 0; m < block_size; ++m) {
            int i = m;
            if (is_backward) i = block_size - m - 1;
            // first reduce l1 dot product.
            // MD: TODO: if thread dot product is implemented it should be
            // called here.
            nnz_scalar_t product = 0;
            Kokkos::parallel_reduce(
                Kokkos::ThreadVectorRange(teamMember, l1_val_size),
                [&](nnz_lno_t colind, nnz_scalar_t& valueToUpdate) {
                  const size_type val_idx = row.value(colind / block_size, i, colind % block_size);
                  valueToUpdate += all_shared_memory[colind] * _adj_vals(val_idx);
                },
                product);
            // l2 dot product.
            // MD: TODO: if thread dot product is implemented, it should be
            // called here again.
            nnz_scalar_t product2 = 0;
            Kokkos::parallel_reduce(
                Kokkos::ThreadVectorRange(teamMember, l2_val_size),
                [&](nnz_lno_t colind2, nnz_scalar_t& valueToUpdate) {
                  nnz_lno_t colind        = colind2 + l1_val_size;
                  const size_type val_idx = row.value(colind / block_size, i, colind % block_size);
                  valueToUpdate += all_global_memory[colind2] * _adj_vals(val_idx);
                },
                product2);

            product += product2;
            // update the new vector entries.
            Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
              nnz_lno_t block_row_index   = ii * block_size + i;
              nnz_scalar_t invDiagonalVal = _permuted_inverse_diagonal(block_row_index);
              _Xvector(block_row_index, vec) += omega * (_Yvector(block_row_index, vec) - product) * invDiagonalVal;

              // we need to update the values of the vector entries if they
              // are already brought to shared memory to sync with global
              // memory.
              if (diagonal_positions[i] != -1) {
                if (diagonal_positions[i] < l1_val_size)
                  all_shared_memory[diagonal_positions[i]] = _Xvector(block_row_index, vec);
                else
                  all_global_memory[diagonal_positions[i] - l1_val_size] = _Xvector(block_row_index, vec);
              }
            });
          }

#if KOKKOSSPARSE_IMPL_PRINTDEBUG
          if (/*i == 0 && ii == 1*/ ii == 0 || (block_size == 1 && ii < 2)) {
            std::cout << "In X/Y column " << vec << std::endl;
            std::cout << "\n\n\nrow:" << ii * block_size + i;
            std::cout << "\nneighbors:";
            for (int z = 0; z < int(row_size); ++z) {
              std::cout << _adj(_xadj(ii) + z) << " ";
            }
            std::cout << "\n\nrow-0:X -- all-shared-memory:";
            for (int z = 0; z < int(row_size * block_size); ++z) {
              std::cout << all_shared_memory[z] << " ";
            }
            std::cout << std::endl << "product:" << product << std::endl;
            std::cout << "diagonal" << _permuted_inverse_diagonal(ii * block_size + i) << std::endl;
            std::cout << "_Yvector: " << _Yvector(ii * block_size + i, vec) << std::endl;
            std::cout << "block_row_index:" << ii * block_size + i << " _Xvector(block_row_index): ";
            _Xvector(ii * block_size + i, vec) << ' ';
            std::cout << std::endl;
          }
#endif
        }
        if (row.size() > num_max_vals_in_l1) {
          Kokkos::single(Kokkos::PerThread(teamMember), [&]() { pool.release_chunk(all_global_memory); });
        }
      });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const BlockTag&, const team_member_t& teamMember) const {
      const nnz_lno_t team_row_begin = teamMember.league_rank() * team_work_size + _color_set_begin;
      const nnz_lno_t team_row_end   = KOKKOSKERNELS_MACRO_MIN(team_row_begin + team_work_size, _color_set_end);

      nnz_scalar_t* all_shared_memory = (nnz_scalar_t*)(teamMember.team_shmem().get_shmem(shared_memory_size));

      all_shared_memory += thread_shared_memory_scalar_size * teamMember.team_rank();

      nnz_lno_t* diagonal_positions = (nnz_lno_t*)all_shared_memory;
      all_shared_memory             = (nnz_scalar_t*)(((nnz_lno_t*)all_shared_memory) + ((block_size / 8) + 1) * 8);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(teamMember, team_row_begin, team_row_end), [&](const nnz_lno_t& ii) {
#if KOKKOSSPARSE_IMPL_PRINTDEBUG
        Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
          for (nnz_lno_t i = 0; i < block_size; diagonal_positions[i++] = -1)
            ;
        });
#endif

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, block_size),
                             [&](nnz_lno_t i) { diagonal_positions[i] = -1; });

        RowIndex row(block_size, _xadj(ii), _xadj(ii + 1));
        // offset in adj_vals of the first row in this block
        // number of scalars in each row of this block
        nnz_lno_t scalar_row_size = row.size() * block_size;

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, scalar_row_size), [&](nnz_lno_t i) {
          size_type adjind   = i / block_size + row.begin();
          nnz_lno_t colIndex = _adj(adjind);
          if (colIndex == ii) {
            diagonal_positions[i % block_size] = i;
          }
        });

        for (nnz_lno_t vec = 0; vec < (nnz_lno_t)_Xvector.extent(1); vec++) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(teamMember, scalar_row_size), [&](nnz_lno_t i) {
            size_type adjind     = i / block_size + row.begin();
            nnz_lno_t colIndex   = _adj(adjind);
            all_shared_memory[i] = _Xvector(colIndex * block_size + i % block_size, vec);
          });

          for (int m = 0; m < block_size; ++m) {
            int i = m;
            if (is_backward) {
              i = block_size - m - 1;
            }
            nnz_scalar_t product = 0;
            Kokkos::parallel_reduce(
                Kokkos::ThreadVectorRange(teamMember, scalar_row_size),
                [&](nnz_lno_t colind, nnz_scalar_t& valueToUpdate) {
                  const size_type val_idx = row.value(colind / block_size, i, colind % block_size);

                  valueToUpdate += all_shared_memory[colind] * _adj_vals(val_idx);
                },
                product);

            Kokkos::single(Kokkos::PerThread(teamMember), [&]() {
              nnz_lno_t block_row_index   = ii * block_size + i;
              nnz_scalar_t invDiagonalVal = _permuted_inverse_diagonal(block_row_index);
              _Xvector(block_row_index, vec) += omega * (_Yvector(block_row_index, vec) - product) * invDiagonalVal;

              if (diagonal_positions[i] != -1) {
                all_shared_memory[diagonal_positions[i]] = _Xvector(block_row_index, vec);
              }
            });

#if KOKKOSSPARSE_IMPL_PRINTDEBUG
            if (!KokkosKernels::Impl::kk_is_gpu_exec_space<typename team_member_t::execution_space>() &&
                (ii == 0 || (block_size == 1 && ii < 2))) {
              std::cout << "\n\n\nrow:" << ii * block_size + i;
              std::cout << "\nneighbors:";
              for (nnz_lno_t z = 0; z < block_row_size; ++z) {
                std::cout << _adj(_xadj(ii) + z) << " ";
              }

              std::cout << "\n\nrow-0:X -- all-shared-memory:";
              for (nnz_lno_t z = 0; z < scalar_row_size; ++z) {
                std::cout << all_shared_memory[z] << " ";
              }
              std::cout << std::endl << "product:" << product << std::endl;
              std::cout << "diagonal" << _permuted_inverse_diagonal(ii * block_size + i) << std::endl;
              std::cout << "_Yvector" << _Yvector(ii * block_size + i, vec) << std::endl;

              std::cout << std::endl
                        << "block_row_index:" << ii * block_size + i
                        << " _Xvector(block_row_index):" << _Xvector(ii * block_size + i, vec) << std::endl
                        << std::endl
                        << std::endl;
            }
#endif
          }
        }
      });
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const LongRowTag&, const team_member_t& teamMember) const {
      nnz_lno_t row         = _color_set_begin + teamMember.league_rank() / _long_row_par;
      nnz_lno_t chunk       = teamMember.league_rank() % _long_row_par;
      size_type row_begin   = _xadj(row);
      size_type row_end     = _xadj(row + 1);
      size_type chunk_begin = row_begin + (row_end - row_begin) * chunk / _long_row_par;
      size_type chunk_end   = row_begin + (row_end - row_begin) * (chunk + 1) / _long_row_par;
      if (chunk_end > row_end) chunk_end = row_end;
      nnz_scalar_t localSum;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(teamMember, chunk_begin, chunk_end),
          [&](size_type j, nnz_scalar_t& lsum) {
            nnz_lno_t colIndex = _adj(j);
            lsum += _adj_vals(j) * _Xvector(colIndex, _long_row_col);
          },
          localSum);
      Kokkos::single(Kokkos::PerTeam(teamMember),
                     [&]() { Kokkos::atomic_add(&_long_row_x(row - _color_set_begin), localSum); });
    }

    size_t team_shmem_size(int /* team_size */) const { return shared_memory_size; }
  };

  struct LongRowComparator {
    KOKKOS_DEFAULTED_FUNCTION LongRowComparator() = default;
    KOKKOS_INLINE_FUNCTION LongRowComparator(const in_lno_row_view_t& xadj_, nnz_lno_t longRowThreshold_)
        : xadj(xadj_), longRowThreshold(longRowThreshold_) {}

    KOKKOS_INLINE_FUNCTION bool operator()(nnz_lno_t lhs, nnz_lno_t rhs) const {
      int lhsLong = xadj(lhs + 1) - xadj(lhs) >= longRowThreshold;
      int rhsLong = xadj(rhs + 1) - xadj(rhs) >= longRowThreshold;
      if (lhsLong < rhsLong)
        return true;
      else if (lhsLong > rhsLong)
        return false;
      // Either both long or both short, just order ascending by ID
      return lhs < rhs;
    }

    in_lno_row_view_t xadj;
    size_type longRowThreshold;
  };

  // Functor to sort each color set - first by whether 'long row', second by ID.
  // Also populates long_rows_per_color.
  struct SortIntoLongRowsFunctor {
    SortIntoLongRowsFunctor(const in_lno_row_view_t& xadj_, nnz_lno_t longRowThreshold_,
                            const nnz_lno_persistent_work_view_t& color_xadj_,
                            const nnz_lno_persistent_work_view_t& color_adj_,
                            const nnz_lno_persistent_work_view_t& long_rows_per_color_,
                            const nnz_lno_persistent_work_view_t& max_row_length_per_color_)
        : xadj(xadj_),
          longRowThreshold(longRowThreshold_),
          color_xadj(color_xadj_),
          color_adj(color_adj_),
          long_rows_per_color(long_rows_per_color_),
          max_row_length_per_color(max_row_length_per_color_) {}

    KOKKOS_INLINE_FUNCTION void operator()(const team_member_t& t, nnz_lno_t& lmostPerColor) const {
      LongRowComparator comp(xadj, longRowThreshold);
      nnz_lno_t color      = t.league_rank();
      nnz_lno_t colorBegin = color_xadj(color);
      nnz_lno_t colorLen   = color_xadj(color + 1) - colorBegin;
      KokkosKernels::TeamBitonicSort(color_adj.data() + colorBegin, colorLen, t, comp);
      t.team_barrier();
      // Now that the color set is sorted, count how many long rows there were
      nnz_lno_t numLongRows;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(t, colorBegin, colorBegin + colorLen),
          [&](nnz_lno_t i, nnz_lno_t& lnumLongRows) {
            nnz_lno_t row = color_adj(i);
            if (xadj(row + 1) - xadj(row) >= longRowThreshold) lnumLongRows++;
          },
          numLongRows);
      Kokkos::single(Kokkos::PerTeam(t), [&]() {
        long_rows_per_color(color) = numLongRows;
        if (numLongRows > lmostPerColor) lmostPerColor = numLongRows;
      });
      nnz_lno_t max_row_length = 0;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(t, colorBegin, colorBegin + colorLen),
          [&](nnz_lno_t i, nnz_lno_t& lmaxLength) {
            nnz_lno_t row = color_adj(i);
            nnz_lno_t len = xadj(row + 1) - xadj(row);
            if (len > lmaxLength) lmaxLength = len;
          },
          Kokkos::Max<nnz_lno_t>(max_row_length));
      Kokkos::single(Kokkos::PerTeam(t), [&]() { max_row_length_per_color(color) = max_row_length; });
    }

    in_lno_row_view_t xadj;
    size_type longRowThreshold;
    nnz_lno_persistent_work_view_t color_xadj;
    nnz_lno_persistent_work_view_t color_adj;
    nnz_lno_persistent_work_view_t long_rows_per_color;
    nnz_lno_persistent_work_view_t max_row_length_per_color;
  };

  /**
   * \brief constructor
   */

  PointGaussSeidel(HandleType* handle_, nnz_lno_t num_rows_, nnz_lno_t num_cols_, const_lno_row_view_t row_map_,
                   const_lno_nnz_view_t entries_, const_scalar_nnz_view_t values_)
      : handle(handle_),
        num_rows(num_rows_),
        num_cols(num_cols_),
        row_map(row_map_),
        entries(entries_),
        values(values_),
        have_diagonal_given(false),
        is_symmetric(true) {}

  PointGaussSeidel(HandleType* handle_, nnz_lno_t num_rows_, nnz_lno_t num_cols_, const_lno_row_view_t row_map_,
                   const_lno_nnz_view_t entries_, bool is_symmetric_ = true)
      : handle(handle_),
        num_rows(num_rows_),
        num_cols(num_cols_),
        row_map(row_map_),
        entries(entries_),
        values(),
        have_diagonal_given(false),
        is_symmetric(is_symmetric_) {}

  /**
   * \brief constructor
   */
  PointGaussSeidel(HandleType* handle_, nnz_lno_t num_rows_, nnz_lno_t num_cols_, const_lno_row_view_t row_map_,
                   const_lno_nnz_view_t entries_, const_scalar_nnz_view_t values_, bool is_symmetric_)
      : handle(handle_),
        num_rows(num_rows_),
        num_cols(num_cols_),
        row_map(row_map_),
        entries(entries_),
        values(values_),
        have_diagonal_given(false),
        is_symmetric(is_symmetric_) {}

  PointGaussSeidel(HandleType* handle_, nnz_lno_t num_rows_, nnz_lno_t num_cols_, const_lno_row_view_t row_map_,
                   const_lno_nnz_view_t entries_, const_scalar_nnz_view_t values_,
                   const_scalar_nnz_view_t given_inverse_diagonal_, bool is_symmetric_)
      : handle(handle_),
        num_rows(num_rows_),
        num_cols(num_cols_),
        row_map(row_map_),
        entries(entries_),
        values(values_),
        given_inverse_diagonal(given_inverse_diagonal_),
        have_diagonal_given(true),
        is_symmetric(is_symmetric_) {}

  void initialize_symbolic() {
    auto gsHandle                    = get_gs_handle();
    const size_type longRowThreshold = gsHandle->get_long_row_threshold();
    const MyExecSpace my_exec_space  = gsHandle->get_execution_space();
    const int num_streams            = gsHandle->get_num_streams();

    // Validate settings
    if (gsHandle->get_block_size() > 1 && longRowThreshold > 0)
      throw std::runtime_error("Can't use MTGS long row algorithm with blocks.");

    const_lno_row_view_t xadj = this->row_map;
    const_lno_nnz_view_t adj  = this->entries;
    size_type nnz             = adj.extent(0);

#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    Kokkos::Timer timer;
#endif
    // TODO: Pass my_exec_space into KokkosGraph kernels
    typename HandleType::GraphColoringHandleType::color_view_t colors;
    color_t numColors;
    {
      HandleType coloringHandle;
      coloringHandle.create_graph_coloring_handle(gsHandle->get_coloring_algorithm());
      auto gchandle = coloringHandle.get_graph_coloring_handle();
      if (!is_symmetric) {
        if (gchandle->get_coloring_algo_type() == KokkosGraph::COLORING_EB) {
          gchandle->symmetrize_and_calculate_lower_diagonal_edge_list(num_rows, xadj, adj);
          KokkosGraph::Experimental::graph_color_symbolic<HandleType, const_lno_row_view_t, const_lno_nnz_view_t>(
              &coloringHandle, num_rows, num_rows, xadj, adj);
        } else {
          row_lno_temp_work_view_t tmp_xadj;
          nnz_lno_temp_work_view_t tmp_adj;
          KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<const_lno_row_view_t, const_lno_nnz_view_t,
                                                                 row_lno_temp_work_view_t, nnz_lno_temp_work_view_t,
                                                                 MyExecSpace>(num_rows, xadj, adj, tmp_xadj, tmp_adj);
          KokkosGraph::Experimental::graph_color_symbolic<HandleType, row_lno_temp_work_view_t,
                                                          nnz_lno_temp_work_view_t>(&coloringHandle, num_rows, num_rows,
                                                                                    tmp_xadj, tmp_adj);
        }
      } else {
        KokkosGraph::Experimental::graph_color_symbolic<HandleType, const_lno_row_view_t, const_lno_nnz_view_t>(
            &coloringHandle, num_rows, num_rows, xadj, adj);
      }
      colors    = gchandle->get_vertex_colors();
      numColors = gchandle->get_num_colors();
    }
    // Wait for coloring to finish on its stream
    using ColoringExecSpace = typename HandleType::HandleExecSpace;
    ColoringExecSpace().fence();
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    std::cout << "COLORING_TIME:" << timer.seconds() << std::endl;
    timer.reset();
#endif

#if KOKKOSSPARSE_IMPL_RUNSEQUENTIAL
    numColors = num_rows;
    KokkosKernels::Impl::print_1Dview(colors);
    std::cout << "numCol:" << numColors << " numRows:" << num_rows << " cols:" << num_cols << " nnz:" << adj.extent(0)
              << std::endl;
    typename HandleType::GraphColoringHandleType::color_view_t::HostMirror h_colors =
        Kokkos::create_mirror_view(colors);
    for (int i = 0; i < num_rows; ++i) {
      h_colors(i) = i + 1;
    }
    Kokkos::deep_copy(my_exec_space, colors, h_colors);
    my_exec_space.fence();
#endif
    nnz_lno_persistent_work_view_t color_xadj;
    nnz_lno_persistent_work_view_t color_adj;
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    timer.reset();
#endif
    KokkosKernels::Impl::create_reverse_map<typename HandleType::GraphColoringHandleType::color_view_t,
                                            nnz_lno_persistent_work_view_t, MyExecSpace>(
        my_exec_space, num_rows, numColors, colors, color_xadj, color_adj);

#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    my_exec_space.fence();
    std::cout << "CREATE_REVERSE_MAP:" << timer.seconds() << std::endl;
    timer.reset();
#endif

    nnz_lno_persistent_work_host_view_t h_color_xadj = Kokkos::create_mirror_view(color_xadj);
    Kokkos::deep_copy(h_color_xadj, color_xadj);

#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    my_exec_space.fence();
    std::cout << "DEEP_COPY:" << timer.seconds() << std::endl;
    timer.reset();
#endif
    if (longRowThreshold > 0) {
      // Count long rows per color set, and sort color sets so that long rows
      // come after regular rows
      nnz_lno_persistent_work_view_t long_rows_per_color(
          Kokkos::view_alloc(my_exec_space, Kokkos::WithoutInitializing, "long_rows_per_color"), numColors);
      nnz_lno_persistent_work_view_t max_row_length_per_color(
          Kokkos::view_alloc(my_exec_space, Kokkos::WithoutInitializing, "max_row_length_per_color"), numColors);
      nnz_lno_t mostLongRowsInColor = 0;
      SortIntoLongRowsFunctor sortIntoLongRowsFunctor(xadj, longRowThreshold, color_xadj, color_adj,
                                                      long_rows_per_color, max_row_length_per_color);
      int sortLongRowsTeamSize = 1;
      {
        team_policy_t temp(my_exec_space, 1, 1);
        sortLongRowsTeamSize = temp.team_size_recommended(sortIntoLongRowsFunctor, Kokkos::ParallelReduceTag());
      }
      Kokkos::parallel_reduce(team_policy_t(my_exec_space, numColors, sortLongRowsTeamSize), sortIntoLongRowsFunctor,
                              Kokkos::Max<nnz_lno_t>(mostLongRowsInColor));
      auto host_long_rows_per_color = Kokkos::create_mirror_view(long_rows_per_color);
      Kokkos::deep_copy(my_exec_space, host_long_rows_per_color, long_rows_per_color);
      my_exec_space.fence();
      gsHandle->set_long_rows_per_color(host_long_rows_per_color);
      auto host_max_row_length_per_color = Kokkos::create_mirror_view(max_row_length_per_color);
      Kokkos::deep_copy(my_exec_space, host_max_row_length_per_color, max_row_length_per_color);
      my_exec_space.fence();
      gsHandle->set_max_row_length_per_color(host_max_row_length_per_color);
      scalar_persistent_work_view_t long_row_x(
          Kokkos::view_alloc(my_exec_space, Kokkos::WithoutInitializing, "long_row_x"), mostLongRowsInColor);
      gsHandle->set_long_row_x(long_row_x);
    } else {
      // Just sort rows by ID.
      KokkosSparse::sort_crs_graph<MyExecSpace, decltype(color_xadj), decltype(color_adj)>(my_exec_space, color_xadj,
                                                                                           color_adj);
    }
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    my_exec_space.fence();
    std::cout << "SORT_TIME:" << timer.seconds() << std::endl;
    timer.reset();
#endif

    row_lno_persistent_work_view_t permuted_xadj("new xadj", num_rows + 1);
    nnz_lno_persistent_work_view_t old_to_new_map("old_to_new_index_", num_rows);
    nnz_lno_persistent_work_view_t permuted_adj("newadj_", nnz);

    Kokkos::parallel_for("KokkosSparse::PointGaussSeidel::create_permuted_xadj",
                         range_policy_t(my_exec_space, 0, num_rows),
                         create_permuted_xadj(color_adj, xadj, permuted_xadj, old_to_new_map));
    // std::cout << "create_permuted_xadj" << std::endl;

#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    my_exec_space.fence();
    std::cout << "CREATE_PERMUTED_XADJ:" << timer.seconds() << std::endl;

    timer.reset();
#endif

    KokkosKernels::Impl::inclusive_parallel_prefix_sum<row_lno_persistent_work_view_t, MyExecSpace>(
        my_exec_space, num_rows + 1, permuted_xadj);

#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    my_exec_space.fence();
    std::cout << "INCLUSIVE_PPS:" << timer.seconds() << std::endl;
    timer.reset();
#endif

    Kokkos::parallel_for("KokkosSparse::PointGaussSeidel::fill_matrix_symbolic",
                         range_policy_t(my_exec_space, 0, num_rows),
                         fill_matrix_symbolic(num_rows, color_adj, xadj, adj,
                                              // adj_vals,
                                              permuted_xadj, permuted_adj,
                                              // newvals_,
                                              old_to_new_map));

#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    my_exec_space.fence();
    std::cout << "SYMBOLIC_FILL:" << timer.seconds() << std::endl;
    timer.reset();
#endif

    nnz_lno_t block_size = get_gs_handle()->get_block_size();

    // MD: if block size is larger than 1;
    // the algorithm copies the vector entries into shared memory and reuses
    // this small shared memory for vector entries.
    if (block_size > 1) {
      // first calculate max row size.
      size_type max_row_size = 0;
      KokkosKernels::Impl::kk_view_reduce_max_row_size<size_type, MyExecSpace>(
          my_exec_space, num_rows, permuted_xadj.data(), permuted_xadj.data() + 1, max_row_size);

      nnz_lno_t brows = permuted_xadj.extent(0) - 1;
      size_type bnnz  = permuted_adj.extent(0) * block_size * block_size;

      int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
      int suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);
      size_t shmem_size_to_use  = this->handle->get_shmem_size();

      // MD: now we calculate how much memory is needed for shared memory.
      // we have two-level vectors: as in spgemm hashmaps.
      // we try to fit everything into shared memory.
      // if they fit, we can use BlockTeam function in Team_SGS functor.
      // on CPUs, we make L1 vector big enough so that it will always hold it.
      // on GPUs, we have a upper bound for shared memory:
      // handle->get_shmem_size(): this is set to 32128 bytes. If things do not
      // fit into shared memory, we allocate vectors in global memory and run
      // BigBlockTeam in Team_SGS functor.
      size_t level_1_mem =
          max_row_size * block_size * sizeof(nnz_scalar_t) + ((block_size / 8) + 1) * 8 * sizeof(nnz_lno_t);
      level_1_mem                = suggested_team_size * level_1_mem;
      size_t level_2_mem         = 0;
      nnz_lno_t num_values_in_l1 = max_row_size;
      nnz_lno_t num_values_in_l2 = 0;
      nnz_lno_t num_big_rows     = 0;

      if (!KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>()) {
        // again, if it is on CPUs, we make L1 as big as we need.
        size_t l1mem = 1;
        while (l1mem < level_1_mem) {
          l1mem *= 2;
        }
        gsHandle->set_level_1_mem(l1mem);
        level_1_mem = l1mem;
        level_2_mem = 0;
      } else {
        // on GPUs set the L1 size to max shmem and calculate how much we need
        // for L2. we try to shift with 8 always because of the errors we
        // experience with memory shifts on GPUs.
        level_1_mem      = shmem_size_to_use;
        num_values_in_l1 = (shmem_size_to_use / suggested_team_size - ((block_size / 8) + 1) * 8 * sizeof(nnz_lno_t)) /
                           sizeof(nnz_scalar_t) / block_size;
        if (((block_size / 8) + 1) * 8 * sizeof(nnz_lno_t) > shmem_size_to_use / suggested_team_size)
          throw "Shared memory size is to small for the given block size\n";
        if (num_values_in_l1 >= (nnz_lno_t)(max_row_size)) {
          num_values_in_l2 = 0;
          level_2_mem      = 0;
          num_big_rows     = 0;
        } else {
          num_values_in_l2 = max_row_size - num_values_in_l1;
          level_2_mem      = num_values_in_l2 * block_size * sizeof(nnz_scalar_t);
          // std::cout << "level_2_mem:" << level_2_mem << std::endl;
          size_t l2mem = 1;
          while (l2mem < level_2_mem) {
            l2mem *= 2;
          }
          level_2_mem = l2mem;
          // std::cout << "level_2_mem:" << level_2_mem << std::endl;

          size_type num_large_rows = 0;
          KokkosSparse::Impl::kk_reduce_numrows_larger_than_threshold<row_lno_persistent_work_view_t, MyExecSpace>(
              my_exec_space, brows, permuted_xadj, num_values_in_l1, num_large_rows);
          num_big_rows =
              KOKKOSKERNELS_MACRO_MIN(num_large_rows, (size_type)(my_exec_space.concurrency() / suggested_vector_size));
          // std::cout << "num_big_rows:" << num_big_rows << std::endl;

          if (KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>()) {
            // check if we have enough memory for this. lower the concurrency if
            // we do not have enugh memory.
            size_t free_byte;
            size_t total_byte;
            KokkosKernels::Impl::kk_get_free_total_memory<typename pool_memory_space::memory_space>(
                free_byte, total_byte, num_streams);
            size_t required_size = size_t(num_big_rows) * level_2_mem;
            if (required_size + num_big_rows * sizeof(int) > free_byte) {
              num_big_rows = ((((free_byte - num_big_rows * sizeof(int)) * 0.8) / 8) * 8) / level_2_mem;
            }
            {
              nnz_lno_t min_chunk_size = 1;
              while (min_chunk_size * 2 <= num_big_rows) {
                min_chunk_size *= 2;
              }
              num_big_rows = min_chunk_size;
            }
          }
        }
      }

      gsHandle->set_level_1_mem(level_1_mem);
      gsHandle->set_level_2_mem(level_2_mem);

      gsHandle->set_num_values_in_l1(num_values_in_l1);
      gsHandle->set_num_values_in_l2(num_values_in_l2);
      gsHandle->set_num_big_rows(num_big_rows);
    }

    gsHandle->set_color_xadj(h_color_xadj);
    gsHandle->set_color_adj(color_adj);
    gsHandle->set_num_colors(numColors);
    gsHandle->set_new_xadj(permuted_xadj);
    gsHandle->set_new_adj(permuted_adj);
    gsHandle->set_old_to_new_map(old_to_new_map);
    gsHandle->set_call_symbolic(true);
    my_exec_space.fence();
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    std::cout << "ALLOC:" << timer.seconds() << std::endl;
#endif
  }

  struct create_permuted_xadj {
    nnz_lno_persistent_work_view_t color_adj;
    const_lno_row_view_t oldxadj;
    row_lno_persistent_work_view_t newxadj;
    nnz_lno_persistent_work_view_t old_to_new_index;
    create_permuted_xadj(nnz_lno_persistent_work_view_t color_adj_, const_lno_row_view_t oldxadj_,
                         row_lno_persistent_work_view_t newxadj_, nnz_lno_persistent_work_view_t old_to_new_index_)
        : color_adj(color_adj_), oldxadj(oldxadj_), newxadj(newxadj_), old_to_new_index(old_to_new_index_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t& i) const {
      nnz_lno_t index         = color_adj(i);
      newxadj(i + 1)          = oldxadj[index + 1] - oldxadj[index];
      old_to_new_index[index] = i;
    }
  };

  struct fill_matrix_symbolic {
    nnz_lno_t num_rows;
    nnz_lno_persistent_work_view_t color_adj;
    const_lno_row_view_t oldxadj;
    const_lno_nnz_view_t oldadj;
    // value_array_type oldadjvals;
    row_lno_persistent_work_view_t newxadj;
    nnz_lno_persistent_work_view_t newadj;
    // value_persistent_work_array_type newadjvals;
    nnz_lno_persistent_work_view_t old_to_new_index;
    fill_matrix_symbolic(nnz_lno_t num_rows_, nnz_lno_persistent_work_view_t color_adj_, const_lno_row_view_t oldxadj_,
                         const_lno_nnz_view_t oldadj_,
                         // value_array_type oldadjvals_,
                         row_lno_persistent_work_view_t newxadj_, nnz_lno_persistent_work_view_t newadj_,
                         // value_persistent_work_array_type newadjvals_,
                         nnz_lno_persistent_work_view_t old_to_new_index_)
        : num_rows(num_rows_),
          color_adj(color_adj_),
          oldxadj(oldxadj_),
          oldadj(oldadj_),  // oldadjvals(oldadjvals_),
          newxadj(newxadj_),
          newadj(newadj_),  // newadjvals(newadjvals_),
          old_to_new_index(old_to_new_index_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t& i) const {
      nnz_lno_t index      = color_adj(i);
      size_type xadj_begin = newxadj(i);

      size_type old_xadj_end = oldxadj[index + 1];
      for (size_type j = oldxadj[index]; j < old_xadj_end; ++j) {
        nnz_lno_t neighbor = oldadj[j];
        if (neighbor < num_rows) neighbor = old_to_new_index[neighbor];
        newadj[xadj_begin++] = neighbor;
        // newadjvals[xadj_begin++] = oldadjvals[j];
      }
    }
  };

  struct fill_matrix_numeric {
    nnz_lno_persistent_work_view_t color_adj;
    const_lno_row_view_t oldxadj;
    const_scalar_nnz_view_t oldadjvals;
    row_lno_persistent_work_view_t newxadj;
    scalar_persistent_work_view_t newadjvals;

    nnz_lno_t num_total_rows;
    nnz_lno_t rows_per_team;
    nnz_lno_t block_matrix_size;
    fill_matrix_numeric(nnz_lno_persistent_work_view_t color_adj_, const_lno_row_view_t oldxadj_,
                        const_scalar_nnz_view_t oldadjvals_, row_lno_persistent_work_view_t newxadj_,
                        scalar_persistent_work_view_t newadjvals_, nnz_lno_t num_total_rows_, nnz_lno_t rows_per_team_,
                        nnz_lno_t block_matrix_size_)
        : color_adj(color_adj_),
          oldxadj(oldxadj_),
          oldadjvals(oldadjvals_),
          newxadj(newxadj_),
          newadjvals(newadjvals_),
          num_total_rows(num_total_rows_),
          rows_per_team(rows_per_team_),
          block_matrix_size(block_matrix_size_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t& i) const {
      nnz_lno_t index        = color_adj(i);
      size_type xadj_begin   = newxadj(i) * block_matrix_size;
      size_type old_xadj_end = oldxadj[index + 1] * block_matrix_size;

      for (size_type j = oldxadj[index] * block_matrix_size; j < old_xadj_end; ++j) {
        newadjvals[xadj_begin++] = oldadjvals[j];
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t& team) const {
      const nnz_lno_t i_begin = team.league_rank() * rows_per_team;
      const nnz_lno_t i_end   = i_begin + rows_per_team <= num_total_rows ? i_begin + rows_per_team : num_total_rows;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, i_begin, i_end), [&](const nnz_lno_t& i) {
        nnz_lno_t index      = color_adj(i);
        size_type xadj_begin = newxadj(i) * block_matrix_size;

        size_type old_xadj_begin = oldxadj[index] * block_matrix_size;
        size_type old_xadj_end   = oldxadj[index + 1] * block_matrix_size;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, old_xadj_end - old_xadj_begin),
                             [&](const nnz_lno_t& j) { newadjvals[xadj_begin + j] = oldadjvals[old_xadj_begin + j]; });
      });
    }
  };

  struct Get_Matrix_Diagonals {
    row_lno_persistent_work_view_t _xadj;
    nnz_lno_persistent_work_view_t _adj;      // CSR storage of the graph.
    scalar_persistent_work_view_t _adj_vals;  // CSR storage of the graph.
    scalar_persistent_work_view_t _diagonals;

    nnz_lno_t num_total_rows;
    nnz_lno_t rows_per_team;
    nnz_lno_t block_size;
    nnz_lno_t block_matrix_size;

    nnz_scalar_t one;

    Get_Matrix_Diagonals(row_lno_persistent_work_view_t xadj_, nnz_lno_persistent_work_view_t adj_,
                         scalar_persistent_work_view_t adj_vals_, scalar_persistent_work_view_t diagonals_,
                         nnz_lno_t num_total_rows_, nnz_lno_t rows_per_team_, nnz_lno_t block_size_,
                         nnz_lno_t block_matrix_size_)
        : _xadj(xadj_),
          _adj(adj_),
          _adj_vals(adj_vals_),
          _diagonals(diagonals_),
          num_total_rows(num_total_rows_),
          rows_per_team(rows_per_team_),
          block_size(block_size_),
          block_matrix_size(block_matrix_size_),
          one(Kokkos::ArithTraits<nnz_scalar_t>::one()) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const nnz_lno_t& row_id) const {
      RowIndex row(block_size, _xadj[row_id], _xadj[row_id + 1]);
      for (nnz_lno_t col_ind = 0; col_ind < row.size(); ++col_ind) {
        const nnz_lno_t column_id = _adj[row.begin() + col_ind];
        if (column_id == row_id) {
          size_type val_index = row.block(col_ind);
          for (nnz_lno_t r = 0; r < block_size; ++r) {
            nnz_scalar_t val                    = _adj_vals[val_index];
            _diagonals[row_id * block_size + r] = one / val;
            val_index += row.block_stride() + 1;
          }
          break;
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const team_member_t& team) const {
      const nnz_lno_t i_begin = team.league_rank() * rows_per_team;
      const nnz_lno_t i_end   = i_begin + rows_per_team <= num_total_rows ? i_begin + rows_per_team : num_total_rows;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, i_begin, i_end), [&](const nnz_lno_t& row_id) {
        RowIndex row(block_size, _xadj[row_id], _xadj[row_id + 1]);

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, row.size()), [&](const nnz_lno_t& col_ind) {
          size_type val_index = col_ind + row.begin();
          nnz_lno_t column_id = _adj[val_index];
          if (column_id == row_id) {
            size_type _val_index = row.block(col_ind);
            for (nnz_lno_t r = 0; r < block_size; ++r) {
              nnz_scalar_t val                    = _adj_vals[_val_index];
              _diagonals[row_id * block_size + r] = one / val;
              _val_index += row.block_stride() + 1;
            }
          }
        });
      });
    }
  };

  void initialize_numeric() {
    auto gsHandle = this->get_gs_handle();
    if (gsHandle->is_symbolic_called() == false) {
      this->initialize_symbolic();
    }

    // Check settings
    if (gsHandle->get_block_size() > 1 && format != KokkosSparse::SparseMatrixFormat::BSR)
      throw std::runtime_error(
          "PointGaussSeidel block size > 1 but format is not "
          "KokkosSparse::SparseMatrixFormat::BSR.\n");
      // else
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    Kokkos::Timer timer;
#endif
    {
      const_lno_row_view_t xadj        = this->row_map;
      const_lno_nnz_view_t adj         = this->entries;
      const_scalar_nnz_view_t adj_vals = this->values;
      MyExecSpace my_exec_space        = gsHandle->get_execution_space();

      size_type nnz = adj_vals.extent(0);

      row_lno_persistent_work_view_t newxadj_       = gsHandle->get_new_xadj();
      nnz_lno_persistent_work_view_t newadj_        = gsHandle->get_new_adj();
      nnz_lno_persistent_work_view_t old_to_new_map = gsHandle->get_old_to_new_map();

      nnz_lno_persistent_work_view_t color_adj = gsHandle->get_color_adj();
      scalar_persistent_work_view_t permuted_adj_vals(
          Kokkos::view_alloc(my_exec_space, Kokkos::WithoutInitializing, "newvals_"), nnz);

      int suggested_vector_size = this->handle->get_suggested_vector_size(num_rows, nnz);
      int suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);
      nnz_lno_t rows_per_team =
          this->handle->get_team_work_size(suggested_team_size, my_exec_space.concurrency(), num_rows);

      nnz_lno_t block_size        = gsHandle->get_block_size();
      nnz_lno_t block_matrix_size = block_size * block_size;

      // MD NOTE: 03/27/2018: below fill matrix operations will work fine with
      // block size 1. If the block size is more than 1, below code assumes that
      // the rows are sorted similar to point crs. for example given a block crs
      // with 3 blocks in a column a,b,c where each of them is 3x3 matrix as
      // below.
      // a11 a12 a13   b11 b12 b13    c11 c12 c13
      // a21 a22 a23   b21 b22 b23    c21 c22 c23
      // a31 a32 a33   b31 b32 b33    c31 c32 c33
      // this copy assumes the storage in the following order
      // a11 a12 a13   b11 b12 b13    c11 c12 c13 a21 a22 a23   b21 b22 b23 c21
      // c22 c23 a31 a32 a33   b31 b32 b33    c31 c32 c33 this is the order that
      // is used in the rest of the algorithm.
      // !!!!!!!!!!!!if the input has different format than
      // this!!!!!!!!!!!!!!!!!! change fill_matrix_numeric so that they store
      // the internal matrix as above. the rest will wok fine.

      if (KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>()) {
        Kokkos::parallel_for("KokkosSparse::GaussSeidel::Team_fill_matrix_numeric",
                             team_policy_t(my_exec_space, (num_rows + rows_per_team - 1) / rows_per_team,
                                           suggested_team_size, suggested_vector_size),
                             fill_matrix_numeric(color_adj, xadj,
                                                 // adj,
                                                 adj_vals, newxadj_,
                                                 // newadj_,
                                                 permuted_adj_vals,
                                                 //,old_to_new_map
                                                 this->num_rows, rows_per_team, block_matrix_size));
      } else {
        Kokkos::parallel_for("KokkosSparse::GaussSeidel::fill_matrix_numeric",
                             range_policy_t(my_exec_space, 0, num_rows),
                             fill_matrix_numeric(color_adj, xadj,
                                                 // adj,
                                                 adj_vals, newxadj_,
                                                 // newadj_,
                                                 permuted_adj_vals,
                                                 //,old_to_new_map
                                                 this->num_rows, rows_per_team, block_matrix_size));
      }
      gsHandle->set_new_adj_val(permuted_adj_vals);

      scalar_persistent_work_view_t permuted_inverse_diagonal(
          Kokkos::view_alloc(my_exec_space, Kokkos::WithoutInitializing, "permuted_inverse_diagonal"),
          num_rows * block_size);
      if (!have_diagonal_given) {
        Get_Matrix_Diagonals gmd(newxadj_, newadj_, permuted_adj_vals, permuted_inverse_diagonal, this->num_rows,
                                 rows_per_team, block_size, block_matrix_size);

        if (KokkosKernels::Impl::kk_is_gpu_exec_space<MyExecSpace>() || block_size > 1) {
          Kokkos::parallel_for("KokkosSparse::GaussSeidel::team_get_matrix_diagonals",
                               team_policy_t(my_exec_space, (num_rows + rows_per_team - 1) / rows_per_team,
                                             suggested_team_size, suggested_vector_size),
                               gmd);
        } else {
          Kokkos::parallel_for("KokkosSparse::GaussSeidel::get_matrix_diagonals",
                               range_policy_t(my_exec_space, 0, num_rows), gmd);
        }

      } else {
        if (block_size > 1)
          KokkosKernels::Impl::permute_block_vector<const_scalar_nnz_view_t, scalar_persistent_work_view_t,
                                                    nnz_lno_persistent_work_view_t, MyExecSpace>(
              my_exec_space, num_rows, block_size, old_to_new_map, given_inverse_diagonal, permuted_inverse_diagonal);
        else
          KokkosKernels::Impl::permute_vector<const_scalar_nnz_view_t, scalar_persistent_work_view_t,
                                              nnz_lno_persistent_work_view_t, MyExecSpace>(
              my_exec_space, num_rows, old_to_new_map, given_inverse_diagonal, permuted_inverse_diagonal);
      }
      gsHandle->set_permuted_inverse_diagonal(permuted_inverse_diagonal);
      gsHandle->set_call_numeric(true);
    }
#ifdef KOKKOSSPARSE_IMPL_TIME_REVERSE
    my_exec_space.fence();
    std::cout << "NUMERIC:" << timer.seconds() << std::endl;
#endif
  }

  // Functor to update unknown entries corresponding to long rows (in the
  // permuted x)
  template <typename x_value_array_type, typename y_value_array_type>
  struct LongRowUpdateFunctor {
    LongRowUpdateFunctor(const x_value_array_type& permuted_x_, const y_value_array_type& permuted_y_,
                         const scalar_persistent_work_view_t& long_row_x_,
                         const scalar_persistent_work_view_t& permuted_inverse_diagonal_, nnz_scalar_t omega_,
                         nnz_lno_t long_row_begin_)
        : permuted_x(permuted_x_),
          permuted_y(permuted_y_),
          long_row_x(long_row_x_),
          permuted_inverse_diagonal(permuted_inverse_diagonal_),
          omega(omega_),
          long_row_begin(long_row_begin_) {}

    KOKKOS_INLINE_FUNCTION void operator()(nnz_lno_t i) const {
      permuted_x(i) += omega * permuted_inverse_diagonal(i) * (permuted_y(i) - long_row_x(i - long_row_begin));
    }

    x_value_array_type permuted_x;
    y_value_array_type permuted_y;
    scalar_persistent_work_view_t long_row_x;
    scalar_persistent_work_view_t permuted_inverse_diagonal;
    nnz_scalar_t omega;
    nnz_lno_t long_row_begin;
  };

  template <typename x_value_array_type, typename y_value_array_type>
  void block_apply(x_value_array_type x_lhs_output_vec, y_value_array_type y_rhs_input_vec,
                   bool init_zero_x_vector = false, int numIter = 1,
                   nnz_scalar_t omega = Kokkos::ArithTraits<nnz_scalar_t>::one(), bool apply_forward = true,
                   bool apply_backward = true, bool update_y_vector = true) {
    auto gsHandle = this->get_gs_handle();
    if (gsHandle->is_numeric_called() == false) {
      this->initialize_numeric();
    }

    nnz_lno_t block_size = gsHandle->get_block_size();
    // nnz_lno_t block_matrix_size = block_size  * block_size ;

    auto Permuted_Xvector = gsHandle->get_permuted_x_vector();
    auto Permuted_Yvector = gsHandle->get_permuted_y_vector();

    row_lno_persistent_work_view_t newxadj                  = gsHandle->get_new_xadj();
    nnz_lno_persistent_work_view_t newadj                   = gsHandle->get_new_adj();
    scalar_persistent_work_view_t newadj_vals               = gsHandle->get_new_adj_val();
    nnz_lno_persistent_work_view_t old_to_new_map           = gsHandle->get_old_to_new_map();
    nnz_lno_persistent_work_view_t color_adj                = gsHandle->get_color_adj();
    scalar_persistent_work_view_t permuted_inverse_diagonal = gsHandle->get_permuted_inverse_diagonal();

    color_t numColors  = gsHandle->get_num_colors();
    auto my_exec_space = gsHandle->get_execution_space();

    if (update_y_vector) {
      KokkosKernels::Impl::permute_block_vector<y_value_array_type, scalar_persistent_work_view2d_t,
                                                nnz_lno_persistent_work_view_t, MyExecSpace>(
          my_exec_space, num_rows, block_size, old_to_new_map, y_rhs_input_vec, Permuted_Yvector);
    }
    if (init_zero_x_vector) {
      KokkosKernels::Impl::zero_vector<MyExecSpace, scalar_persistent_work_view2d_t>(
          my_exec_space, num_cols * block_size, Permuted_Xvector);
    } else {
      KokkosKernels::Impl::permute_block_vector<x_value_array_type, scalar_persistent_work_view2d_t,
                                                nnz_lno_persistent_work_view_t, MyExecSpace>(
          my_exec_space, num_cols, block_size, old_to_new_map, x_lhs_output_vec, Permuted_Xvector);
    }

#if KOKKOSSPARSE_IMPL_PRINTDEBUG
    std::cout << "Y:";
    KokkosKernels::Impl::print_1Dview(Permuted_Yvector);
    std::cout << "Original Y:";
    KokkosKernels::Impl::print_1Dview(y_rhs_input_vec);

    std::cout << "X:";
    KokkosKernels::Impl::print_1Dview(Permuted_Xvector);

    std::cout << "permuted_xadj:";
    KokkosKernels::Impl::print_1Dview(newxadj);
    std::cout << "permuted_adj:";
    KokkosKernels::Impl::print_1Dview(newadj);
    std::cout << "permuted_adj_vals:";
    KokkosKernels::Impl::print_1Dview(newadj_vals);
    std::cout << "permuted_diagonals:";
    KokkosKernels::Impl::print_1Dview(permuted_inverse_diagonal);
#endif
    nnz_lno_persistent_work_host_view_t h_color_xadj = gsHandle->get_color_xadj();

    nnz_lno_t brows = newxadj.extent(0) - 1;
    size_type bnnz  = newadj_vals.extent(0);

    int suggested_vector_size = this->handle->get_suggested_vector_size(brows, bnnz);
    int suggested_team_size   = this->handle->get_suggested_team_size(suggested_vector_size);
    nnz_lno_t team_row_chunk_size =
        this->handle->get_team_work_size(suggested_team_size, my_exec_space.concurrency(), brows);

    // size_t shmem_size_to_use = this->handle->get_shmem_size();
    size_t l1_shmem_size       = gsHandle->get_level_1_mem();
    nnz_lno_t num_values_in_l1 = gsHandle->get_num_values_in_l1();

    size_t level_2_mem         = gsHandle->get_level_2_mem();
    nnz_lno_t num_values_in_l2 = gsHandle->get_num_values_in_l2();
    nnz_lno_t num_chunks       = gsHandle->get_num_big_rows();

    pool_memory_space m_space(num_chunks, level_2_mem / sizeof(nnz_scalar_t), 0,
                              KokkosKernels::Impl::ManyThread2OneChunk, false);

#if KOKKOSSPARSE_IMPL_PRINTDEBUG
    std::cout << "l1_shmem_size:" << l1_shmem_size << " num_values_in_l1:" << num_values_in_l1
              << " level_2_mem:" << level_2_mem << " num_values_in_l2:" << num_values_in_l2
              << " num_chunks:" << num_chunks << std::endl;
#endif

    Team_PSGS gs(newxadj, newadj, newadj_vals, Permuted_Xvector, Permuted_Yvector, 0, 0, permuted_inverse_diagonal,
                 m_space, num_values_in_l1, num_values_in_l2, omega, block_size, team_row_chunk_size, l1_shmem_size,
                 suggested_team_size, suggested_vector_size);

    this->IterativePSGS(gs, numColors, h_color_xadj, numIter, apply_forward, apply_backward);

    KokkosKernels::Impl::permute_block_vector<scalar_persistent_work_view2d_t, x_value_array_type,
                                              nnz_lno_persistent_work_view_t, MyExecSpace>(
        my_exec_space, num_cols, block_size, color_adj, Permuted_Xvector, x_lhs_output_vec);
#if KOKKOSSPARSE_IMPL_PRINTDEBUG
    std::cout << "After X:";
    KokkosKernels::Impl::print_1Dview(Permuted_Xvector);
    std::cout << "Result X:";
    KokkosKernels::Impl::print_1Dview(x_lhs_output_vec);
    std::cout << "Y:";
    KokkosKernels::Impl::print_1Dview(Permuted_Yvector);
#endif
  }

  template <typename x_value_array_type, typename y_value_array_type>
  void point_apply(x_value_array_type x_lhs_output_vec, y_value_array_type y_rhs_input_vec,
                   bool init_zero_x_vector = false, int numIter = 1,
                   nnz_scalar_t omega = Kokkos::ArithTraits<nnz_scalar_t>::one(), bool apply_forward = true,
                   bool apply_backward = true, bool update_y_vector = true) {
    auto gsHandle      = get_gs_handle();
    auto my_exec_space = gsHandle->get_execution_space();

    auto Permuted_Xvector = gsHandle->get_permuted_x_vector();
    auto Permuted_Yvector = gsHandle->get_permuted_y_vector();

    row_lno_persistent_work_view_t newxadj                  = gsHandle->get_new_xadj();
    nnz_lno_persistent_work_view_t newadj                   = gsHandle->get_new_adj();
    scalar_persistent_work_view_t newadj_vals               = gsHandle->get_new_adj_val();
    nnz_lno_persistent_work_view_t old_to_new_map           = gsHandle->get_old_to_new_map();
    nnz_lno_persistent_work_view_t color_adj                = gsHandle->get_color_adj();
    scalar_persistent_work_view_t permuted_inverse_diagonal = gsHandle->get_permuted_inverse_diagonal();

    color_t numColors = gsHandle->get_num_colors();

    if (update_y_vector) {
      KokkosKernels::Impl::permute_vector<y_value_array_type, scalar_persistent_work_view2d_t,
                                          nnz_lno_persistent_work_view_t, MyExecSpace>(
          my_exec_space, num_rows, old_to_new_map, y_rhs_input_vec, Permuted_Yvector);
    }
    if (init_zero_x_vector) {
      KokkosKernels::Impl::zero_vector<MyExecSpace, scalar_persistent_work_view2d_t>(my_exec_space, num_cols,
                                                                                     Permuted_Xvector);
    } else {
      KokkosKernels::Impl::permute_vector<x_value_array_type, scalar_persistent_work_view2d_t,
                                          nnz_lno_persistent_work_view_t, MyExecSpace>(
          my_exec_space, num_cols, old_to_new_map, x_lhs_output_vec, Permuted_Xvector);
    }

#if KOKKOSSPARSE_IMPL_PRINTDEBUG
    std::cout << "--point Before X:";
    KokkosKernels::Impl::print_1Dview(Permuted_Xvector, true);
    std::cout << "--point Before Y:";
    KokkosKernels::Impl::print_1Dview(Permuted_Yvector, true);
#endif

    nnz_lno_persistent_work_host_view_t h_color_xadj = gsHandle->get_color_xadj();
    if (gsHandle->get_algorithm_type() == GS_PERMUTED) {
      PSGS gs(newxadj, newadj, newadj_vals, Permuted_Xvector, Permuted_Yvector, color_adj, omega,
              permuted_inverse_diagonal);
      this->IterativePSGS(gs, numColors, h_color_xadj, numIter, apply_forward, apply_backward);
    } else {
      pool_memory_space m_space(0, 0, 0, KokkosKernels::Impl::ManyThread2OneChunk, false);

      Team_PSGS gs(newxadj, newadj, newadj_vals, Permuted_Xvector, Permuted_Yvector, 0, 0, permuted_inverse_diagonal,
                   m_space, 0, 0, omega);

      this->IterativePSGS(gs, numColors, h_color_xadj, numIter, apply_forward, apply_backward);
    }

    KokkosKernels::Impl::permute_vector<scalar_persistent_work_view2d_t, x_value_array_type,
                                        nnz_lno_persistent_work_view_t, MyExecSpace>(
        my_exec_space, num_cols, color_adj, Permuted_Xvector, x_lhs_output_vec);
#if KOKKOSSPARSE_IMPL_PRINTDEBUG
    Kokkos::fence();
    std::cout << "--point After X:";
    KokkosKernels::Impl::print_1Dview(Permuted_Xvector);
    std::cout << "--point Result X:";
    KokkosKernels::Impl::print_1Dview(x_lhs_output_vec);
#endif
  }

  template <typename x_value_array_type, typename y_value_array_type>
  void apply(x_value_array_type x_lhs_output_vec, y_value_array_type y_rhs_input_vec, bool init_zero_x_vector = false,
             int numIter = 1, nnz_scalar_t omega = Kokkos::ArithTraits<nnz_scalar_t>::one(), bool apply_forward = true,
             bool apply_backward = true, bool update_y_vector = true) {
    auto gsHandle = get_gs_handle();
    if (gsHandle->is_numeric_called() == false) {
      this->initialize_numeric();
    }

    // Check settings
    if (gsHandle->get_block_size() > 1 && format != KokkosSparse::SparseMatrixFormat::BSR)
      throw std::runtime_error(
          "PointGaussSeidel block size > 1 but format is not "
          "KokkosSparse::SparseMatrixFormat::BSR.\n");

    // make sure x and y have been allocated with the correct dimensions
    nnz_lno_t block_size = gsHandle->get_block_size();
    gsHandle->allocate_x_y_vectors(this->num_rows * block_size, this->num_cols * block_size,
                                   x_lhs_output_vec.extent(1));
    if (block_size == 1) {
      this->point_apply(x_lhs_output_vec, y_rhs_input_vec, init_zero_x_vector, numIter, omega, apply_forward,
                        apply_backward, update_y_vector);
    } else {
      this->block_apply(x_lhs_output_vec, y_rhs_input_vec, init_zero_x_vector, numIter, omega, apply_forward,
                        apply_backward, update_y_vector);
    }
  }

  void IterativePSGS(Team_PSGS& gs, color_t numColors, nnz_lno_persistent_work_host_view_t h_color_xadj,
                     int num_iteration, bool apply_forward, bool apply_backward) {
    auto gsHandle             = this->get_gs_handle();
    MyExecSpace my_exec_space = gsHandle->get_execution_space();
    nnz_lno_persistent_work_host_view_t long_rows_per_color;
    nnz_lno_persistent_work_host_view_t max_row_length_per_color;
    scalar_persistent_work_view_t long_row_x;
    bool haveLongRows   = false;
    int longRowTeamSize = 1;
    if (gsHandle->get_long_row_threshold() > 0) {
      long_rows_per_color      = gsHandle->get_long_rows_per_color();
      max_row_length_per_color = gsHandle->get_max_row_length_per_color();
      long_row_x               = gsHandle->get_long_row_x();
      haveLongRows             = true;
      longrow_apply_team_policy_t tempPolicy(my_exec_space, 1, 1);
      longRowTeamSize = tempPolicy.team_size_recommended(gs, Kokkos::ParallelForTag());
    }

    for (int iter = 0; iter < num_iteration; ++iter) {
      nnz_lno_t suggested_team_size = gs.suggested_team_size;
      nnz_lno_t team_row_chunk_size = gs.team_work_size;
      int vector_size               = gs.vector_size;
      nnz_lno_t block_size          = gsHandle->get_block_size();

      for (int doingBackward = 0; doingBackward < 2; doingBackward++) {
        const char* labelRegular  = doingBackward ? "KokkosSparse::GaussSeidel::Team_PSGS::backward"
                                                  : "KokkosSparse::GaussSeidel::Team_PSGS::forward";
        const char* labelBlock    = doingBackward ? "KokkosSparse::GaussSeidel::BLOCK_Team_PSGS::backward"
                                                  : "KokkosSparse::GaussSeidel::BLOCK_Team_PSGS::forward";
        const char* labelBigBlock = doingBackward ? "KokkosSparse::GaussSeidel::BIGBLOCK_Team_PSGS::backward"
                                                  : "KokkosSparse::GaussSeidel::BIGBLOCK_Team_PSGS::forward";
        const char* labelLong     = doingBackward ? "KokkosSparse::GaussSeidel::Team_PSGS::backwardLongRows"
                                                  : "KokkosSparse::GaussSeidel::Team_PSGS::forwardLongRows";

        if (!doingBackward && !apply_forward) continue;
        if (doingBackward && !apply_backward) continue;
        gs.is_backward = doingBackward;

        for (color_t colorIter = 0; colorIter < numColors; ++colorIter) {
          // i is just the color set now being processed
          color_t i                   = doingBackward ? (numColors - colorIter - 1) : colorIter;
          nnz_lno_t color_index_begin = h_color_xadj(i);
          nnz_lno_t color_index_end   = h_color_xadj(i + 1);
          nnz_lno_t numLongRows       = haveLongRows ? long_rows_per_color(i) : 0;
          nnz_lno_t numRegularRows    = color_index_end - color_index_begin - numLongRows;

          gs._color_set_begin = color_index_begin;
          gs._color_set_end   = color_index_end - numLongRows;

          if (numRegularRows) {
            if (block_size == 1) {
              Kokkos::parallel_for(
                  labelRegular,
                  Kokkos::Experimental::require(
                      team_policy_t(my_exec_space, (numRegularRows + team_row_chunk_size - 1) / team_row_chunk_size,
                                    suggested_team_size, vector_size),
                      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  gs);
            } else if (gs.num_max_vals_in_l2 == 0) {
              Kokkos::parallel_for(
                  labelBlock,
                  Kokkos::Experimental::require(
                      block_apply_team_policy_t(my_exec_space,
                                                (numRegularRows + team_row_chunk_size - 1) / team_row_chunk_size,
                                                suggested_team_size, vector_size),
                      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  gs);
            } else {
              Kokkos::parallel_for(
                  labelBigBlock,
                  Kokkos::Experimental::require(
                      bigblock_apply_team_policy_t(my_exec_space,
                                                   (numRegularRows + team_row_chunk_size - 1) / team_row_chunk_size,
                                                   suggested_team_size, vector_size),
                      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  gs);
            }
          }
          if (numLongRows) {
            gs._color_set_begin     = color_index_end - numLongRows;
            gs._color_set_end       = color_index_end;
            gs._long_row_x          = long_row_x;
            nnz_lno_t max_par       = max_row_length_per_color(i);
            nnz_lno_t teams_per_row = ((max_par + 3) / 4 + longRowTeamSize - 1) / longRowTeamSize;
            gs._long_row_par        = teams_per_row;
            for (nnz_lno_t long_row_col = 0; long_row_col < gs._Xvector.extent_int(1); long_row_col++) {
              auto Xcol        = Kokkos::subview(gs._Xvector, Kokkos::ALL(), long_row_col);
              auto Ycol        = Kokkos::subview(gs._Yvector, Kokkos::ALL(), long_row_col);
              gs._long_row_col = long_row_col;
              Kokkos::deep_copy(my_exec_space, long_row_x, nnz_scalar_t());
              Kokkos::parallel_for(
                  labelLong,
                  Kokkos::Experimental::require(
                      longrow_apply_team_policy_t(my_exec_space, numLongRows * teams_per_row, longRowTeamSize),
                      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  gs);
              Kokkos::parallel_for(
                  "KokkosSparse::GaussSeidel::LongRows::x_update",
                  Kokkos::Experimental::require(
                      range_policy_t(my_exec_space, color_index_end - numLongRows, color_index_end),
                      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  LongRowUpdateFunctor<decltype(Xcol), decltype(Ycol)>(
                      Xcol, Ycol, long_row_x, gs._permuted_inverse_diagonal, gs.omega, color_index_end - numLongRows));
            }
          }
        }
      }
    }
  }

  void IterativePSGS(PSGS& gs, color_t numColors, nnz_lno_persistent_work_host_view_t h_color_xadj, int num_iteration,
                     bool apply_forward, bool apply_backward) {
    auto gsHandle             = this->get_gs_handle();
    MyExecSpace my_exec_space = gsHandle->get_execution_space();
    nnz_lno_persistent_work_host_view_t long_rows_per_color;
    nnz_lno_persistent_work_host_view_t max_row_length_per_color;
    scalar_persistent_work_view_t long_row_x;
    bool haveLongRows = false;
    if (gsHandle->get_long_row_threshold() > 0) {
      long_rows_per_color      = gsHandle->get_long_rows_per_color();
      max_row_length_per_color = gsHandle->get_max_row_length_per_color();
      long_row_x               = gsHandle->get_long_row_x();
      gs._long_row_x           = long_row_x;
      haveLongRows             = true;
    }

    for (int iter = 0; iter < num_iteration; ++iter) {
      for (int doingBackward = 0; doingBackward < 2; doingBackward++) {
        if (!doingBackward && !apply_forward) continue;
        if (doingBackward && !apply_backward) continue;

        for (color_t colorIter = 0; colorIter < numColors; ++colorIter) {
          // i is just the color set now being processed
          color_t i = doingBackward ? (numColors - colorIter - 1) : colorIter;
          const char* labelShort =
              doingBackward ? "KokkosSparse::GaussSeidel::PSGS::backward" : "KokkosSparse::GaussSeidel::PSGS::forward";
          const char* labelLong       = doingBackward ? "KokkosSparse::GaussSeidel::PSGS::backwardLongRows"
                                                      : "KokkosSparse::GaussSeidel::PSGS::forwardLongRows";
          nnz_lno_t color_index_begin = h_color_xadj(i);
          nnz_lno_t color_index_end   = h_color_xadj(i + 1);
          nnz_lno_t numLongRows       = haveLongRows ? long_rows_per_color(i) : 0;
          nnz_lno_t numRegularRows    = color_index_end - color_index_begin - numLongRows;
          if (numRegularRows) {
            Kokkos::parallel_for(labelShort,
                                 Kokkos::Experimental::require(
                                     range_policy_t(my_exec_space, color_index_begin, color_index_end - numLongRows),
                                     Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                                 gs);
          }
          if (numLongRows) {
            gs._color_set_begin   = color_index_end - numLongRows;
            nnz_lno_t max_par     = max_row_length_per_color(i);
            nnz_lno_t par_per_row = (max_par + 1023) / 1024;
            gs._long_row_par      = par_per_row;
            for (nnz_lno_t long_row_col = 0; long_row_col < gs._Xvector.extent_int(1); long_row_col++) {
              auto Xcol        = Kokkos::subview(gs._Xvector, Kokkos::ALL(), long_row_col);
              auto Ycol        = Kokkos::subview(gs._Yvector, Kokkos::ALL(), long_row_col);
              gs._long_row_col = long_row_col;
              Kokkos::deep_copy(my_exec_space, long_row_x, nnz_scalar_t());
              Kokkos::parallel_for(
                  labelLong,
                  Kokkos::Experimental::require(
                      Kokkos::RangePolicy<MyExecSpace, LongRowTag>(my_exec_space, 0, numLongRows * par_per_row),
                      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  gs);
              Kokkos::parallel_for(
                  "KokkosSparse::GaussSeidel::LongRows::x_update",
                  Kokkos::Experimental::require(
                      range_policy_t(my_exec_space, color_index_end - numLongRows, color_index_end),
                      Kokkos::Experimental::WorkItemProperty::HintLightWeight),
                  LongRowUpdateFunctor<decltype(Xcol), decltype(Ycol)>(
                      Xcol, Ycol, long_row_x, gs._permuted_inverse_diagonal, gs.omega, color_index_end - numLongRows));
            }
          }
        }
      }
    }
  }
};
}  // namespace Impl
}  // namespace KokkosSparse
#endif
