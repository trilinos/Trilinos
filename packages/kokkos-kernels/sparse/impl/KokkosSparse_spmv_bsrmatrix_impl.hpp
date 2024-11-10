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

#ifndef KOKKOSSPARSE_IMPL_SPMV_BSRMATRIX_IMPL_HPP_
#define KOKKOSSPARSE_IMPL_SPMV_BSRMATRIX_IMPL_HPP_

#include "KokkosKernels_Error.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"

#if defined(KOKKOS_ENABLE_CUDA) && (defined(KOKKOS_ARCH_VOLTA) || defined(KOKKOS_ARCH_AMPERE))

#include <type_traits>
#include <mma.h>

namespace KokkosSparse {
namespace Impl {

struct BsrMatrixSpMVTensorCoreFunctorParams {
  int teamsPerBlockM;
  int teamsPerBlockN;
  int leagueDim_x;
  int leagueDim_y;
};

/*! \brief Can the tensor core impl be used in ExecutionSpace to operate on
    AMatrix, XMatrix, and YMatrix?
*/
template <typename ExecutionSpace, typename AMatrix, typename XMatrix, typename YMatrix>
class TensorCoresAvailable {
#if defined(KOKKOS_ENABLE_CUDA)
  using AScalar = typename AMatrix::non_const_value_type;
  using YScalar = typename YMatrix::non_const_value_type;
  using XScalar = typename XMatrix::non_const_value_type;

  using a_mem_space = typename AMatrix::memory_space;
  using x_mem_space = typename XMatrix::memory_space;
  using y_mem_space = typename YMatrix::memory_space;

  template <typename T>
  constexpr static bool is_scalar() {
    return std::is_scalar_v<T> || std::is_same_v<std::remove_cv_t<T>, Kokkos::Experimental::half_t>;
  }

 public:
  constexpr static inline bool value = Kokkos::SpaceAccessibility<ExecutionSpace, a_mem_space>::accessible &&
                                       Kokkos::SpaceAccessibility<ExecutionSpace, x_mem_space>::accessible &&
                                       Kokkos::SpaceAccessibility<ExecutionSpace, y_mem_space>::accessible &&
                                       is_scalar<AScalar>() && is_scalar<XScalar>() && is_scalar<YScalar>() &&
                                       std::is_same_v<ExecutionSpace, Kokkos::Cuda>;
#else
 public:
  constexpr static inline bool value = false;
#endif
};

/// \brief Functor for the BsrMatrix SpMV multivector implementation utilizing
/// tensor cores.
///
/// \tparam AMatrix The type of the A matrix (a BsrMatrix)
/// \tparam AFragScalar The type of the CUDA wmma fragment that will be loaded
/// from the A matrix. The scalar type of the wmma fragment may be different
/// that that of the A matrix. \tparam FRAG_M (with FRAG_N and FRAG_K), the
/// m-n-k size of the CUDA wmma fragment type. \tparam LEAGUE_DIM_X (with
/// TEAMS_PER_BLOCK_M and TEAMS_PER_BLOCK_N) if non-zero, statically-known
/// launch parameters to reduce the cost of divmod operations on the GPU. If 0,
/// provided runtime values will be used instead.
template <typename execution_space, typename AMatrix,
          typename AFragScalar,  // input matrix type and fragment scalar type
          typename XMatrix, typename XFragScalar, typename YMatrix, typename YFragScalar, unsigned FRAG_M,
          unsigned FRAG_N,
          unsigned FRAG_K,  // fragment sizes
          unsigned LEAGUE_DIM_X = 0, unsigned TEAMS_PER_BLOCK_M = 0, unsigned TEAMS_PER_BLOCK_N = 0>
struct BsrMatrixSpMVTensorCoreFunctor {
  typedef nvcuda::wmma::accumulator accumulator;
  typedef nvcuda::wmma::row_major row_major;
  typedef nvcuda::wmma::col_major col_major;
  typedef nvcuda::wmma::matrix_a matrix_a;
  typedef nvcuda::wmma::matrix_b matrix_b;
  using FragA = nvcuda::wmma::fragment<matrix_a, FRAG_M, FRAG_N, FRAG_K, AFragScalar, row_major>;
  using FragX = nvcuda::wmma::fragment<matrix_b, FRAG_M, FRAG_N, FRAG_K, XFragScalar, row_major>;
  using FragY = nvcuda::wmma::fragment<accumulator, FRAG_M, FRAG_N, FRAG_K, YFragScalar>;

  typedef typename AMatrix::device_type Device;
  typedef Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef typename AMatrix::value_type AScalar;
  typedef typename YMatrix::value_type YScalar;
  typedef typename XMatrix::value_type XScalar;
  typedef typename AMatrix::non_const_ordinal_type AOrdinal;
  typedef typename AMatrix::non_const_size_type AOffset;

  // views of the shared memory used in the functor to cast types to the CUDA
  // wmma types A matrix is MxK X matrix is KxN Y matrix is MxN
  typedef typename Kokkos::View<AFragScalar *[FRAG_M][FRAG_K],  // one fragment per warp in the team (2D
                                                                // grid of warps in team)
                                Kokkos::LayoutRight, typename Device::execution_space::scratch_memory_space,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      AScratchView;
  typedef typename Kokkos::View<XFragScalar *[FRAG_K][FRAG_N],
                                typename Kokkos::LayoutRight,  // so that [FRAG_K][FRAG_N] part is
                                                               // contiguous in memory
                                typename Device::execution_space::scratch_memory_space,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      XScratchView;
  typedef typename Kokkos::View<YFragScalar **[FRAG_M][FRAG_N],
                                typename Kokkos::LayoutRight,  // so that [FRAG_M][FRAG_N] part is
                                                               // contiguous in memory
                                typename Device::execution_space::scratch_memory_space,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      YScratchView;

  YScalar alpha;
  AMatrix a;
  XMatrix x;
  YScalar beta;
  YMatrix y;

  BsrMatrixSpMVTensorCoreFunctorParams params;

  // a team is a 2D grid of warps
  static constexpr int WARPS_PER_TEAM_X = 2;
  static constexpr int WARPS_PER_TEAM_Y = 2;
  static constexpr int THREADS_PER_WARP = 32;

  BsrMatrixSpMVTensorCoreFunctor() = delete;  // need all runtime parameters

  // the launch parameters should be generated by a call to ::launch_parameters
  BsrMatrixSpMVTensorCoreFunctor(YScalar _alpha, AMatrix _a, XMatrix _x, YScalar _beta, YMatrix _y,
                                 const BsrMatrixSpMVTensorCoreFunctorParams &_params)
      : alpha(_alpha), a(_a), x(_x), beta(_beta), y(_y), params(_params) {}

  size_t league_size() const { return params.leagueDim_x * params.leagueDim_y; }

  size_t team_size() const { return THREADS_PER_WARP * WARPS_PER_TEAM_X * WARPS_PER_TEAM_Y; }

  // single column of fragments from A
  KOKKOS_INLINE_FUNCTION size_t a_scratch_size() const {
    return WARPS_PER_TEAM_Y * FRAG_M * FRAG_K * sizeof(AFragScalar);
  }
  // single row of fragments from X
  KOKKOS_INLINE_FUNCTION size_t x_scratch_size() const {
    return WARPS_PER_TEAM_X * FRAG_K * FRAG_N * sizeof(XFragScalar);
  }
  // one fragment per warp in the team
  KOKKOS_INLINE_FUNCTION size_t y_scratch_size() const {
    return WARPS_PER_TEAM_X * WARPS_PER_TEAM_Y * FRAG_M * FRAG_N * sizeof(YFragScalar);
  }

  size_t team_scratch_size() const { return a_scratch_size() + x_scratch_size() + y_scratch_size(); }

  /// \brief determine the mapping parameters for the 1D Kokkos::parallel_for
  /// space to the hierarchical 2D space of the functor kernel. This should be
  /// called to determine what arguments to pass to the constructor
  static BsrMatrixSpMVTensorCoreFunctorParams launch_parameters(const YScalar & /*alpha*/, const AMatrix &a,
                                                                const XMatrix & /*x*/, const YScalar & /*beta*/,
                                                                const YMatrix &y) {
    BsrMatrixSpMVTensorCoreFunctorParams params;

    // compute how many blocks there are in each dimension of the product MV
    int blocksPerYM = (y.extent(0) + a.blockDim() - 1) / a.blockDim();
    int blocksPerYN = (y.extent(1) + a.blockDim() - 1) / a.blockDim();

    // compute how many fragments are needed to cover each block
    int fragsPerBlockM = (a.blockDim() + FRAG_M - 1) / FRAG_M;
    int fragsPerBlockN = (a.blockDim() + FRAG_N - 1) / FRAG_N;

    // determine how many teams will need to cover each block (Y in M direction,
    // X in N direction)
    params.teamsPerBlockM = (fragsPerBlockM + WARPS_PER_TEAM_Y - 1) / WARPS_PER_TEAM_Y;
    params.teamsPerBlockN = (fragsPerBlockN + WARPS_PER_TEAM_X - 1) / WARPS_PER_TEAM_X;

    // determine how many teams will be needed co cover the product vector
    int yTeamsM = params.teamsPerBlockM * blocksPerYM;
    int yTeamsN = params.teamsPerBlockN * blocksPerYN;

    // Y dimension to M, X dimension to N
    params.leagueDim_x = yTeamsN;
    params.leagueDim_y = yTeamsM;

    return params;
  }

  // execute the functor with provided launch parameters
  void dispatch(const execution_space &exec) {
    typename BsrMatrixSpMVTensorCoreFunctor::team_policy policy(exec, league_size(), team_size());
    policy.set_scratch_size(0, Kokkos::PerTeam(team_scratch_size()));
    Kokkos::parallel_for("KokkosSparse::Experimental::BsrMatrixSpMVTensorCoreFunctor", policy, *this);
  }

  /*
     Consider the product vector as being made up of blocks that are the
     same size as the blocks in the input sparse matrix.
     teams are tiled across each block
     the size of each team is determined by the 2D grid of warps in the team,
     and the shape of each warp's fragment

     The number of warps per team is static:
     WARPS_PER_TEAM_X * WARPS_PER_TEAM_Y

     Based on its position in the product vector, each team steps through
     corresponding block-sized tiles of A and X. Those tiles are loaded into
     shared memory
     Warps in the team iterate over the shared-memory tile and perform the
     accumuation
     then the fragments write the results back to shared memory, and then
     global memory
  */

  KOKKOS_INLINE_FUNCTION void operator()(const team_member &mbr) const {
    using nvcuda::wmma::fill_fragment;
    using nvcuda::wmma::load_matrix_sync;
    using nvcuda::wmma::mma_sync;
    using nvcuda::wmma::store_matrix_sync;

    FragA fa;
    FragX fx;
    FragY fy;

    // override with template params if given
    const int ld_x = LEAGUE_DIM_X > 0 ? LEAGUE_DIM_X : params.leagueDim_x;
    const int tpbn = TEAMS_PER_BLOCK_N > 0 ? TEAMS_PER_BLOCK_N : params.teamsPerBlockN;
    const int tpbm = TEAMS_PER_BLOCK_M > 0 ? TEAMS_PER_BLOCK_M : params.teamsPerBlockM;

    // which team I am in the league
    const int teamIdx_x = mbr.league_rank() % ld_x;
    const int teamIdx_y = mbr.league_rank() / ld_x;

    // which block I contribute to in the product vector
    const int blockIdx_i = teamIdx_y / tpbm;
    const int blockIdx_j = teamIdx_x / tpbn;

    // which team am I in the block
    const int teamIdx_i = teamIdx_y % tpbm;
    const int teamIdx_j = teamIdx_x % tpbn;

    // which warp I am in the team
    const int warpIdx_x = (mbr.team_rank() / 32) % WARPS_PER_TEAM_X;
    const int warpIdx_y = (mbr.team_rank() / 32) / WARPS_PER_TEAM_X;

    // which lane I am in the warp
    const int lx = mbr.team_rank() % THREADS_PER_WARP;

    // which row of a/y the fragment this warp contributes to starts at
    const AOrdinal ay_i = blockIdx_i * a.blockDim()                // offset due to block
                          + teamIdx_i * WARPS_PER_TEAM_Y * FRAG_M  // offset of team within block
                          + warpIdx_y * FRAG_M;                    // offset of warp within team

    // which column of x/y the fragments warp will read from/contribute to
    // starts at
    const AOrdinal xy_j = blockIdx_j * a.blockDim() + teamIdx_j * WARPS_PER_TEAM_X * FRAG_N + warpIdx_x * FRAG_N;

    AFragScalar *_sa = (AFragScalar *)mbr.team_shmem().get_shmem(a_scratch_size());
    XFragScalar *_sx = (XFragScalar *)mbr.team_shmem().get_shmem(x_scratch_size());
    YFragScalar *_sy = (YFragScalar *)mbr.team_shmem().get_shmem(y_scratch_size());

    AScratchView sa(_sa, WARPS_PER_TEAM_Y);
    XScratchView sx(_sx, WARPS_PER_TEAM_X);
    YScratchView sy(_sy, WARPS_PER_TEAM_Y, WARPS_PER_TEAM_X);

    // team loads its fragments of Y that make up part or all of the block of Y
    // it's responsible for. each warp loads the part corresponding to its y
    // fragment stage through shared memory to convert to fragment type

    // no need for a team barrier because each warp uses an individual part of
    // shared memory
    for (unsigned i = lx; i < FRAG_M * FRAG_N; i += THREADS_PER_WARP) {
      const unsigned fi = i / FRAG_N;  // position in fragment of Y
      const unsigned fj = i % FRAG_N;
      const AOrdinal bi = teamIdx_i * WARPS_PER_TEAM_Y * FRAG_M + warpIdx_y * FRAG_M + fi;  // position in block of Y
      const AOrdinal bj = teamIdx_j * WARPS_PER_TEAM_X * FRAG_N + warpIdx_x * FRAG_N + fj;

      // load 0 outside of the block boundary and y vector boundary
      // load 0 outside of the vector boundary
      if (bi < a.blockDim() && bj < a.blockDim() && xy_j + fj < y.extent(1)) {
        sy(warpIdx_y, warpIdx_x, fi, fj) = YFragScalar(beta * y(ay_i + fi, xy_j + fj));
      } else {
        sy(warpIdx_y, warpIdx_x, fi, fj) = YFragScalar(0);
      }
    }
    // no barrier - each warp uses independent shared memory

    // load from the shared memory
    load_matrix_sync(fy, &sy(warpIdx_y, warpIdx_x, 0, 0), FRAG_N, nvcuda::wmma::mem_row_major);

    auto rowView = a.block_row_Const(blockIdx_i);

    // team loops through all blocks in the row
    for (AOffset ci = a.graph.row_map(blockIdx_i); ci < a.graph.row_map(blockIdx_i + 1); ++ci) {
      AOrdinal j = a.graph.entries(ci);

      // pointer to the beginning of the block
      const AScalar *ap = nullptr;
      {
        size_t off = ci - a.graph.row_map(blockIdx_i);    // which block in this row
        ap         = rowView.local_row_in_block(off, 0);  // offset of this block
      }

      // the block may be bigger than a single team,
      // each team is only one fragment long in the K direction
      // so will need to iterate fragments in the K direction across the block
      // the team will collaboratively load the fragments from A and X

      // and require multiple loads and accumulates
      // for mxn grid of fragments in the product vector, we need m rows of
      // fragments from A and n cols of fragments from X. only hold part of a
      // single column of fragments (from A) or part of a single row (from X) at
      // once
      for (AOrdinal bk = 0; bk < a.blockDim(); bk += FRAG_K /*M*/) {
        // team collaborative load of A
        // the footprint is one fragment wide in K direction
        mbr.team_barrier();
        for (unsigned i = mbr.team_rank(); i < WARPS_PER_TEAM_Y * FRAG_M * FRAG_K; i += mbr.team_size()) {
          const unsigned ti = i / FRAG_K;  // offset inside the fragments
          const unsigned tj = i % FRAG_K;
          // add in offset within block
          const AOrdinal bi = teamIdx_i * WARPS_PER_TEAM_Y * FRAG_M + ti;
          const AOrdinal bj = bk + tj;

          // fill shmem with 0 outside of the block boundary
          if (bi < a.blockDim() && bj < a.blockDim()) {
            sa(ti / FRAG_M, ti % FRAG_M, tj) = AFragScalar(alpha * ap[bi * a.blockDim() + bj]);
          } else {
            sa(ti / FRAG_M, ti % FRAG_M, tj) = AFragScalar(0);
          }
        }

        // collaborative load of X fragments into shared memory
        // entire team loads fragment footprint
        for (unsigned i = mbr.team_rank(); i < WARPS_PER_TEAM_X * FRAG_N * FRAG_K; i += mbr.team_size()) {
          const unsigned ti = i / (WARPS_PER_TEAM_X * FRAG_N);  // position in combined tiles
          const unsigned tj = i % (WARPS_PER_TEAM_X * FRAG_N);

          // add in offset within block
          const AOrdinal bi = bk + ti;
          const AOrdinal bj = teamIdx_j * WARPS_PER_TEAM_X * FRAG_N + tj;

          // load 0 outside of the block boundary
          // x is not necessarily a multiple of block size, so make sure access
          // is in bounds
          if (bi < a.blockDim() && bj < a.blockDim() && unsigned(blockIdx_j * a.blockDim() + bj) < x.extent(1)) {
            // tile is some fragments in the j/n direction that are frag_n wide
            sx(tj / FRAG_N, ti, tj % FRAG_N) = XFragScalar(x(j * a.blockDim() + bi, blockIdx_j * a.blockDim() + bj));
          } else {
            sx(tj / FRAG_N, ti, tj % FRAG_N) = XFragScalar(0);
          }
        }
        mbr.team_barrier();

        // load correct fragment from shared memory and accumulate
        // only need to do any math if our fragment will write a result back to
        // Y
        if (ay_i < static_cast<AOrdinal>(y.extent(0)) && xy_j < static_cast<AOrdinal>(y.extent(1))) {
          load_matrix_sync(fa, &sa(warpIdx_y, 0, 0), FRAG_K);
          load_matrix_sync(fx, &sx(warpIdx_x, 0, 0), FRAG_N);
          mma_sync(fy, fa, fx, fy);
        }
      }
    }  // loop through blocks in row of A

    // store Y fragments into shared memory
    store_matrix_sync(&sy(warpIdx_y, warpIdx_x, 0, 0), fy, FRAG_N, nvcuda::wmma::mem_row_major);
    // team loads its fragments of Y that make up part or all of the block of Y
    // it's responsible for. each warp loads the part corresponding to its y
    // fragment
    mbr.team_barrier();
    for (unsigned i = lx; i < FRAG_M * FRAG_N; i += THREADS_PER_WARP) {
      const unsigned fi = i / FRAG_N;  // position in fragment of Y
      const unsigned fj = i % FRAG_N;
      const AOrdinal bi = teamIdx_i * WARPS_PER_TEAM_Y * FRAG_M + warpIdx_y * FRAG_M + fi;  // position in block of Y
      const AOrdinal bj = teamIdx_j * WARPS_PER_TEAM_X * FRAG_N + warpIdx_x * FRAG_N + fj;

      // only store inside the block boundary
      // FIXME: what if Y is not wide enough? check y(_, j)
      if (bi < a.blockDim() && bj < a.blockDim() && xy_j + fj < y.extent(1)) {
        y(ay_i + fi, xy_j + fj) = sy(warpIdx_y, warpIdx_x, fi, fj);
      }
    }
    mbr.team_barrier();
  }
};

/// \brief Avoid instantiating tensor core functor for unsupported types
///
/// Instantiate some common template parameter values
/// for BsrMatrixSpMVTensorCoreFunctor.
/// This is a struct instead of a function for template...using shorthand
/// Discriminates between non-complex/on-GPU (supported) and otherwise
/// (unsupported) scalar types, and throws a runtime error for unsupported types
template <typename execution_space, typename AMatrix,
          typename AFragScalar,  // input matrix type and fragment scalar type
          typename XMatrix, typename XFragScalar, typename YMatrix, typename YFragScalar, unsigned FRAG_M,
          unsigned FRAG_N, unsigned FRAG_K>
struct BsrMatrixSpMVTensorCoreDispatcher {
  typedef typename AMatrix::value_type AScalar;
  typedef typename YMatrix::value_type YScalar;
  typedef typename XMatrix::value_type XScalar;

  template <unsigned X, unsigned Y, unsigned Z>
  using Dyn = BsrMatrixSpMVTensorCoreFunctor<execution_space, AMatrix, AFragScalar, XMatrix, XFragScalar, YMatrix,
                                             YFragScalar, FRAG_M, FRAG_N, FRAG_K, X, Y, Z>;

  // to be used when the various matrix types are supported
  static void tag_dispatch(std::true_type, const execution_space &exec, const YScalar alpha, AMatrix a, XMatrix x,
                           YScalar beta, YMatrix y) {
    BsrMatrixSpMVTensorCoreFunctorParams params = Dyn<0, 0, 0>::launch_parameters(alpha, a, x, beta, y);

    if (false) {  // consistency of formatting for next sections
    } else if (1 == params.leagueDim_x && 1 == params.teamsPerBlockM && 1 == params.teamsPerBlockN) {
      Dyn<1, 1, 1>(alpha, a, x, beta, y, params).dispatch(exec);
    } else if (1 == params.leagueDim_x && 2 == params.teamsPerBlockM && 2 == params.teamsPerBlockN) {
      Dyn<1, 2, 2>(alpha, a, x, beta, y, params).dispatch(exec);
    } else if (1 == params.leagueDim_x && 4 == params.teamsPerBlockM && 4 == params.teamsPerBlockN) {
      Dyn<1, 4, 4>(alpha, a, x, beta, y, params).dispatch(exec);
    } else if (1 == params.leagueDim_x && 8 == params.teamsPerBlockM && 8 == params.teamsPerBlockN) {
      Dyn<1, 8, 8>(alpha, a, x, beta, y, params).dispatch(exec);
    } else if (2 == params.leagueDim_x && 1 == params.teamsPerBlockM && 1 == params.teamsPerBlockN) {
      Dyn<2, 1, 1>(alpha, a, x, beta, y, params).dispatch(exec);
    } else if (2 == params.leagueDim_x && 2 == params.teamsPerBlockM && 2 == params.teamsPerBlockN) {
      Dyn<2, 2, 2>(alpha, a, x, beta, y, params).dispatch(exec);
    } else if (2 == params.leagueDim_x && 4 == params.teamsPerBlockM && 4 == params.teamsPerBlockN) {
      Dyn<2, 4, 4>(alpha, a, x, beta, y, params).dispatch(exec);
    } else if (2 == params.leagueDim_x && 8 == params.teamsPerBlockM && 8 == params.teamsPerBlockN) {
      Dyn<2, 8, 8>(alpha, a, x, beta, y, params).dispatch(exec);
    } else {
      Dyn<0, 0, 0>(alpha, a, x, beta, y, params).dispatch(exec);
    }
  }

  // to be used to avoid instantiating on unsupported types
  static void tag_dispatch(std::false_type, const execution_space &, YScalar, AMatrix, XMatrix, YScalar, YMatrix) {
    KokkosKernels::Impl::throw_runtime_exception(
        "Tensor core SpMV is only supported for non-complex types in GPU "
        "execution spaces");
  }
  static void dispatch(const execution_space &exec, YScalar alpha, AMatrix a, XMatrix x, YScalar beta, YMatrix y) {
    // tag will be false unless all conditions are met
    using tag = std::integral_constant<bool, TensorCoresAvailable<execution_space, AMatrix, XMatrix, YMatrix>::value>;
    tag_dispatch(tag{}, exec, alpha, a, x, beta, y);
  }
};

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // #if CUDA && (VOLTA || AMPERE)

//
//
//

#include "KokkosBlas.hpp"
#include "KokkosBlas2_serial_gemv_internal.hpp"
#include "KokkosBlas2_team_gemv_impl.hpp"
#include "KokkosBatched_Gemm_Serial_Internal.hpp"
#include "KokkosBatched_Gemm_TeamVector_Internal.hpp"
#include "KokkosBlas1_team_scal_impl.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"

namespace KokkosSparse {
namespace Impl {
namespace Bsr {

template <class AMatrix, class XVector, class YVector>
struct BSR_GEMV_Functor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;

  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename AMatrix::non_const_size_type size_type;

  const value_type alpha;
  const value_type beta;

  AMatrix m_A;
  XVector m_x;
  YVector m_y;

  const int block_dim;
  const bool conjugate;

  BSR_GEMV_Functor(const value_type alpha_, const AMatrix &m_A_, const XVector &m_x_, const value_type beta_,
                   YVector &m_y_, const int block_dim_, const bool conj_)
      : alpha(alpha_), beta(beta_), m_A(m_A_), m_x(m_x_), m_y(m_y_), block_dim(block_dim_), conjugate(conj_) {
    static_assert(static_cast<int>(XVector::rank) == 1, "XVector must be a rank 1 View.");
    static_assert(static_cast<int>(YVector::rank) == 1, "YVector must be a rank 1 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iBlock) const {
    const auto ystart        = iBlock * block_dim;
    const auto start         = m_A.graph.row_map(iBlock);
    const ordinal_type count = static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
    const auto row           = m_A.block_row_Const(iBlock);
    const auto beta1         = static_cast<value_type>(1);
    //
    if (conjugate) {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto xstart = row.block_colidx(ic) * block_dim;
        for (ordinal_type ii = 0; ii < block_dim; ++ii) {
          value_type t(0);
          for (ordinal_type jj = 0; jj < block_dim; ++jj) {
            const auto aval = Kokkos::ArithTraits<value_type>::conj(Aview(ii, jj));
            t += aval * m_x(xstart + jj);
          }
          m_y(ystart + ii) += alpha * t;
        }
      }
    } else {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto xstart = row.block_colidx(ic) * block_dim;
        KokkosBlas::Impl::SerialGemvInternal<KokkosBlas::Algo::Gemv::Blocked>::invoke<value_type, value_type>(
            block_dim, block_dim, alpha, Aview.data(), block_dim, 1, &m_x(xstart), static_cast<int>(m_x.stride_0()),
            beta1, &m_y(ystart), static_cast<int>(m_y.stride_0()));
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &dev) const {
    using y_value_type        = typename YVector::non_const_value_type;
    const ordinal_type iBlock = static_cast<ordinal_type>(dev.league_rank());

    const size_type Y_ptBeg = iBlock * block_dim;
    const size_type Y_ptEnd = Y_ptBeg + block_dim;
    auto Y_cur              = Kokkos::subview(m_y, ::Kokkos::make_pair(Y_ptBeg, Y_ptEnd));

    const y_value_type val_one = Kokkos::ArithTraits<y_value_type>::one();
    ;
    if (beta != val_one) {
      KokkosBlas::Impl::TeamVectorScaleInternal::invoke(dev, block_dim, beta, Y_cur.data(),
                                                        static_cast<int>(Y_cur.stride_0()));
    }

    dev.team_barrier();

    const auto myRow = m_A.block_row_Const(iBlock);
    const auto count = myRow.length;

    if (conjugate) {
      for (ordinal_type jBlock = 0; jBlock < count; ++jBlock) {
        const auto A_cur    = myRow.block(jBlock);
        const auto X_blkCol = myRow.block_colidx(jBlock);
        const auto X_ptBeg  = X_blkCol * block_dim;
        const auto X_cur    = Kokkos::subview(m_x, ::Kokkos::make_pair(X_ptBeg, X_ptBeg + block_dim));
        KokkosBlas::Impl::TeamVectorGemvInternal<KokkosBlas::Algo::Gemv::Unblocked>::invoke(
            dev, KokkosBlas::Impl::OpConj{}, A_cur.extent(0), A_cur.extent(1), alpha, A_cur.data(), A_cur.stride_0(),
            A_cur.stride_1(), X_cur.data(), X_cur.stride_0(), val_one, Y_cur.data(), Y_cur.stride_0());
      }
    } else {
      for (ordinal_type jBlock = 0; jBlock < count; ++jBlock) {
        const auto A_cur    = myRow.block(jBlock);
        const auto X_blkCol = myRow.block_colidx(jBlock);
        const auto X_ptBeg  = X_blkCol * block_dim;
        const auto X_cur    = Kokkos::subview(m_x, ::Kokkos::make_pair(X_ptBeg, X_ptBeg + block_dim));
        KokkosBlas::Impl::TeamVectorGemvInternal<KokkosBlas::Algo::Gemv::Unblocked>::invoke(
            dev, block_dim, block_dim, alpha, A_cur.data(), static_cast<int>(A_cur.stride_0()),
            static_cast<int>(A_cur.stride_1()), X_cur.data(), static_cast<int>(X_cur.stride_0()), val_one, Y_cur.data(),
            static_cast<int>(Y_cur.stride_0()));
      }
    }
  }
};

/* ******************* */

//
// spMatVec_no_transpose: version for CPU execution spaces
// (RangePolicy or trivial serial impl used)
//
template <class Handle, class AT, class AO, class AD, class AS, class AlphaType, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<typename YVector::execution_space>()>::type
              * = nullptr>
void spMatVec_no_transpose(
    const typename AD::execution_space &exec, Handle *handle, const AlphaType &alpha,
    const KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  // This is required to maintain semantics of KokkosKernels native SpMV:
  // if y contains NaN but beta = 0, the result y should be filled with 0.
  // For example, this is useful for passing in uninitialized y and beta=0.
  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(exec, y, Kokkos::ArithTraits<BetaType>::zero());
  else if (beta != Kokkos::ArithTraits<BetaType>::one())
    KokkosBlas::scal(exec, y, beta, y);

  //
  // Treat the case y <- alpha * A * x + beta * y
  //

  typedef KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;

  BSR_GEMV_Functor<AMatrix_Internal, XVector, YVector> func(alpha, A, x, beta, y, A.blockDim(), useConjugate);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::parallel_for(
        "KokkosSparse::bspmv<NoTranspose,Dynamic>",
        Kokkos::RangePolicy<typename AMatrix_Internal::device_type::execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(
            exec, 0, A.numRows()),
        func);
  } else {
    Kokkos::parallel_for(
        "KokkosSparse::bspmv<NoTranspose,Static>",
        Kokkos::RangePolicy<typename AMatrix_Internal::device_type::execution_space, Kokkos::Schedule<Kokkos::Static>>(
            exec, 0, A.numRows()),
        func);
  }
}

/* ******************* */

//
// spMatVec_no_transpose: version for GPU execution spaces (TeamPolicy used)
//
template <class Handle, class AT, class AO, class AD, class AS, class AlphaType, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<typename YVector::execution_space>()>::type
              * = nullptr>
void spMatVec_no_transpose(
    const typename AD::execution_space &exec, Handle *handle, const AlphaType &alpha,
    const KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  if (A.numRows() <= static_cast<AO>(0)) {
    return;
  }

  typedef KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;
  typedef typename AMatrix_Internal::execution_space execution_space;

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;

  int team_size        = -1;
  int vector_length    = -1;
  const auto block_dim = A.blockDim();

  team_size = 8;
  if (block_dim <= 4) {
    vector_length = 4;
    team_size     = 64;
  } else if (block_dim <= 8) {
    vector_length = 8;
    team_size     = 32;
  } else if (block_dim <= 16) {
    vector_length = 16;
    team_size     = 16;
  } else {
    vector_length = 32;
    team_size     = 8;
  }
  int64_t worksets = A.numRows();

  //
  // Use the handle to allow the user to pass in some tuning parameters.
  //
  if (handle->team_size != -1) team_size = handle->team_size;
  if (handle->vector_length != -1) vector_length = handle->vector_length;

  BSR_GEMV_Functor<AMatrix_Internal, XVector, YVector> func(alpha, A, x, beta, y, block_dim, useConjugate);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>> policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, Kokkos::AUTO,
                                                                                      vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, team_size,
                                                                                      vector_length);
    Kokkos::parallel_for("KokkosSparse::bspmv<NoTranspose,Dynamic>", policy, func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>> policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, Kokkos::AUTO,
                                                                                     vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, team_size,
                                                                                     vector_length);
    Kokkos::parallel_for("KokkosSparse::bspmv<NoTranspose, Static>", policy, func);
  }
}

/* ******************* */

template <class AMatrix, class XVector, class YVector>
struct BSR_GEMV_Transpose_Functor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;

  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename AMatrix::non_const_size_type size_type;

  const value_type alpha;

  AMatrix m_A;
  XVector m_x;
  YVector m_y;

  const int block_dim;
  const bool conjugate;

  BSR_GEMV_Transpose_Functor(const value_type alpha_, const AMatrix &m_A_, const XVector &m_x_, const YVector &m_y_,
                             const bool conj_)
      : alpha(alpha_), m_A(m_A_), m_x(m_x_), m_y(m_y_), block_dim(static_cast<int>(m_A_.blockDim())), conjugate(conj_) {
    static_assert(static_cast<int>(XVector::rank) == 1, "XVector must be a rank 1 View.");
    static_assert(static_cast<int>(YVector::rank) == 1, "YVector must be a rank 1 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iBlock) const {
    //
    // Assume that alpha is not zero
    //
    const auto xstart        = iBlock * block_dim;
    const auto xview         = Kokkos::subview(m_x, Kokkos::make_pair(xstart, xstart + block_dim));
    const auto start         = m_A.graph.row_map(iBlock);
    const ordinal_type count = static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
    const auto row           = m_A.block_row_Const(iBlock);
    const auto beta1         = static_cast<value_type>(1);
    const auto alpha1        = beta1;
    if (conjugate) {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto ystart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jj = 0; jj < block_dim; ++jj) {
          value_type t(0);
          for (ordinal_type ii = 0; ii < block_dim; ++ii) {
            const auto aval = Kokkos::ArithTraits<value_type>::conj(Aview(ii, jj));
            t += aval * xview(ii);
          }
          t *= alpha;
          Kokkos::atomic_add(&m_y(ystart + jj), t);
        }
      }
    } else {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto ystart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jj = 0; jj < block_dim; ++jj) {
          value_type t(0);
          KokkosBlas::Impl::SerialGemvInternal<KokkosBlas::Algo::Gemv::Blocked>::invoke<value_type, value_type>(
              1, block_dim, alpha1, Aview.data() + jj, Aview.stride_1(), Aview.stride_0(), xview.data(),
              xview.stride_0(), beta1, &t, 1);
          t *= alpha;
          Kokkos::atomic_add(&m_y(ystart + jj), t);
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &dev) const {
    using y_value_type        = typename YVector::non_const_value_type;
    const ordinal_type iBlock = static_cast<ordinal_type>(dev.league_rank());

    const size_type X_ptBeg = iBlock * block_dim;
    const size_type X_ptEnd = X_ptBeg + block_dim;
    const auto X_cur        = Kokkos::subview(m_x, ::Kokkos::make_pair(X_ptBeg, X_ptEnd));

    const auto myRow = m_A.block_row_Const(iBlock);
    const auto count = myRow.length;

    const y_value_type val_zero = Kokkos::ArithTraits<y_value_type>::zero();
    y_value_type *shared_y      = (y_value_type *)dev.team_shmem().get_shmem(block_dim * sizeof(y_value_type));

    if (conjugate) {
      Kokkos::View<y_value_type *, typename AMatrix::device_type, Kokkos::MemoryUnmanaged> shared_view(shared_y,
                                                                                                       block_dim);
      for (ordinal_type jBlock = 0; jBlock < count; ++jBlock) {
        const auto A_cur = myRow.block(jBlock);
        //
        KokkosBlas::TeamVectorGemv<team_member, KokkosBlas::Trans::ConjTranspose,
                                   KokkosBlas::Algo::Gemv::Default>::invoke(dev, alpha, A_cur, X_cur, val_zero,
                                                                            shared_view);
        //
        dev.team_barrier();
        //
        const auto Y_blkCol = myRow.block_colidx(jBlock);
        const auto Y_ptBeg  = Y_blkCol * block_dim;
        auto Y_cur          = Kokkos::subview(m_y, ::Kokkos::make_pair(Y_ptBeg, Y_ptBeg + block_dim));
        Kokkos::parallel_for(Kokkos::TeamVectorRange(dev, 0, block_dim),
                             [&](const ordinal_type &ijk) { Kokkos::atomic_add(&Y_cur(ijk), shared_view(ijk)); });
        //
        dev.team_barrier();
      }
    } else {
      for (ordinal_type jBlock = 0; jBlock < count; ++jBlock) {
        const auto A_cur = myRow.block(jBlock);
        //
        KokkosBlas::Impl::TeamVectorGemvInternal<KokkosBlas::Algo::Gemv::Unblocked>::invoke(
            dev, block_dim, block_dim, alpha, A_cur.data(), static_cast<int>(A_cur.stride_1()),
            static_cast<int>(A_cur.stride_0()), X_cur.data(), static_cast<int>(X_cur.stride_0()), val_zero, shared_y,
            1);
        //
        dev.team_barrier();
        //
        const auto Y_blkCol = myRow.block_colidx(jBlock);
        const auto Y_ptBeg  = Y_blkCol * block_dim;
        auto Y_cur          = Kokkos::subview(m_y, ::Kokkos::make_pair(Y_ptBeg, Y_ptBeg + block_dim));
        Kokkos::parallel_for(Kokkos::TeamVectorRange(dev, 0, block_dim),
                             [&](const ordinal_type &ijk) { Kokkos::atomic_add(&Y_cur(ijk), shared_y[ijk]); });
        //
        dev.team_barrier();
      }
    }
  }
};

/* ******************* */

/// \brief  spMatVec_transpose: version for CPU execution spaces (RangePolicy or
/// trivial serial impl used)
template <class Handle, class AT, class AO, class AD, class AS, class AlphaType, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<typename YVector::execution_space>()>::type
              * = nullptr>
void spMatVec_transpose(
    const typename AD::execution_space &exec, Handle *handle, const AlphaType &alpha,
    const KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  // This is required to maintain semantics of KokkosKernels native SpMV:
  // if y contains NaN but beta = 0, the result y should be filled with 0.
  // For example, this is useful for passing in uninitialized y and beta=0.
  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(exec, y, Kokkos::ArithTraits<BetaType>::zero());
  else if (beta != Kokkos::ArithTraits<BetaType>::one())
    KokkosBlas::scal(exec, y, beta, y);

  if (alpha == Kokkos::ArithTraits<AlphaType>::zero()) return;

  //
  // Treat the case y <- alpha * A^T * x + beta * y
  //

  typedef KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;

  BSR_GEMV_Transpose_Functor<AMatrix_Internal, XVector, YVector> func(alpha, A, x, y, useConjugate);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::parallel_for(
        "KokkosSparse::bspmv<Transpose,Dynamic>",
        Kokkos::RangePolicy<typename AMatrix_Internal::device_type::execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(
            0, A.numRows()),
        func);
  } else {
    Kokkos::parallel_for(
        "KokkosSparse::bspmv<Transpose,Static>",
        Kokkos::RangePolicy<typename AMatrix_Internal::device_type::execution_space, Kokkos::Schedule<Kokkos::Static>>(
            0, A.numRows()),
        func);
  }
}

//
// spMatVec_transpose: version for GPU execution spaces (TeamPolicy used)
//
template <class Handle, class AMatrix, class AlphaType, class XVector, class BetaType, class YVector,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<typename YVector::execution_space>()>::type
              * = nullptr>
void spMatVec_transpose(const typename AMatrix::execution_space &exec, Handle *handle, const AlphaType &alpha,
                        const AMatrix &A, const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  if (A.numRows() <= 0) {
    return;
  }

  typedef typename AMatrix::execution_space execution_space;

  const auto block_dim = A.blockDim();

  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(exec, y, Kokkos::ArithTraits<BetaType>::zero());
  else if (beta != Kokkos::ArithTraits<BetaType>::one())
    KokkosBlas::scal(exec, y, beta, y);

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;
  int team_size             = -1;
  int vector_length         = -1;

  int64_t worksets = A.numRows();

  if (block_dim <= 4) {
    vector_length = 4;
    team_size     = 64;
  }
  if (block_dim <= 8) {
    vector_length = 8;
    team_size     = 32;
  }
  if (block_dim <= 16) {
    vector_length = 16;
    team_size     = 16;
  } else {
    vector_length = 32;
    team_size     = 8;
  }

  //
  // Use the handle to allow the user to pass in some tuning parameters.
  //
  if (handle->team_size != -1) team_size = handle->team_size;
  if (handle->vector_length != -1) vector_length = handle->vector_length;

  BSR_GEMV_Transpose_Functor<AMatrix, XVector, YVector> func(alpha, A, x, y, useConjugate);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>> policy(exec, 1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, Kokkos::AUTO,
                                                                                      vector_length)
                   .set_scratch_size(0, Kokkos::PerTeam(block_dim * sizeof(typename YVector::non_const_value_type)));
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, team_size,
                                                                                      vector_length)
                   .set_scratch_size(0, Kokkos::PerTeam(block_dim * sizeof(typename YVector::non_const_value_type)));
    Kokkos::parallel_for("KokkosSparse::bspmv<Transpose,Dynamic>", policy, func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>> policy(exec, 1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, Kokkos::AUTO,
                                                                                     vector_length)
                   .set_scratch_size(0, Kokkos::PerTeam(block_dim * sizeof(typename YVector::non_const_value_type)));
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, team_size,
                                                                                     vector_length)
                   .set_scratch_size(0, Kokkos::PerTeam(block_dim * sizeof(typename YVector::non_const_value_type)));
    Kokkos::parallel_for("KokkosSparse::bspmv<Transpose, Static>", policy, func);
  }
}

/* ******************* */

template <class AMatrix, class XVector, class YVector>
struct BSR_GEMM_Functor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;

  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename AMatrix::non_const_size_type size_type;

  const value_type alpha;
  const value_type beta;

  AMatrix m_A;
  XVector m_x;
  YVector m_y;

  const int block_dim;
  const int num_rhs;
  const bool conjugate;

  BSR_GEMM_Functor(const value_type alpha_, const AMatrix m_A_, const XVector m_x_, const value_type beta_,
                   const YVector m_y_, const bool conj_)
      : alpha(alpha_),
        beta(beta_),
        m_A(m_A_),
        m_x(m_x_),
        m_y(m_y_),
        block_dim(static_cast<int>(m_A_.blockDim())),
        num_rhs(static_cast<int>(m_x_.extent(1))),
        conjugate(conj_) {
    static_assert(static_cast<int>(XVector::rank) == 2, "XVector must be a rank 2 View.");
    static_assert(static_cast<int>(YVector::rank) == 2, "YVector must be a rank 2 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iBlock) const {
    //
    const auto ystart        = iBlock * block_dim;
    const auto start         = m_A.graph.row_map(iBlock);
    const ordinal_type count = static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
    const auto row           = m_A.block_row_Const(iBlock);
    const auto beta1         = static_cast<value_type>(1);
    const auto ldx           = m_x.stride_1();
    const auto ldy           = m_y.stride_1();
    //
    if (conjugate) {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto xstart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jr = 0; jr < num_rhs; ++jr) {
          for (ordinal_type ii = 0; ii < block_dim; ++ii) {
            value_type t(0);
            for (ordinal_type jj = 0; jj < block_dim; ++jj) {
              const auto aval = Kokkos::ArithTraits<value_type>::conj(Aview(ii, jj));
              t += aval * m_x(xstart + jj, jr);
            }
            m_y(ystart + ii, jr) += alpha * t;
          }
        }
      }
    } else {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto xstart = row.block_colidx(ic) * block_dim;
        KokkosBatched::SerialGemmInternal<KokkosBatched::Algo::Gemm::Blocked>::invoke<value_type, value_type>(
            static_cast<ordinal_type>(block_dim), static_cast<ordinal_type>(num_rhs),
            static_cast<ordinal_type>(block_dim), alpha, Aview.data(), Aview.stride_0(), Aview.stride_1(),
            &m_x(xstart, 0), m_x.stride_0(), ldx, beta1, &m_y(ystart, 0), m_y.stride_0(), ldy);
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &dev) const {
    using y_value_type        = typename YVector::non_const_value_type;
    const ordinal_type iBlock = static_cast<ordinal_type>(dev.league_rank());

    const size_type Y_ptBeg = iBlock * block_dim;
    const size_type Y_ptEnd = Y_ptBeg + block_dim;
    auto Y_cur              = Kokkos::subview(m_y, ::Kokkos::make_pair(Y_ptBeg, Y_ptEnd), Kokkos::ALL());

    const y_value_type val_one = Kokkos::ArithTraits<y_value_type>::one();
    if (beta != val_one) {
      KokkosBlas::Impl::TeamVectorScaleInternal::invoke(dev, block_dim, num_rhs, beta, Y_cur.data(),
                                                        static_cast<int>(Y_cur.stride_0()),
                                                        static_cast<int>(Y_cur.stride_1()));
    }

    dev.team_barrier();

    const auto myRow = m_A.block_row_Const(iBlock);
    const auto count = myRow.length;

    if (conjugate) {
      for (ordinal_type jBlock = 0; jBlock < count; ++jBlock) {
        const auto A_cur    = myRow.block(jBlock);
        const auto X_blkCol = myRow.block_colidx(jBlock);
        const auto X_ptBeg  = X_blkCol * block_dim;
        const auto X_cur    = Kokkos::subview(m_x, ::Kokkos::make_pair(X_ptBeg, X_ptBeg + block_dim), Kokkos::ALL());
        KokkosBatched::TeamVectorGemmInternal<KokkosBatched::Algo::Gemm::Unblocked, true>::invoke(
            dev, static_cast<int>(block_dim), static_cast<int>(num_rhs), static_cast<int>(block_dim), alpha,
            A_cur.data(), static_cast<int>(A_cur.stride_0()), static_cast<int>(A_cur.stride_1()), X_cur.data(),
            static_cast<int>(X_cur.stride_0()), static_cast<int>(X_cur.stride_1()), val_one, Y_cur.data(),
            static_cast<int>(Y_cur.stride_0()), static_cast<int>(Y_cur.stride_1()));
      }
    } else {
      for (ordinal_type jBlock = 0; jBlock < count; ++jBlock) {
        const auto A_cur    = myRow.block(jBlock);
        const auto X_blkCol = myRow.block_colidx(jBlock);
        const auto X_ptBeg  = X_blkCol * block_dim;
        const auto X_cur    = Kokkos::subview(m_x, ::Kokkos::make_pair(X_ptBeg, X_ptBeg + block_dim), Kokkos::ALL());
        KokkosBatched::TeamVectorGemmInternal<KokkosBatched::Algo::Gemm::Unblocked, false>::invoke(
            dev, block_dim, num_rhs, block_dim, alpha, A_cur.data(), static_cast<int>(A_cur.stride_0()),
            static_cast<int>(A_cur.stride_1()), X_cur.data(), static_cast<int>(X_cur.stride_0()),
            static_cast<int>(X_cur.stride_1()), val_one, Y_cur.data(), static_cast<int>(Y_cur.stride_0()),
            static_cast<int>(Y_cur.stride_1()));
      }
    }
  }
};

/* ******************* */

//
// spMatMultiVec_no_transpose: version for CPU execution spaces
// (RangePolicy or trivial serial impl used)
//
template <class Handle, class AT, class AO, class AD, class AS, class AlphaType, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<typename YVector::execution_space>()>::type
              * = nullptr>
void spMatMultiVec_no_transpose(
    const typename AD::execution_space &exec, Handle *handle, const AlphaType &alpha,
    const KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  // This is required to maintain semantics of KokkosKernels native SpMV:
  // if y contains NaN but beta = 0, the result y should be filled with 0.
  // For example, this is useful for passing in uninitialized y and beta=0.
  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(exec, y, Kokkos::ArithTraits<BetaType>::zero());
  else if (beta != Kokkos::ArithTraits<BetaType>::one())
    KokkosBlas::scal(exec, y, beta, y);
  //
  // Treat the case y <- alpha * A * x + beta * y
  //
  typedef KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;

  BSR_GEMM_Functor<AMatrix_Internal, XVector, YVector> func(alpha, A, x, beta, y, useConjugate);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::parallel_for(
        "KokkosSparse::bsr_spm_mv<NoTranspose,Dynamic>",
        Kokkos::RangePolicy<typename AMatrix_Internal::device_type::execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(
            0, A.numRows()),
        func);
  } else {
    Kokkos::parallel_for(
        "KokkosSparse::bsr_spm_mv<NoTranspose,Static>",
        Kokkos::RangePolicy<typename AMatrix_Internal::device_type::execution_space, Kokkos::Schedule<Kokkos::Static>>(
            0, A.numRows()),
        func);
  }
}

/* ******************* */

//
// spMatMultiVec_no_transpose: version for GPU execution spaces (TeamPolicy
// used)
//
template <class Handle, class AT, class AO, class AD, class AS, class AlphaType, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<typename YVector::execution_space>()>::type
              * = nullptr>
void spMatMultiVec_no_transpose(
    const typename AD::execution_space &exec, Handle *handle, const AlphaType &alpha,
    const KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  if (A.numRows() <= static_cast<AO>(0)) {
    return;
  }

  typedef KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;
  typedef typename AMatrix_Internal::execution_space execution_space;

  bool use_dynamic_schedule = handle->force_dynamic_schedule;  // Forces the use of a dynamic schedule
  bool use_static_schedule  = handle->force_static_schedule;   // Forces the use of a static schedule

  int team_size     = -1;
  int vector_length = -1;
  int64_t worksets  = A.numRows();

  const auto block_dim = A.blockDim();
  if (block_dim <= 4) {
    vector_length = 4;
    team_size     = 64;
  } else if (block_dim <= 8) {
    vector_length = 8;
    team_size     = 32;
  } else if (block_dim <= 16) {
    vector_length = 16;
    team_size     = 16;
  } else {
    vector_length = 32;
    team_size     = 8;
  }

  //
  // Use the handle to allow the user to pass in some tuning parameters.
  //
  if (handle->team_size != -1) team_size = handle->team_size;
  if (handle->vector_length != -1) vector_length = handle->vector_length;

  BSR_GEMM_Functor<AMatrix_Internal, XVector, YVector> func(alpha, A, x, beta, y, useConjugate);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>> policy(exec, 1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, Kokkos::AUTO,
                                                                                      vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, team_size,
                                                                                      vector_length);
    Kokkos::parallel_for("KokkosSparse::bsr_spm_mv<NoTranspose,Dynamic>", policy, func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>> policy(1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, Kokkos::AUTO,
                                                                                     vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, team_size,
                                                                                     vector_length);
    Kokkos::parallel_for("KokkosSparse::bsr_spm_mv<NoTranspose, Static>", policy, func);
  }
}

/* ******************* */
template <class execution_space, class AMatrix, class XVector, class YVector>
struct BSR_GEMM_Transpose_Functor {
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;

  //! Nonconst version of the type of column indices in the sparse matrix.
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  //! Nonconst version of the type of row offsets in the sparse matrix.
  typedef typename AMatrix::non_const_size_type size_type;

  const value_type alpha;
  AMatrix m_A;
  XVector m_x;
  YVector m_y;

  const int block_dim;
  const int num_rhs;
  const bool conjugate;

  BSR_GEMM_Transpose_Functor(const value_type alpha_, const AMatrix &m_A_, const XVector &m_x_, YVector &m_y_,
                             const bool conj_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        m_y(m_y_),
        block_dim(static_cast<int>(m_A_.blockDim())),
        num_rhs(static_cast<int>(m_x_.extent(1))),
        conjugate(conj_) {
    static_assert(static_cast<int>(XVector::rank) == 2, "XVector must be a rank 2 View.");
    static_assert(static_cast<int>(YVector::rank) == 2, "YVector must be a rank 2 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const ordinal_type iBlock) const {
    //
    const auto xstart        = iBlock * block_dim;
    const auto xview         = Kokkos::subview(m_x, Kokkos::make_pair(xstart, xstart + block_dim), Kokkos::ALL());
    const auto start         = m_A.graph.row_map(iBlock);
    const ordinal_type count = static_cast<ordinal_type>(m_A.graph.row_map(iBlock + 1) - start);
    const auto row           = m_A.block_row_Const(iBlock);
    const auto beta1         = static_cast<value_type>(1);
    const auto alpha1        = beta1;
    const auto ldx           = m_x.stride_1();
    //
    if (conjugate) {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto ystart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jr = 0; jr < num_rhs; ++jr) {
          for (ordinal_type jj = 0; jj < block_dim; ++jj) {
            value_type t(0);
            for (ordinal_type ii = 0; ii < block_dim; ++ii) {
              const auto aval = Kokkos::ArithTraits<value_type>::conj(Aview(ii, jj));
              t += aval * xview(ii, jr);
            }
            t *= alpha;
            Kokkos::atomic_add(&m_y(ystart + jj, jr), t);
          }
        }
      }
    } else {
      for (ordinal_type ic = 0; ic < count; ++ic) {
        const auto Aview  = row.block(ic);
        const auto ystart = row.block_colidx(ic) * block_dim;
        for (ordinal_type jr = 0; jr < num_rhs; ++jr) {
          for (ordinal_type jj = 0; jj < block_dim; ++jj) {
            value_type t(0);
            KokkosBlas::Impl::SerialGemvInternal<KokkosBlas::Algo::Gemv::Blocked>::invoke<value_type, value_type>(
                1, block_dim, alpha1, Aview.data() + jj, Aview.stride_1(), Aview.stride_0(), xview.data() + jr * ldx,
                xview.stride_0(), beta1, &t, 1);
            t *= alpha;
            Kokkos::atomic_add(&m_y(ystart + jj, jr), t);
          }
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const team_member &dev) const {
    using y_value_type        = typename YVector::non_const_value_type;
    const ordinal_type iBlock = static_cast<ordinal_type>(dev.league_rank());

    const size_type X_ptBeg = iBlock * block_dim;
    const size_type X_ptEnd = X_ptBeg + block_dim;
    const auto X_cur        = Kokkos::subview(m_x, ::Kokkos::make_pair(X_ptBeg, X_ptEnd), Kokkos::ALL());

    const auto myRow = m_A.block_row_Const(iBlock);
    const auto count = myRow.length;

    const y_value_type val_zero = Kokkos::ArithTraits<y_value_type>::zero();
    y_value_type *shared_y = (y_value_type *)dev.team_shmem().get_shmem(block_dim * num_rhs * sizeof(y_value_type));

    if (conjugate) {
      for (ordinal_type jBlock = 0; jBlock < count; ++jBlock) {
        const auto A_cur = myRow.block(jBlock);
        //
        KokkosBatched::TeamVectorGemmInternal<KokkosBatched::Algo::Gemm::Unblocked, true>::invoke(
            dev, block_dim, num_rhs, block_dim, alpha, A_cur.data(), static_cast<int>(A_cur.stride_1()),
            static_cast<int>(A_cur.stride_0()), X_cur.data(), static_cast<int>(X_cur.stride_0()),
            static_cast<int>(X_cur.stride_1()), val_zero, shared_y, 1, block_dim);
        //
        dev.team_barrier();
        //
        const auto Y_blkCol = myRow.block_colidx(jBlock);
        const auto Y_ptBeg  = Y_blkCol * block_dim;
        auto Y_cur          = Kokkos::subview(m_y, ::Kokkos::make_pair(Y_ptBeg, Y_ptBeg + block_dim), Kokkos::ALL());
        Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, num_rhs), [&](const ordinal_type &kc) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, 0, block_dim), [&](const ordinal_type &kr) {
            Kokkos::atomic_add(&Y_cur(kr, kc), shared_y[kr + kc * block_dim]);
          });
        });
        dev.team_barrier();
      }
    } else {
      for (ordinal_type jBlock = 0; jBlock < count; ++jBlock) {
        const auto A_cur = myRow.block(jBlock);
        //
        KokkosBatched::TeamVectorGemmInternal<KokkosBatched::Algo::Gemm::Unblocked, false>::invoke(
            dev, block_dim, num_rhs, block_dim, alpha, A_cur.data(), static_cast<int>(A_cur.stride_1()),
            static_cast<int>(A_cur.stride_0()), X_cur.data(), static_cast<int>(X_cur.stride_0()),
            static_cast<int>(X_cur.stride_1()), val_zero, shared_y, 1, block_dim);
        //
        dev.team_barrier();
        //
        const auto Y_blkCol = myRow.block_colidx(jBlock);
        const auto Y_ptBeg  = Y_blkCol * block_dim;
        auto Y_cur          = Kokkos::subview(m_y, ::Kokkos::make_pair(Y_ptBeg, Y_ptBeg + block_dim), Kokkos::ALL());
        Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, num_rhs), [&](const ordinal_type &kc) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, 0, block_dim), [&](const ordinal_type &kr) {
            Kokkos::atomic_add(&Y_cur(kr, kc), shared_y[kr + kc * block_dim]);
          });
        });
        dev.team_barrier();
      }
    }
  }
};

/* ******************* */

/// \brief  spMatMultiVec_transpose: version for CPU execution spaces
/// (RangePolicy or trivial serial impl used)
template <class execution_space, class Handle, class AT, class AO, class AD, class AS, class AlphaType, class XVector,
          class BetaType, class YVector,
          typename std::enable_if<!KokkosKernels::Impl::kk_is_gpu_exec_space<typename YVector::execution_space>()>::type
              * = nullptr>
void spMatMultiVec_transpose(
    const execution_space &exec, Handle *handle, const AlphaType &alpha,
    const KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS> &A,
    const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  // This is required to maintain semantics of KokkosKernels native SpMV:
  // if y contains NaN but beta = 0, the result y should be filled with 0.
  // For example, this is useful for passing in uninitialized y and beta=0.
  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(exec, y, Kokkos::ArithTraits<BetaType>::zero());
  else if (beta != Kokkos::ArithTraits<BetaType>::one())
    KokkosBlas::scal(exec, y, beta, y);
  //
  // Treat the case y <- alpha * A^T * x + beta * y
  //
  typedef KokkosSparse::Experimental::BsrMatrix<AT, AO, AD, Kokkos::MemoryTraits<Kokkos::Unmanaged>, AS>
      AMatrix_Internal;

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;

  BSR_GEMM_Transpose_Functor<execution_space, AMatrix_Internal, XVector, YVector> func(alpha, A, x, y, useConjugate);
  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::parallel_for("KokkosSparse::bsr_spm_mv<Transpose,Dynamic>",
                         Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, 0, A.numRows()),
                         func);
  } else {
    Kokkos::parallel_for("KokkosSparse::bsr_spm_mv<Transpose,Static>",
                         Kokkos::RangePolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, 0, A.numRows()),
                         func);
  }
}

//
// spMatMultiVec_transpose: version for GPU execution spaces (TeamPolicy used)
//
template <class execution_space, class Handle, class AMatrix, class AlphaType, class XVector, class BetaType,
          class YVector,
          typename std::enable_if<KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()>::type * = nullptr>
void spMatMultiVec_transpose(const execution_space &exec, Handle *handle, const AlphaType &alpha, const AMatrix &A,
                             const XVector &x, const BetaType &beta, YVector &y, bool useConjugate) {
  if (A.numRows() <= 0) {
    return;
  }

  if (beta == Kokkos::ArithTraits<BetaType>::zero())
    Kokkos::deep_copy(exec, y, Kokkos::ArithTraits<BetaType>::zero());
  else if (beta != Kokkos::ArithTraits<BetaType>::one())
    KokkosBlas::scal(exec, y, beta, y);

  bool use_dynamic_schedule = handle->force_dynamic_schedule;
  bool use_static_schedule  = handle->force_static_schedule;
  int team_size             = -1;
  int vector_length         = -1;
  int64_t worksets          = A.numRows();

  const auto block_dim = A.blockDim();
  if (block_dim <= 4) {
    vector_length = 4;
    team_size     = 64;
  } else if (block_dim <= 8) {
    vector_length = 8;
    team_size     = 32;
  } else if (block_dim <= 8) {
    vector_length = 16;
    team_size     = 16;
  } else {
    vector_length = 32;
    team_size     = 8;
  }

  //
  // Use the handle to allow the user to pass in some tuning parameters.
  //
  if (handle->team_size != -1) team_size = handle->team_size;
  if (handle->vector_length != -1) vector_length = handle->vector_length;

  BSR_GEMM_Transpose_Functor<execution_space, AMatrix, XVector, YVector> func(alpha, A, x, y, useConjugate);

  if (((A.nnz() > 10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>> policy(exec, 1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, Kokkos::AUTO,
                                                                                      vector_length)
                   .set_scratch_size(
                       0, Kokkos::PerTeam(block_dim * x.extent(1) * sizeof(typename YVector::non_const_value_type)));
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic>>(exec, worksets, team_size,
                                                                                      vector_length)
                   .set_scratch_size(
                       0, Kokkos::PerTeam(block_dim * x.extent(1) * sizeof(typename YVector::non_const_value_type)));
    Kokkos::parallel_for("KokkosSparse::bsr_spm_mv<Transpose,Dynamic>", policy, func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>> policy(exec, 1, 1);
    if (team_size < 0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, Kokkos::AUTO,
                                                                                     vector_length)
                   .set_scratch_size(
                       0, Kokkos::PerTeam(block_dim * x.extent(1) * sizeof(typename YVector::non_const_value_type)));
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static>>(exec, worksets, team_size,
                                                                                     vector_length)
                   .set_scratch_size(
                       0, Kokkos::PerTeam(block_dim * x.extent(1) * sizeof(typename YVector::non_const_value_type)));
    Kokkos::parallel_for("KokkosSparse::bsr_spm_mv<Transpose, Static>", policy, func);
  }
}

/* ******************* */

}  // namespace Bsr
}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_IMPL_SPMV_BSRMATRIX_IMPL_HPP_
