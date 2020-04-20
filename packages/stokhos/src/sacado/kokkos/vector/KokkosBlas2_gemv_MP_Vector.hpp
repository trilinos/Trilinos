// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef KOKKOSBLAS2_GEMV_MP_VECTOR_HPP
#define KOKKOSBLAS2_GEMV_MP_VECTOR_HPP

#include <type_traits>
#include "Sacado_ConfigDefs.h"

#include "Stokhos_ViewStorage.hpp"
#include "Sacado_MP_Vector.hpp"
#include "Kokkos_View_MP_Vector.hpp"
#include "Kokkos_ArithTraits_MP_Vector.hpp"
#include "KokkosBlas.hpp"

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Kokkos_Core.hpp"

#include "Stokhos_config.h"

#define Sacado_MP_Vector_GEMV_Tile_Size(size) (STOKHOS_GEMV_CACHE_SIZE / size)

// Functor for the update
template <class AViewType,
          class XViewType,
          class YViewType,
          class IndexType = typename AViewType::size_type>
struct updateF
{
    using AlphaCoeffType = typename AViewType::non_const_value_type;
    using BetaCoeffType = typename YViewType::non_const_value_type;

    updateF(const AlphaCoeffType &alpha,
            const AViewType &A,
            const XViewType &x,
            const BetaCoeffType &beta,
            const YViewType &y,
            const IndexType m_c) : alpha_(alpha), A_(A), x_(x), beta_(beta), y_(y), m_c_(m_c)
    {
    }

public:
    // i_tile is the current tile of the input matrix A
    KOKKOS_INLINE_FUNCTION void
    operator()(const IndexType &i_tile) const
    {
        const IndexType m = y_.extent(0);
        const IndexType n = x_.extent(0);

        IndexType i_min = m_c_ * i_tile;
        bool last_tile = (i_min + m_c_ >= m);
        IndexType i_max = (last_tile) ? m : (i_min + m_c_);

#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        if (beta_ == BetaCoeffType(0.))
            for (IndexType i = i_min; i < i_max; ++i)
                y_(i) = beta_;
        else
            for (IndexType i = i_min; i < i_max; ++i)
                y_(i) *= beta_;

        for (IndexType j = 0; j < n; ++j)
        {
            AlphaCoeffType alphab = alpha_ * x_(j);

            for (IndexType i = i_min; i < i_max; ++i)
                y_(i) += alphab * A_(i, j);
        }
    }

private:
    AlphaCoeffType alpha_;
    typename AViewType::const_type A_;
    typename XViewType::const_type x_;
    BetaCoeffType beta_;
    YViewType y_;
    const IndexType m_c_;
};

// Functor for the inner products
template <class AViewType,
          class XViewType,
          class YViewType,
          class IndexType = typename AViewType::size_type>
struct innerF
{
    using execution_space = typename AViewType::execution_space;
    using policy_type = Kokkos::TeamPolicy<execution_space>;
    using member_type = typename policy_type::member_type;

    using AlphaCoeffType = typename AViewType::non_const_value_type;
    using Scalar = AlphaCoeffType;

    innerF(const AlphaCoeffType &alpha,
           const AViewType &A,
           const XViewType &x,
           const YViewType &y,
           const IndexType n_c) : alpha_(alpha), A_(A), x_(x), y_(y), n_c_(n_c)
    {
    }

public:
    KOKKOS_INLINE_FUNCTION void
    operator()(const member_type &team) const
    {
        const IndexType m = y_.extent(0);
        const IndexType n = x_.extent(0);

        const int j = team.league_rank();
        const IndexType j_min = n_c_ * j;
        const IndexType nj = (j_min + n_c_ > n) ? (n - j_min) : n_c_;
        const IndexType i_min = j % m;

        for (IndexType i = i_min; i < m; ++i)
        {
            Scalar tmp = 0.;
            Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange(team, nj), [=](int jj, Scalar &tmp_sum) {
                    tmp_sum += A_(jj + j_min, i) * x_(jj + j_min);
                },
                tmp);
            if (team.team_rank() == 0)
            {
                tmp *= alpha_;
                Kokkos::atomic_add<Scalar>(&y_(i), tmp);
            }
        }
        for (IndexType i = 0; i < i_min; ++i)
        {
            Scalar tmp = 0.;
            Kokkos::parallel_reduce(
                Kokkos::TeamThreadRange(team, nj), [=](int jj, Scalar &tmp_sum) {
                    tmp_sum += A_(jj + j_min, i) * x_(jj + j_min);
                },
                tmp);
            if (team.team_rank() == 0)
            {
                tmp *= alpha_;
                Kokkos::atomic_add<Scalar>(&y_(i), tmp);
            }
        }
    }

private:
    AlphaCoeffType alpha_;
    typename AViewType::const_type A_;
    typename XViewType::const_type x_;
    YViewType y_;
    const IndexType n_c_;
};

template <class Scalar,
          class VA,
          class VX,
          class VY>
void update_MP(
    typename VA::const_value_type &alpha,
    const VA &A,
    const VX &x,
    typename VY::const_value_type &beta,
    const VY &y)
{
    using execution_space = typename VA::execution_space;
    using IndexType = typename VA::size_type;
    using policy_type = Kokkos::RangePolicy<execution_space, IndexType>;

    // Get the dimensions
    const size_t m = y.extent(0);

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    const size_t N = execution_space::thread_pool_size();
#else
    const size_t N = execution_space::impl_thread_pool_size();
#endif
    const size_t m_c_star = Sacado_MP_Vector_GEMV_Tile_Size(sizeof(Scalar));
    const size_t n_tiles_per_thread = ceil(((double)m) / (N * m_c_star));
    const size_t m_c = ceil(((double)m) / (N * n_tiles_per_thread));
    const size_t n_tiles = N * n_tiles_per_thread;

    policy_type range(0, n_tiles);

    using functor_type = updateF<VA, VX, VY, IndexType>;
    functor_type functor(alpha, A, x, beta, y, m_c);

    Kokkos::parallel_for("KokkosBlas::gemv[Update]", range, functor);
}

template <class Scalar,
          class VA,
          class VX,
          class VY>
void inner_products_MP(
    typename VA::const_value_type &alpha,
    const VA &A,
    const VX &x,
    typename VY::const_value_type &beta,
    const VY &y)
{
    using execution_space = typename VA::execution_space;
    using IndexType = typename VA::size_type;
    using team_policy_type = Kokkos::TeamPolicy<execution_space>;

    // Get the dimensions
    const size_t m = y.extent(0);
    const size_t n = x.extent(0);

    const size_t team_size = STOKHOS_GEMV_TEAM_SIZE;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    const size_t N = execution_space::thread_pool_size();
#else
    const size_t N = execution_space::impl_thread_pool_size();
#endif
    const size_t m_c_star = Sacado_MP_Vector_GEMV_Tile_Size(sizeof(Scalar));
    const size_t n_tiles_per_thread = ceil(((double)n) / (N * m_c_star));
    const size_t m_c = ceil(((double)n) / (N * n_tiles_per_thread));
    const size_t n_per_tile2 = m_c * team_size;

    const size_t n_i2 = ceil(((double)n) / n_per_tile2);

    team_policy_type team(n_i2, team_size);

    if (beta == Scalar(0.))
        Kokkos::parallel_for(
            m, KOKKOS_LAMBDA(const int i) {
                y(i) = beta;
            });
    else
        Kokkos::parallel_for(
            m, KOKKOS_LAMBDA(const int i) {
                y(i) *= beta;
            });

    using functor_type = innerF<VA, VX, VY, IndexType>;
    functor_type functor(alpha, A, x, y, n_per_tile2);

    Kokkos::parallel_for("KokkosBlas::gemv[InnerProducts]", team, functor);
}

namespace KokkosBlas
{
template <typename DA, typename... PA,
          typename DX, typename... PX,
          typename DY, typename... PY>
typename std::enable_if<Kokkos::is_view_mp_vector<Kokkos::View<DA, PA...>>::value &&
                        Kokkos::is_view_mp_vector<Kokkos::View<DX, PX...>>::value &&
                        Kokkos::is_view_mp_vector<Kokkos::View<DY, PY...>>::value>::type
gemv(const char trans[],
     typename Kokkos::View<DA, PA...>::const_value_type &alpha,
     const Kokkos::View<DA, PA...> &A,
     const Kokkos::View<DX, PX...> &x,
     typename Kokkos::View<DY, PY...>::const_value_type &beta,
     const Kokkos::View<DY, PY...> &y)
{
    typedef typename Kokkos::View<DA, PA...>::value_type Scalar;
    typedef Kokkos::View<DA, PA...> VA;
    typedef Kokkos::View<DX, PX...> VX;
    typedef Kokkos::View<DY, PY...> VY;

    static_assert(VA::rank == 2, "GEMM: A must have rank 2 (be a matrix).");
    static_assert(VX::rank == 1, "GEMM: x must have rank 1 (be a vector).");
    static_assert(VY::rank == 1, "GEMM: y must have rank 1 (be a vector).");

    if (trans[0] == 'n' || trans[0] == 'N')
        update_MP<Scalar, VA, VX, VY>(alpha, A, x, beta, y);
    else
        inner_products_MP<Scalar, VA, VX, VY>(alpha, A, x, beta, y);
}

} // namespace KokkosBlas
#endif
