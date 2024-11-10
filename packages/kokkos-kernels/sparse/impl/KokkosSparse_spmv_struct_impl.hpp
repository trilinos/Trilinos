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

#ifndef KOKKOSSPARSE_IMPL_SPMV_STRUCT_DEF_HPP_
#define KOKKOSSPARSE_IMPL_SPMV_STRUCT_DEF_HPP_

#include <sstream>

#include "Kokkos_InnerProductSpaceTraits.hpp"
#include "KokkosKernels_ExecSpaceUtils.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosKernels_Error.hpp"

namespace KokkosSparse {
namespace Impl {

enum { FD, FE };

// This TransposeFunctor is functional, but not necessarily performant.
template <class execution_space, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate>
struct SPMV_Struct_Transpose_Functor {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;
  typedef typename YVector::non_const_value_type coefficient_type;
  typedef typename YVector::non_const_value_type y_value_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  const coefficient_type beta;
  YVector m_y;
  const ordinal_type rows_per_thread;

  SPMV_Struct_Transpose_Functor(const coefficient_type& alpha_, const AMatrix& m_A_, const XVector& m_x_,
                                const coefficient_type& beta_, const YVector& m_y_, const ordinal_type rows_per_thread_)
      : alpha(alpha_), m_A(m_A_), m_x(m_x_), beta(beta_), m_y(m_y_), rows_per_thread(rows_per_thread_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    const ordinal_type teamWorkStart =
        (static_cast<ordinal_type>(dev.league_rank() * dev.team_size() + dev.team_rank())) * rows_per_thread;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, rows_per_thread), [&](ordinal_type loop) {
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.
      ordinal_type iRow = teamWorkStart + loop;
      if (iRow >= m_A.numRows()) {
        return;
      }

      const auto row                = m_A.rowConst(iRow);
      const ordinal_type row_length = row.length;

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, row_length), [&](ordinal_type iEntry) {
        const value_type val   = conjugate ? ATV::conj(row.value(iEntry)) : row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        Kokkos::atomic_add(&m_y(ind), static_cast<y_value_type>(alpha * val * m_x(iRow)));
      });
    });
  }
};

template <class execution_space, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate>
struct SPMV_Struct_Functor {
  typedef typename AMatrix::non_const_size_type size_type;
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename execution_space::scratch_memory_space scratch_space;
  typedef typename KokkosSparse::SparseRowViewConst<AMatrix> row_view_const;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;
  typedef Kokkos::View<ordinal_type*, scratch_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > shared_ordinal_1d;
  using y_value_type = typename YVector::non_const_value_type;

  // Tags to perform SPMV on interior and boundaries
  struct interior3ptTag {};   // 1D FD and FE discretization
  struct interior5ptTag {};   // 2D FD discretization
  struct interior9ptTag {};   // 2D FE discretization
  struct interior7ptTag {};   // 3D FD discretization
  struct interior27ptTag {};  // 3D FE discretization
  struct exterior1DTag {};
  struct exterior2DTag {};
  struct exterior3DTag {};

  // Classic spmv params
  const y_value_type alpha;
  AMatrix m_A;
  XVector m_x;
  const y_value_type beta;
  YVector m_y;

  // Additional structured spmv params
  int numDimensions = 0;
  ordinal_type ni = 0, nj = 0, nk = 0;
  const int stencil_type;
  ordinal_type numInterior, numExterior;
  const int64_t rows_per_team;
  const int64_t rows_per_team_ext;

  SPMV_Struct_Functor(const Kokkos::View<ordinal_type*, Kokkos::HostSpace> structure_, const int stencil_type_,
                      const y_value_type alpha_, const AMatrix m_A_, const XVector m_x_, const y_value_type beta_,
                      const YVector m_y_, const int64_t rows_per_team_, const int64_t rows_per_team_ext_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        beta(beta_),
        m_y(m_y_),
        stencil_type(stencil_type_),
        numExterior(0),
        rows_per_team(rows_per_team_),
        rows_per_team_ext(rows_per_team_ext_) {
    static_assert(static_cast<int>(XVector::rank) == 1, "XVector must be a rank 1 View.");
    static_assert(static_cast<int>(YVector::rank) == 1, "YVector must be a rank 1 View.");

    numDimensions = structure_.extent(0);
    if (numDimensions == 1) {
      ni = static_cast<ordinal_type>(structure_(0));
    } else if (numDimensions == 2) {
      ni = static_cast<ordinal_type>(structure_(0));
      nj = static_cast<ordinal_type>(structure_(1));
    } else if (numDimensions == 3) {
      ni = static_cast<ordinal_type>(structure_(0));
      nj = static_cast<ordinal_type>(structure_(1));
      nk = static_cast<ordinal_type>(structure_(2));
    }
  }

  void compute_interior(const execution_space& exec, const int64_t worksets, const int team_size,
                        const int vector_length) {
    if (numDimensions == 1) {
      // Treat interior points using structured algorithm
      numInterior = ni - 2;
      if (numInterior > 0) {
        size_t shared_size = shared_ordinal_1d::shmem_size(3);
        Kokkos::TeamPolicy<interior3ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> > policy(exec, 1, 1);
        if (team_size < 0) {
          policy = Kokkos::TeamPolicy<interior3ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                       exec, worksets, Kokkos::AUTO, vector_length)
                       .set_scratch_size(0, Kokkos::PerTeam(shared_size));
        } else {
          policy = Kokkos::TeamPolicy<interior3ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                       exec, worksets, team_size, vector_length)
                       .set_scratch_size(0, Kokkos::PerTeam(shared_size));
        }
        Kokkos::parallel_for("KokkosSparse::spmv_struct<NoTranspose,Static>: interior", policy, *this);
      }

    } else if (numDimensions == 2) {
      // Treat interior points using structured algorithm
      numInterior = (ni - 2) * (nj - 2);
      if (numInterior > 0) {
        if (stencil_type == 1) {
          size_t shared_size = shared_ordinal_1d::shmem_size(5);
          Kokkos::TeamPolicy<interior5ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> > policy(exec, 1, 1);
          if (team_size < 0) {
            policy = Kokkos::TeamPolicy<interior5ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                         exec, worksets, Kokkos::AUTO, vector_length)
                         .set_scratch_size(0, Kokkos::PerTeam(shared_size));
          } else {
            policy = Kokkos::TeamPolicy<interior5ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                         exec, worksets, team_size, vector_length)
                         .set_scratch_size(0, Kokkos::PerTeam(shared_size));
          }
          Kokkos::parallel_for("KokkosSparse::spmv_struct<NoTranspose,Static>: interior", policy, *this);
        } else if (stencil_type == 2) {
          size_t shared_size = shared_ordinal_1d::shmem_size(9);
          Kokkos::TeamPolicy<interior9ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> > policy(exec, 1, 1);
          if (team_size < 0) {
            policy = Kokkos::TeamPolicy<interior9ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                         exec, worksets, Kokkos::AUTO, vector_length)
                         .set_scratch_size(0, Kokkos::PerTeam(shared_size));
          } else {
            policy = Kokkos::TeamPolicy<interior9ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                         exec, worksets, team_size, vector_length)
                         .set_scratch_size(0, Kokkos::PerTeam(shared_size));
          }
          Kokkos::parallel_for("KokkosSparse::spmv_struct<NoTranspose,Static>: interior", policy, *this);
        }
      }

    } else if (numDimensions == 3) {
      // Treat interior points using structured algorithm
      numInterior = (ni - 2) * (nj - 2) * (nk - 2);
      if (numInterior > 0) {
        if (stencil_type == 1) {
          size_t shared_size = shared_ordinal_1d::shmem_size(7);
          Kokkos::TeamPolicy<interior7ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> > policy(exec, 1, 1);
          if (team_size < 0) {
            policy = Kokkos::TeamPolicy<interior7ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                         exec, worksets, Kokkos::AUTO, vector_length)
                         .set_scratch_size(0, Kokkos::PerTeam(shared_size));
          } else {
            policy = Kokkos::TeamPolicy<interior7ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                         exec, worksets, team_size, vector_length)
                         .set_scratch_size(0, Kokkos::PerTeam(shared_size));
          }
          Kokkos::parallel_for("KokkosSparse::spmv_struct<NoTranspose,Static>: interior", policy, *this);
        } else if (stencil_type == 2) {
          size_t shared_size = shared_ordinal_1d::shmem_size(27);
          Kokkos::TeamPolicy<interior27ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> > policy(exec, 1, 1);
          if (team_size < 0) {
            policy = Kokkos::TeamPolicy<interior27ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                         exec, worksets, Kokkos::AUTO, vector_length)
                         .set_scratch_size(0, Kokkos::PerTeam(shared_size));
          } else {
            policy = Kokkos::TeamPolicy<interior27ptTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
                         exec, worksets, team_size, vector_length)
                         .set_scratch_size(0, Kokkos::PerTeam(shared_size));
          }
          Kokkos::parallel_for("KokkosSparse::spmv_struct<NoTranspose,Static>: interior", policy, *this);
        }
      }
    }
  }  // compute_interior

  KOKKOS_INLINE_FUNCTION
  void operator()(const interior3ptTag&, const team_member& dev) const {
    // Allocate and initialize columnOffsets array for the team
    shared_ordinal_1d columnOffsets(dev.team_scratch(0), 3);
    Kokkos::single(Kokkos::PerTeam(dev), [&]() {
      columnOffsets(0) = -1;
      columnOffsets(1) = 0;
      columnOffsets(2) = 1;
    });
    dev.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, rows_per_team), [&](const ordinal_type& loop) {
      const ordinal_type interiorIdx = static_cast<ordinal_type>(dev.league_rank()) * rows_per_team + loop;
      if (interiorIdx >= numInterior) {
        return;
      }

      ordinal_type rowIdx;
      rowIdx = interiorIdx + 1;

      const size_type rowOffset = m_A.graph.row_map(rowIdx);
      y_value_type sum          = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(dev, 3),
          [&](const ordinal_type& idx, y_value_type& lclSum) {
            lclSum += (conjugate ? ATV::conj(m_A.values(rowOffset + idx)) : m_A.values(rowOffset + idx)) *
                      m_x(rowIdx + columnOffsets(idx));
          },
          sum);

      Kokkos::single(Kokkos::PerThread(dev), [&]() { m_y(rowIdx) = beta * m_y(rowIdx) + alpha * sum; });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const interior5ptTag&, const team_member& dev) const {
    // Allocate and initialize columnOffsets array for the team
    shared_ordinal_1d columnOffsets(dev.team_scratch(0), 5);
    Kokkos::single(Kokkos::PerTeam(dev), [&]() {
      columnOffsets(0) = -ni;
      columnOffsets(1) = -1;
      columnOffsets(2) = 0;
      columnOffsets(3) = 1;
      columnOffsets(4) = ni;
    });
    dev.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, rows_per_team), [&](const ordinal_type& loop) {
      const ordinal_type interiorIdx = static_cast<ordinal_type>(dev.league_rank()) * rows_per_team + loop;
      if (interiorIdx >= numInterior) {
        return;
      }

      ordinal_type i, j, rowIdx;
      j      = interiorIdx / (ni - 2);
      i      = interiorIdx % (ni - 2);
      rowIdx = (j + 1) * ni + i + 1;

      const size_type rowOffset = m_A.graph.row_map(rowIdx);
      y_value_type sum          = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(dev, 5),
          [&](const ordinal_type& idx, y_value_type& lclSum) {
            lclSum += (conjugate ? ATV::conj(m_A.values(rowOffset + idx)) : m_A.values(rowOffset + idx)) *
                      m_x(rowIdx + columnOffsets(idx));
          },
          sum);

      Kokkos::single(Kokkos::PerThread(dev), [&]() { m_y(rowIdx) = beta * m_y(rowIdx) + alpha * sum; });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const interior9ptTag&, const team_member& dev) const {
    // Allocate and initialize columnOffsets array for the team
    shared_ordinal_1d columnOffsets(dev.team_scratch(0), 9);
    Kokkos::single(Kokkos::PerTeam(dev), [&]() {
      columnOffsets(0) = -ni - 1;
      columnOffsets(1) = -ni;
      columnOffsets(2) = -ni + 1;
      columnOffsets(3) = -1;
      columnOffsets(4) = 0;
      columnOffsets(5) = 1;
      columnOffsets(6) = ni - 1;
      columnOffsets(7) = ni;
      columnOffsets(8) = ni + 1;
    });
    dev.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, rows_per_team), [&](const ordinal_type& loop) {
      const ordinal_type interiorIdx = static_cast<ordinal_type>(dev.league_rank()) * rows_per_team + loop;
      if (interiorIdx >= numInterior) {
        return;
      }

      ordinal_type i, j, rowIdx;
      j      = interiorIdx / (ni - 2);
      i      = interiorIdx % (ni - 2);
      rowIdx = (j + 1) * ni + i + 1;

      const size_type rowOffset = m_A.graph.row_map(rowIdx);
      y_value_type sum          = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(dev, 9),
          [&](const ordinal_type& idx, y_value_type& lclSum) {
            lclSum += (conjugate ? ATV::conj(m_A.values(rowOffset + idx)) : m_A.values(rowOffset + idx)) *
                      m_x(rowIdx + columnOffsets(idx));
          },
          sum);

      Kokkos::single(Kokkos::PerThread(dev), [&]() { m_y(rowIdx) = beta * m_y(rowIdx) + alpha * sum; });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const interior7ptTag&, const team_member& dev) const {
    // Allocate and initialize columnOffsets array for the team
    shared_ordinal_1d columnOffsets(dev.team_scratch(0), 7);
    Kokkos::single(Kokkos::PerTeam(dev), [&]() {
      columnOffsets(0) = -ni * nj;
      columnOffsets(1) = -ni;
      columnOffsets(2) = -1;
      columnOffsets(3) = 0;
      columnOffsets(4) = 1;
      columnOffsets(5) = ni;
      columnOffsets(6) = ni * nj;
    });
    dev.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, rows_per_team), [&](const ordinal_type& loop) {
      const ordinal_type interiorIdx = static_cast<ordinal_type>(dev.league_rank()) * rows_per_team + loop;
      if (interiorIdx >= numInterior) {
        return;
      }

      ordinal_type i, j, k, rowIdx, rem;
      k      = interiorIdx / ((ni - 2) * (nj - 2));
      rem    = interiorIdx % ((ni - 2) * (nj - 2));
      j      = rem / (ni - 2);
      i      = rem % (ni - 2);
      rowIdx = (k + 1) * nj * ni + (j + 1) * ni + (i + 1);

      const size_type rowOffset = m_A.graph.row_map(rowIdx);
      y_value_type sum          = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(dev, 7),
          [&](const ordinal_type& idx, y_value_type& lclSum) {
            lclSum += (conjugate ? ATV::conj(m_A.values(rowOffset + idx)) : m_A.values(rowOffset + idx)) *
                      m_x(rowIdx + columnOffsets(idx));
          },
          sum);

      Kokkos::single(Kokkos::PerThread(dev), [&]() { m_y(rowIdx) = beta * m_y(rowIdx) + alpha * sum; });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const interior27ptTag&, const team_member& dev) const {
    // Allocate and initialize columnOffsets array for the team
    shared_ordinal_1d columnOffsets(dev.team_scratch(0), 27);
    Kokkos::single(Kokkos::PerTeam(dev), [&]() {
      columnOffsets(0)  = -ni * nj - ni - 1;
      columnOffsets(1)  = -ni * nj - ni;
      columnOffsets(2)  = -ni * nj - ni + 1;
      columnOffsets(3)  = -ni * nj - 1;
      columnOffsets(4)  = -ni * nj;
      columnOffsets(5)  = -ni * nj + 1;
      columnOffsets(6)  = -ni * nj + ni - 1;
      columnOffsets(7)  = -ni * nj + ni;
      columnOffsets(8)  = -ni * nj + ni + 1;
      columnOffsets(9)  = -ni - 1;
      columnOffsets(10) = -ni;
      columnOffsets(11) = -ni + 1;
      columnOffsets(12) = -1;
      columnOffsets(13) = 0;
      columnOffsets(14) = 1;
      columnOffsets(15) = ni - 1;
      columnOffsets(16) = ni;
      columnOffsets(17) = ni + 1;
      columnOffsets(18) = ni * nj - ni - 1;
      columnOffsets(19) = ni * nj - ni;
      columnOffsets(20) = ni * nj - ni + 1;
      columnOffsets(21) = ni * nj - 1;
      columnOffsets(22) = ni * nj;
      columnOffsets(23) = ni * nj + 1;
      columnOffsets(24) = ni * nj + ni - 1;
      columnOffsets(25) = ni * nj + ni;
      columnOffsets(26) = ni * nj + ni + 1;
    });
    dev.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, rows_per_team), [&](const ordinal_type& loop) {
      const ordinal_type interiorIdx = static_cast<ordinal_type>(dev.league_rank()) * rows_per_team + loop;
      if (interiorIdx >= numInterior) {
        return;
      }

      ordinal_type i, j, k, rowIdx, rem;
      k      = interiorIdx / ((ni - 2) * (nj - 2));
      rem    = interiorIdx % ((ni - 2) * (nj - 2));
      j      = rem / (ni - 2);
      i      = rem % (ni - 2);
      rowIdx = (k + 1) * nj * ni + (j + 1) * ni + (i + 1);

      const size_type rowOffset = m_A.graph.row_map(rowIdx);

      y_value_type sum(0.0);

      // clang-format off
          KOKKOS_IF_ON_HOST((
          for (ordinal_type idx = 0; idx < 27; ++idx) {
            sum +=
                m_A.values(rowOffset + idx) * m_x(rowIdx + columnOffsets(idx));
          }
          ))

          KOKKOS_IF_ON_DEVICE((
          Kokkos::parallel_reduce(
              Kokkos::ThreadVectorRange(dev, 27),
              [&](const ordinal_type& idx, y_value_type& lclSum) {
                lclSum += (conjugate ? ATV::conj(m_A.values(rowOffset + idx))
                                     : m_A.values(rowOffset + idx)) *
                          m_x(rowIdx + columnOffsets(idx));
              },
              sum);
          ))
      // clang-format on

      Kokkos::single(Kokkos::PerThread(dev), [&]() { m_y(rowIdx) = beta * m_y(rowIdx) + alpha * sum; });
    });
  }

  void compute_exterior(const execution_space& exec, const int64_t worksets, const int team_size,
                        const int vector_length) {
    // Treat exterior points using unstructured algorithm
    if (numDimensions == 1) {
      numExterior = 2;
      if (numExterior > 0) {
        Kokkos::RangePolicy<exterior1DTag, execution_space, Kokkos::Schedule<Kokkos::Static> > policy(exec, 0,
                                                                                                      numExterior);
        Kokkos::parallel_for("KokkosSparse::spmv_struct<NoTranspose,Static>: exterior", policy, *this);
      }

    } else if (numDimensions == 2) {
      numExterior = 2 * (nj + ni - 2);
      if (numExterior > 0) {
        Kokkos::TeamPolicy<exterior2DTag, execution_space, Kokkos::Schedule<Kokkos::Static> > policy(exec, 1, 1);
        if (team_size < 0) {
          policy = Kokkos::TeamPolicy<exterior2DTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
              exec, worksets, Kokkos::AUTO, vector_length);
        } else {
          policy = Kokkos::TeamPolicy<exterior2DTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
              exec, worksets, team_size, vector_length);
        }
        Kokkos::parallel_for("KokkosSparse::spmv_struct<NoTranspose,Static>: exterior", policy, *this);
      }

    } else if (numDimensions == 3) {
      numExterior = ni * nj * nk - (ni - 2) * (nj - 2) * (nk - 2);
      if (numExterior > 0) {
        Kokkos::TeamPolicy<exterior3DTag, execution_space, Kokkos::Schedule<Kokkos::Static> > policy(exec, 1, 1);
        if (team_size < 0) {
          policy = Kokkos::TeamPolicy<exterior3DTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
              exec, worksets, Kokkos::AUTO, vector_length);
        } else {
          policy = Kokkos::TeamPolicy<exterior3DTag, execution_space, Kokkos::Schedule<Kokkos::Static> >(
              exec, worksets, team_size, vector_length);
        }

        Kokkos::parallel_for("KokkosSparse::spmv_struct<NoTranspose,Static>: exterior", policy, *this);
      }
    }
  }  // compute_exterior

  KOKKOS_INLINE_FUNCTION
  void operator()(const exterior1DTag&, const ordinal_type& exteriorIdx) const {
    if (exteriorIdx >= numExterior) {
      return;
    }
    ordinal_type rowIdx = exteriorIdx * (ni - 1);

    const size_type rowOffset     = m_A.graph.row_map(rowIdx);
    const ordinal_type row_length = static_cast<ordinal_type>(m_A.graph.row_map(rowIdx + 1) - rowOffset);
    y_value_type sum              = 0.0;
    for (ordinal_type idx = 0; idx < row_length; ++idx) {
      sum += (conjugate ? ATV::conj(m_A.values(rowOffset + idx)) : m_A.values(rowOffset + idx)) *
             m_x(m_A.graph.entries(rowOffset + idx));
    }
    m_y(rowIdx) = beta * m_y(rowIdx) + alpha * sum;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const exterior2DTag&, const team_member& dev) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, rows_per_team_ext), [&](const ordinal_type& loop) {
      const ordinal_type exteriorIdx = static_cast<ordinal_type>(dev.league_rank()) * rows_per_team_ext + loop;
      if (exteriorIdx >= numExterior) {
        return;
      }
      const ordinal_type topFlag    = exteriorIdx / (ni + 2 * nj - 4);
      const ordinal_type bottomFlag = static_cast<ordinal_type>((exteriorIdx / ni) == 0);

      ordinal_type rowIdx = 0;
      if (bottomFlag == 1) {
        rowIdx = exteriorIdx;
      } else if (topFlag == 1) {
        rowIdx = exteriorIdx - (ni + 2 * nj - 4) + ni * (nj - 1);
      } else {
        ordinal_type edgeIdx = (exteriorIdx - ni) / 2;
        ordinal_type edgeFlg = (exteriorIdx - ni) % 2;
        rowIdx               = (edgeIdx + 1) * ni + edgeFlg * (ni - 1);
      }

      const size_type rowOffset     = m_A.graph.row_map(rowIdx);
      const ordinal_type row_length = static_cast<ordinal_type>(m_A.graph.row_map(rowIdx + 1) - rowOffset);
      y_value_type sum              = 0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(dev, row_length),
          [&](const ordinal_type& idx, y_value_type& lclSum) {
            lclSum += (conjugate ? ATV::conj(m_A.values(rowOffset + idx)) : m_A.values(rowOffset + idx)) *
                      m_x(m_A.graph.entries(rowOffset + idx));
          },
          sum);
      Kokkos::single(Kokkos::PerThread(dev), [&]() { m_y(rowIdx) = beta * m_y(rowIdx) + alpha * sum; });
    });
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const exterior3DTag&, const team_member& dev) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, 0, rows_per_team_ext), [&](const ordinal_type& loop) {
      const ordinal_type exteriorIdx = static_cast<ordinal_type>(dev.league_rank()) * rows_per_team_ext + loop;
      if (exteriorIdx >= numExterior) {
        return;
      }
      const ordinal_type topFlag    = static_cast<ordinal_type>(numExterior - exteriorIdx - 1 < ni * nj);
      const ordinal_type bottomFlag = static_cast<ordinal_type>(exteriorIdx / (ni * nj) == 0);

      ordinal_type rowIdx = 0;
      if (bottomFlag == 1) {
        rowIdx = exteriorIdx;
      } else if (topFlag == 1) {
        rowIdx = exteriorIdx - ni * nj - 2 * (nk - 2) * (nj + ni - 2) + (nk - 1) * ni * nj;
      } else {
        ordinal_type k, rem;
        k   = (exteriorIdx - ni * nj) / (2 * (ni - 1 + nj - 1));
        rem = (exteriorIdx - ni * nj) % (2 * (ni - 1 + nj - 1));
        // ordinal_type frontFlg = static_cast<ordinal_type>(rem < ni);
        // ordinal_type backFlg = static_cast<ordinal_type>(rem - ni - 2*(nj
        // - 1) - 1 > 0);
        if (rem < ni) {
          rowIdx = (k + 1) * ni * nj + rem;
        } else if (rem < ni + 2 * (nj - 2)) {
          ordinal_type edgeIdx = (rem - ni) / 2;
          ordinal_type edgeFlg = (rem - ni) % 2;
          if (edgeFlg == 0) {
            rowIdx = (k + 1) * ni * nj + (edgeIdx + 1) * ni;
          } else if (edgeFlg == 1) {
            rowIdx = (k + 1) * ni * nj + (edgeIdx + 2) * ni - 1;
          }
        } else {
          rowIdx = (k + 1) * ni * nj + rem - ni - 2 * (nj - 2) + (nj - 1) * ni;
        }
      }

      const size_type rowOffset     = m_A.graph.row_map(rowIdx);
      const ordinal_type row_length = static_cast<ordinal_type>(m_A.graph.row_map(rowIdx + 1) - rowOffset);
      y_value_type sum              = 0;
      Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(dev, row_length),
          [&](const ordinal_type& idx, y_value_type& lclSum) {
            lclSum += (conjugate ? ATV::conj(m_A.values(rowOffset + idx)) : m_A.values(rowOffset + idx)) *
                      m_x(m_A.graph.entries(rowOffset + idx));
          },
          sum);
      Kokkos::single(Kokkos::PerThread(dev), [&]() { m_y(rowIdx) = beta * m_y(rowIdx) + alpha * sum; });
    });
  }
};  // SPMV_Struct_Functor

template <class execution_space>
int64_t spmv_struct_launch_parameters(int64_t numInterior, int64_t nnz, int nnz_per_row, int64_t rows_per_thread,
                                      int& team_size, int& vector_length) {
  int64_t rows_per_team;

  if (nnz_per_row < 1) nnz_per_row = 1;

  // Determine rows per thread
  if (rows_per_thread < 1) {
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>())
      rows_per_thread = 1;
    else {
      if (nnz_per_row < 20 && numInterior * nnz_per_row > 5000000) {
        rows_per_thread = 256;
      } else
        rows_per_thread = 64;
    }
  }

  if (team_size < 1) {
    if (KokkosKernels::Impl::kk_is_gpu_exec_space<execution_space>()) {
      team_size = 128 / vector_length;
    } else {
      team_size = 1;
    }
  }

  rows_per_team = rows_per_thread * team_size;

  if (rows_per_team < 0) {
    int64_t nnz_per_team = 4096;
    int64_t conc         = execution_space().concurrency();
    while ((conc * nnz_per_team * 4 > nnz) && (nnz_per_team > 256)) nnz_per_team /= 2;
    rows_per_team = (nnz_per_team + nnz_per_row - 1) / nnz_per_row;
  }

  return rows_per_team;
}  // spmv_struct_launch_parameters

template <class execution_space, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate>
static void spmv_struct_beta_no_transpose(
    const execution_space& exec, const int stencil_type,
    const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
    typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
    typename YVector::const_value_type& beta, const YVector& y) {
  typedef typename AMatrix::ordinal_type ordinal_type;
  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  int vector_length           = -1;
  int nnzPerRow               = -1;
  int team_size_int           = -1;
  int team_size_ext           = -1;
  int64_t rows_per_thread_int = -1;
  int64_t rows_per_thread_ext = -1;
  int64_t numInteriorPts      = 0;
  int64_t numExteriorPts      = 0;

  if (structure.extent(0) == 1) {
    numInteriorPts = structure(0) - 2;
    numExteriorPts = 2;
    vector_length  = 1;
  } else if (structure.extent(0) == 2) {
    numInteriorPts = (structure(1) - 2) * (structure(0) - 2);
    numExteriorPts = 2 * (structure(1) + structure(0) - 2);
    if (stencil_type == 1) {
      vector_length = 2;
    } else if (stencil_type == 2) {
      vector_length = 2;
    }
  } else if (structure.extent(0) == 3) {
    numInteriorPts = (structure(2) - 2) * (structure(1) - 2) * (structure(0) - 2);
    numExteriorPts = structure(2) * structure(1) * structure(0) - numInteriorPts;
    if (stencil_type == 1) {
      vector_length = 2;
    } else if (stencil_type == 2) {
      vector_length = 8;
    }
  }

  int64_t rows_per_team_int = spmv_struct_launch_parameters<execution_space>(
      numInteriorPts, A.nnz(), nnzPerRow, rows_per_thread_int, team_size_int, vector_length);
  int64_t worksets_interior = (numInteriorPts + rows_per_team_int - 1) / rows_per_team_int;

  int64_t rows_per_team_ext = spmv_struct_launch_parameters<execution_space>(
      numExteriorPts, A.nnz(), nnzPerRow, rows_per_thread_ext, team_size_ext, vector_length);
  int64_t worksets_exterior = (numExteriorPts + rows_per_team_ext - 1) / rows_per_team_ext;

  SPMV_Struct_Functor<execution_space, AMatrix, XVector, YVector, dobeta, conjugate> spmv_struct(
      structure, stencil_type, alpha, A, x, beta, y, rows_per_team_int, rows_per_team_ext);

  spmv_struct.compute_interior(exec, worksets_interior, team_size_int, vector_length);
  spmv_struct.compute_exterior(exec, worksets_exterior, team_size_ext, vector_length);
}  // spmv_struct_beta_no_transpose

template <class execution_space, class AMatrix, class XVector, class YVector, int dobeta, bool conjugate>
static void spmv_struct_beta_transpose(
    const execution_space& exec, const int /*stencil_type*/,
    const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& /*structure*/,
    typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
    typename YVector::const_value_type& beta, const YVector& y) {
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  if (dobeta != 1) {
    KokkosBlas::scal(exec, y, beta, y);
  }

  typedef typename AMatrix::size_type size_type;

  // Assuming that no row contains duplicate entries, NNZPerRow
  // cannot be more than the number of columns of the matrix.  Thus,
  // the appropriate type is ordinal_type.
  const ordinal_type NNZPerRow = static_cast<ordinal_type>(A.nnz() / A.numRows());

  int vector_length = 1;
  while ((static_cast<ordinal_type>(vector_length * 2 * 3) <= NNZPerRow) && (vector_length < 32)) vector_length *= 2;

  typedef SPMV_Struct_Transpose_Functor<execution_space, AMatrix, XVector, YVector, dobeta, conjugate> OpType;

  typename AMatrix::const_ordinal_type nrow = A.numRows();

  OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow));

  const int rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
  const int team_size       = Kokkos::TeamPolicy<execution_space>(rows_per_thread, Kokkos::AUTO, vector_length)
                            .team_size_recommended(op, Kokkos::ParallelForTag());
  const int rows_per_team = rows_per_thread * team_size;
  const size_type nteams  = (nrow + rows_per_team - 1) / rows_per_team;
  Kokkos::parallel_for("KokkosSparse::spmv_struct<Transpose>",
                       Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);
}

template <class execution_space, class AMatrix, class XVector, class YVector, int dobeta>
static void spmv_struct_beta(
    const execution_space& exec, const char mode[], const int stencil_type,
    const Kokkos::View<typename AMatrix::non_const_ordinal_type*, Kokkos::HostSpace>& structure,
    typename YVector::const_value_type& alpha, const AMatrix& A, const XVector& x,
    typename YVector::const_value_type& beta, const YVector& y) {
  if (mode[0] == NoTranspose[0]) {
    spmv_struct_beta_no_transpose<execution_space, AMatrix, XVector, YVector, dobeta, false>(
        exec, stencil_type, structure, alpha, A, x, beta, y);
  } else if (mode[0] == Conjugate[0]) {
    spmv_struct_beta_no_transpose<execution_space, AMatrix, XVector, YVector, dobeta, true>(
        exec, stencil_type, structure, alpha, A, x, beta, y);
  } else if (mode[0] == Transpose[0]) {
    spmv_struct_beta_transpose<execution_space, AMatrix, XVector, YVector, dobeta, false>(exec, stencil_type, structure,
                                                                                          alpha, A, x, beta, y);
  } else if (mode[0] == ConjugateTranspose[0]) {
    spmv_struct_beta_transpose<execution_space, AMatrix, XVector, YVector, dobeta, true>(exec, stencil_type, structure,
                                                                                         alpha, A, x, beta, y);
  } else {
    std::stringstream ss;
    ss << __FILE__ << ":" << __LINE__ << " Invalid transpose mode " << mode << " for KokkosSparse::spmv_struct()";
    KokkosKernels::Impl::throw_runtime_exception(ss.str());
  }
}

// Functor for implementing transpose and conjugate transpose sparse
// matrix-vector multiply with multivector (2-D View) input and
// output.  This functor works, but is not necessarily performant.
template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate>
struct SPMV_MV_Struct_Transpose_Functor {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type A_value_type;
  typedef typename YVector::non_const_value_type y_value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef typename YVector::non_const_value_type coefficient_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  const coefficient_type beta;
  YVector m_y;

  const ordinal_type n;
  const ordinal_type rows_per_thread;

  SPMV_MV_Struct_Transpose_Functor(const coefficient_type& alpha_, const AMatrix& m_A_, const XVector& m_x_,
                                   const coefficient_type& beta_, const YVector& m_y_,
                                   const ordinal_type rows_per_thread_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        beta(beta_),
        m_y(m_y_),
        n(m_x_.extent(1)),
        rows_per_thread(rows_per_thread_) {}

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    // This should be a thread loop as soon as we can use C++11
    const ordinal_type teamWorkStart =
        (static_cast<ordinal_type>(dev.league_rank() * dev.team_size() + dev.team_rank())) * rows_per_thread;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev, rows_per_thread), [&](ordinal_type loop) {
      const ordinal_type iRow = teamWorkStart + loop;
      if (iRow >= m_A.numRows()) {
        return;
      }
      const auto row                = m_A.rowConst(iRow);
      const ordinal_type row_length = row.length;

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, row_length), [&](ordinal_type iEntry) {
        const A_value_type val =
            conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type k = 0; k < n; ++k) {
            Kokkos::atomic_add(&m_y(ind, k), static_cast<y_value_type>(alpha * val * m_x(iRow, k)));
          }
        } else {
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
          for (ordinal_type k = 0; k < n; ++k) {
            Kokkos::atomic_add(&m_y(ind, k), static_cast<y_value_type>(val * m_x(iRow, k)));
          }
        }
      });
    });
  }
};

template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate>
struct SPMV_MV_Struct_LayoutLeft_Functor {
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type A_value_type;
  typedef typename YVector::non_const_value_type y_value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef typename YVector::non_const_value_type coefficient_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  const coefficient_type beta;
  YVector m_y;
  //! The number of columns in the input and output MultiVectors.
  ordinal_type n;
  ordinal_type rows_per_thread;
  int vector_length;

  SPMV_MV_Struct_LayoutLeft_Functor(const coefficient_type& alpha_, const AMatrix& m_A_, const XVector& m_x_,
                                    const coefficient_type& beta_, const YVector& m_y_,
                                    const ordinal_type rows_per_thread_, int vector_length_)
      : alpha(alpha_),
        m_A(m_A_),
        m_x(m_x_),
        beta(beta_),
        m_y(m_y_),
        n(m_x_.extent(1)),
        rows_per_thread(rows_per_thread_),
        vector_length(vector_length_) {}

  template <int UNROLL>
  KOKKOS_INLINE_FUNCTION void strip_mine(const team_member& dev, const ordinal_type& iRow,
                                         const ordinal_type& kk) const {
    y_value_type sum[UNROLL];

#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      sum[k] = Kokkos::ArithTraits<y_value_type>::zero();
    }

    const auto row = m_A.rowConst(iRow);

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, row.length), [&](ordinal_type iEntry) {
      const A_value_type val =
          conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);

#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        sum[k] += val * m_x(ind, kk + k);
      }
    });

    if (doalpha == -1) {
      for (int ii = 0; ii < UNROLL; ++ii) {
        y_value_type sumt;
        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(dev, vector_length),
            [&](ordinal_type, y_value_type& lsum) {
              // in this context, sum[ii] is a partial sum ii on one of the
              // vector lanes.
              lsum -= sum[ii];
            },
            sumt);
        sum[ii] = sumt;
        // that was an all-reduce, so sum[ii] is the same on every vector lane
      }
    } else {
      for (int ii = 0; ii < UNROLL; ++ii) {
        y_value_type sumt;
        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(dev, vector_length),
            [&](ordinal_type, y_value_type& lsum) {
              // in this context, sum[ii] is a partial sum ii on one of the
              // vector lanes.
              lsum += sum[ii];
            },
            sumt);
        if (doalpha == 1)
          sum[ii] = sumt;
        else
          sum[ii] = sumt * alpha;
      }
    }

    Kokkos::single(Kokkos::PerThread(dev), [&]() {
      if (dobeta == 0) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, UNROLL),
                             [&](ordinal_type k) { m_y(iRow, kk + k) = sum[k]; });
      } else if (dobeta == 1) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, UNROLL),
                             [&](ordinal_type k) { m_y(iRow, kk + k) += sum[k]; });
      } else if (dobeta == -1) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, UNROLL),
                             [&](ordinal_type k) { m_y(iRow, kk + k) = -m_y(iRow, kk + k) + sum[k]; });
      } else {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(dev, UNROLL),
                             [&](ordinal_type k) { m_y(iRow, kk + k) = beta * m_y(iRow, kk + k) + sum[k]; });
      }
    });
  }

  KOKKOS_INLINE_FUNCTION void strip_mine_1(const team_member& dev, const ordinal_type& iRow) const {
    const auto row = m_A.rowConst(iRow);

    y_value_type sum;
    Kokkos::parallel_reduce(
        Kokkos::ThreadVectorRange(dev, row.length),
        [&](ordinal_type iEntry, y_value_type& lsum) {
          const A_value_type val =
              conjugate ? Kokkos::ArithTraits<A_value_type>::conj(row.value(iEntry)) : row.value(iEntry);
          lsum += val * m_x(row.colidx(iEntry), 0);
        },
        sum);

    Kokkos::single(Kokkos::PerThread(dev), [&]() {
      if (doalpha == -1) {
        sum = -sum;
      } else if (doalpha * doalpha != 1) {
        sum *= alpha;
      }

      if (dobeta == 0) {
        m_y(iRow, 0) = sum;
      } else if (dobeta == 1) {
        m_y(iRow, 0) += sum;
      } else if (dobeta == -1) {
        m_y(iRow, 0) = -m_y(iRow, 0) + sum;
      } else {
        m_y(iRow, 0) = beta * m_y(iRow, 0) + sum;
      }
    });
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    for (ordinal_type loop = 0; loop < rows_per_thread; ++loop) {
      // iRow indexes over (local) rows of the matrix, so its correct
      // type is ordinal_type.

      const ordinal_type iRow = (dev.league_rank() * dev.team_size() + dev.team_rank()) * rows_per_thread + loop;
      if (iRow >= m_A.numRows()) {
        return;
      }

      // mfh 20 Mar 2015, 07 Jun 2016: This is ordinal_type because it
      // needs to have the same type as n.
      ordinal_type kk = 0;

#ifdef KOKKOS_FAST_COMPILE
      for (; kk + 4 <= n; kk += 4) {
        strip_mine<4>(dev, iRow, kk);
      }
      for (; kk < n; ++kk) {
        strip_mine<1>(dev, iRow, kk);
      }
#else
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
      if ((n > 8) && (n % 8 == 1)) {
        strip_mine<9>(dev, iRow, kk);
        kk += 9;
      }
      for (; kk + 8 <= n; kk += 8) strip_mine<8>(dev, iRow, kk);
      if (kk < n) {
        switch (n - kk) {
#else   // NOT a CUDA device
      if ((n > 16) && (n % 16 == 1)) {
        strip_mine<17>(dev, iRow, kk);
        kk += 17;
      }

      for (; kk + 16 <= n; kk += 16) {
        strip_mine<16>(dev, iRow, kk);
      }

      if (kk < n) {
        switch (n - kk) {
          case 15: strip_mine<15>(dev, iRow, kk); break;

          case 14: strip_mine<14>(dev, iRow, kk); break;

          case 13: strip_mine<13>(dev, iRow, kk); break;

          case 12: strip_mine<12>(dev, iRow, kk); break;

          case 11: strip_mine<11>(dev, iRow, kk); break;

          case 10: strip_mine<10>(dev, iRow, kk); break;

          case 9: strip_mine<9>(dev, iRow, kk); break;

          case 8: strip_mine<8>(dev, iRow, kk); break;
#endif  // __CUDA_ARCH__ or __HIP_DEVICE_COMPILE__
          case 7: strip_mine<7>(dev, iRow, kk); break;

          case 6: strip_mine<6>(dev, iRow, kk); break;

          case 5: strip_mine<5>(dev, iRow, kk); break;

          case 4: strip_mine<4>(dev, iRow, kk); break;

          case 3: strip_mine<3>(dev, iRow, kk); break;

          case 2: strip_mine<2>(dev, iRow, kk); break;

          case 1: strip_mine_1(dev, iRow); break;
        }
      }
#endif  // KOKKOS_FAST_COMPILE
    }
  }
};

template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate>
static void spmv_alpha_beta_mv_struct_no_transpose(const execution_space& exec,
                                                   const typename YVector::non_const_value_type& alpha,
                                                   const AMatrix& A, const XVector& x,
                                                   const typename YVector::non_const_value_type& beta,
                                                   const YVector& y) {
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }
  if (doalpha == 0) {
    if (dobeta != 1) {
      KokkosBlas::scal(exec, y, beta, y);
    }
    return;
  } else {
    typedef typename AMatrix::size_type size_type;

    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type>(A.nnz() / A.numRows());

    int vector_length = 1;
    while ((static_cast<ordinal_type>(vector_length * 2 * 3) <= NNZPerRow) && (vector_length < 8)) vector_length *= 2;

#ifndef KOKKOS_FAST_COMPILE  // This uses templated functions on doalpha and
                             // dobeta and will produce 16 kernels

    typedef SPMV_MV_Struct_LayoutLeft_Functor<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, conjugate>
        OpType;
    OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow), vector_length);

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
    const int team_size       = Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
                              .team_size_recommended(op, Kokkos::ParallelForTag());
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams  = (nrow + rows_per_team - 1) / rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv_struct<MV,NoTranspose>",
                         Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);

#else   // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for
        // alpha/beta

    typedef SPMV_MV_Struct_LayoutLeft_Functor<execution_space, AMatrix, XVector, YVector, 2, 2, conjugate> OpType;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow), vector_length);

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
    const int team_size       = Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
                              .team_size_recommended(op, Kokkos::ParallelForTag());
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams  = (nrow + rows_per_team - 1) / rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv_struct<MV,NoTranspose>",
                         Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);
#endif  // KOKKOS_FAST_COMPILE
  }
}

template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta, bool conjugate>
static void spmv_alpha_beta_mv_struct_transpose(const execution_space& exec,
                                                const typename YVector::non_const_value_type& alpha, const AMatrix& A,
                                                const XVector& x, const typename YVector::non_const_value_type& beta,
                                                const YVector& y) {
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows() <= static_cast<ordinal_type>(0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  if (dobeta != 1) {
    KokkosBlas::scal(exec, y, beta, y);
  }

  if (doalpha != 0) {
    typedef typename AMatrix::size_type size_type;

    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type>(A.nnz() / A.numRows());

    int vector_length = 1;
    while ((static_cast<ordinal_type>(vector_length * 2 * 3) <= NNZPerRow) && (vector_length < 8)) vector_length *= 2;

#ifndef KOKKOS_FAST_COMPILE  // This uses templated functions on doalpha and
                             // dobeta and will produce 16 kernels

    typedef SPMV_MV_Struct_Transpose_Functor<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, conjugate>
        OpType;
    OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow));

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
    const int team_size       = Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
                              .team_size_recommended(op, Kokkos::ParallelForTag());
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams  = (nrow + rows_per_team - 1) / rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv_struct<MV,Transpose>",
                         Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);

#else  // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for
       // alpha/beta

    typedef SPMV_MV_Struct_Transpose_Functor<execution_space, AMatrix, XVector, YVector, 2, 2, conjugate, SizeType>
        OpType;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha, A, x, beta, y, RowsPerThread<execution_space>(NNZPerRow));

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<execution_space>(NNZPerRow);
    const int team_size       = Kokkos::TeamPolicy<execution_space>(exec, rows_per_thread, Kokkos::AUTO, vector_length)
                              .team_size_recommended(op, Kokkos::ParallelForTag());
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams  = (nrow + rows_per_team - 1) / rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv_struct<MV,Transpose>",
                         Kokkos::TeamPolicy<execution_space>(exec, nteams, team_size, vector_length), op);

#endif  // KOKKOS_FAST_COMPILE
  }
}

template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha, int dobeta>
static void spmv_alpha_beta_mv_struct(const execution_space& exec, const char mode[],
                                      const typename YVector::non_const_value_type& alpha, const AMatrix& A,
                                      const XVector& x, const typename YVector::non_const_value_type& beta,
                                      const YVector& y) {
  if (mode[0] == NoTranspose[0]) {
    spmv_alpha_beta_mv_struct_no_transpose<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, false>(
        exec, alpha, A, x, beta, y);
  } else if (mode[0] == Conjugate[0]) {
    spmv_alpha_beta_mv_struct_no_transpose<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, true>(
        exec, alpha, A, x, beta, y);
  } else if (mode[0] == Transpose[0]) {
    spmv_alpha_beta_mv_struct_transpose<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, false>(
        exec, alpha, A, x, beta, y);
  } else if (mode[0] == ConjugateTranspose[0]) {
    spmv_alpha_beta_mv_struct_transpose<execution_space, AMatrix, XVector, YVector, doalpha, dobeta, true>(
        exec, alpha, A, x, beta, y);
  } else {
    std::stringstream ss;
    ss << __FILE__ << ":" << __LINE__ << " Invalid transpose mode " << mode << " for KokkosSparse::spmv_struct()";
    KokkosKernels::Impl::throw_runtime_exception(ss.str());
  }
}

template <class execution_space, class AMatrix, class XVector, class YVector, int doalpha>
void spmv_alpha_mv_struct(const execution_space& exec, const char mode[],
                          const typename YVector::non_const_value_type& alpha, const AMatrix& A, const XVector& x,
                          const typename YVector::non_const_value_type& beta, const YVector& y) {
  typedef typename YVector::non_const_value_type coefficient_type;
  typedef Kokkos::ArithTraits<coefficient_type> KAT;

  if (beta == KAT::zero()) {
    spmv_alpha_beta_mv_struct<execution_space, AMatrix, XVector, YVector, doalpha, 0>(exec, mode, alpha, A, x, beta, y);
  } else if (beta == KAT::one()) {
    spmv_alpha_beta_mv_struct<execution_space, AMatrix, XVector, YVector, doalpha, 1>(exec, mode, alpha, A, x, beta, y);
  } else if (beta == -KAT::one()) {
    spmv_alpha_beta_mv_struct<execution_space, AMatrix, XVector, YVector, doalpha, -1>(exec, mode, alpha, A, x, beta,
                                                                                       y);
  } else {
    spmv_alpha_beta_mv_struct<execution_space, AMatrix, XVector, YVector, doalpha, 2>(exec, mode, alpha, A, x, beta, y);
  }
}

}  // namespace Impl
}  // namespace KokkosSparse

#endif  // KOKKOSSPARSE_IMPL_SPMV_STRUCT_DEF_HPP_
