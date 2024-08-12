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

#ifndef KOKKOS_SPMV_INSPECTOR_HPP_
#define KOKKOS_SPMV_INSPECTOR_HPP_

#include "Kokkos_SPMV.hpp"

template <class AMatrix, class XVector, class YVector, int dobeta, bool conjugate>
struct SPMV_Inspector_Functor {
  typedef typename AMatrix::execution_space execution_space;
  typedef typename AMatrix::non_const_ordinal_type ordinal_type;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename AMatrix::non_const_size_type size_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type team_member;
  typedef Kokkos::ArithTraits<value_type> ATV;

  const value_type alpha;
  AMatrix m_A;
  XVector m_x;
  Kokkos::View<const ordinal_type*> m_workset_offsets;

  const value_type beta;
  YVector m_y;

  SPMV_Inspector_Functor(const value_type alpha_, const AMatrix m_A_, const XVector m_x_,
                         const Kokkos::View<const ordinal_type*> m_workset_offsets_, const value_type beta_,
                         const YVector m_y_)
      : alpha(alpha_), m_A(m_A_), m_x(m_x_), m_workset_offsets(m_workset_offsets_), beta(beta_), m_y(m_y_) {
    static_assert(static_cast<int>(XVector::rank) == 1, "XVector must be a rank 1 View.");
    static_assert(static_cast<int>(YVector::rank) == 1, "YVector must be a rank 1 View.");
  }

  KOKKOS_INLINE_FUNCTION void operator()(const team_member& dev) const {
    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(dev, m_workset_offsets(dev.league_rank()), m_workset_offsets(dev.league_rank() + 1)),
        [&](const ordinal_type& iRow) {
          // const ordinal_type iRow = static_cast<ordinal_type> (
          // dev.league_rank() ) * rows_per_team + loop;
          if (iRow >= m_A.numRows()) {
            return;
          }
          const KokkosSparse::SparseRowViewConst<AMatrix> row = m_A.rowConst(iRow);
          const ordinal_type row_length                       = static_cast<ordinal_type>(row.length);
          value_type sum                                      = 0;

          Kokkos::parallel_reduce(
              Kokkos::ThreadVectorRange(dev, row_length),
              [&](const ordinal_type& iEntry, value_type& lsum) {
                const value_type val = conjugate ? ATV::conj(row.value(iEntry)) : row.value(iEntry);
                lsum += val * m_x(row.colidx(iEntry));
              },
              sum);

          Kokkos::single(Kokkos::PerThread(dev), [&]() {
            sum *= alpha;

            if (dobeta == 0) {
              m_y(iRow) = sum;
            } else {
              m_y(iRow) = beta * m_y(iRow) + sum;
            }
          });
        });
  }
};

template <typename AType, typename XType, typename YType, class Schedule>
void kk_inspector_matvec(AType A, XType x, YType y, int team_size, int vector_length) {
  typedef typename AType::execution_space execution_space;
  typedef typename AType::device_type::memory_space memory_space;
  typedef typename AType::non_const_size_type size_type;
  typedef typename AType::non_const_ordinal_type lno_t;
  typedef typename AType::non_const_value_type scalar_t;

  int const concurrency = execution_space().concurrency();
  static int worksets =
      std::is_same<Schedule, Kokkos::Static>::value ? team_size > 0 ? concurrency / team_size : concurrency :  // static
          team_size > 0 ? concurrency * 32 / team_size
                        : concurrency * 32;  // dynamic
  static Kokkos::View<lno_t*, memory_space> workset_offsets;
  if (workset_offsets.extent(0) == 0) {
    workset_offsets       = Kokkos::View<lno_t*>("WorksetOffsets", worksets + 1);
    const size_type nnz   = A.nnz();
    lno_t nnz_per_workset = (nnz + worksets - 1) / worksets;
    workset_offsets(0)    = 0;
    lno_t ws              = 1;
    for (lno_t row = 0; row < A.numRows(); row++) {
      if (A.graph.row_map(row) > size_type(ws) * nnz_per_workset) {
        workset_offsets(ws) = row;
        ws++;
      }
    }
    if (workset_offsets(ws - 1) < A.numRows()) {
      workset_offsets(ws) = A.numRows();
    }
    printf("Worksets: %i %i\n", worksets, ws);
    worksets = ws;
  }
  scalar_t s_a(1.0);
  scalar_t s_b(0.0);
  SPMV_Inspector_Functor<AType, XType, YType, 0, false> func(s_a, A, x, workset_offsets, s_b, y);

  Kokkos::TeamPolicy<Kokkos::Schedule<Schedule> > policy(1, 1);

  if (team_size > 0)
    policy = Kokkos::TeamPolicy<Kokkos::Schedule<Schedule> >(worksets, team_size, vector_length);
  else
    policy = Kokkos::TeamPolicy<Kokkos::Schedule<Schedule> >(worksets, Kokkos::AUTO, vector_length);

  Kokkos::parallel_for("KokkosSparse::PerfTest::SpMV_Inspector", policy, func);
}

#endif /* KOKKOS_SPMV_HPP_ */
