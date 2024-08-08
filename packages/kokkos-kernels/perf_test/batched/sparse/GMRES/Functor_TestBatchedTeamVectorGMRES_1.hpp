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

template <typename DeviceType, typename ValuesViewType, typename IntView, typename VectorViewType,
          typename KrylovHandleType, bool UsePrec>
struct Functor_TestBatchedTeamVectorGMRES_1 {
  const ValuesViewType _D;
  const ValuesViewType _diag;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;
  const int _N_team, _team_size, _vector_length;
  const int _N_iteration;
  const double _tol;
  const int _ortho_strategy;
  const int _arnoldi_level;
  const int _other_level;
  KrylovHandleType _handle;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVectorGMRES_1(const ValuesViewType &D, const IntView &r, const IntView &c,
                                       const VectorViewType &X, const VectorViewType &B, const int N_team,
                                       const int team_size, const int vector_length, const int N_iteration,
                                       const double tol, const int ortho_strategy, const int arnoldi_level,
                                       const int other_level, KrylovHandleType &handle)
      : _D(D),
        _r(r),
        _c(c),
        _X(X),
        _B(B),
        _N_team(N_team),
        _team_size(team_size),
        _vector_length(vector_length),
        _N_iteration(N_iteration),
        _tol(tol),
        _ortho_strategy(ortho_strategy),
        _arnoldi_level(arnoldi_level),
        _other_level(other_level),
        _handle(handle) {}

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVectorGMRES_1(const ValuesViewType &D, const ValuesViewType &diag, const IntView &r,
                                       const IntView &c, const VectorViewType &X, const VectorViewType &B,
                                       const int N_team, const int team_size, const int vector_length,
                                       const int N_iteration, const double tol, int ortho_strategy,
                                       const int arnoldi_level, const int other_level, KrylovHandleType &handle)
      : _D(D),
        _diag(diag),
        _r(r),
        _c(c),
        _X(X),
        _B(B),
        _N_team(N_team),
        _team_size(team_size),
        _vector_length(vector_length),
        _N_iteration(N_iteration),
        _tol(tol),
        _ortho_strategy(ortho_strategy),
        _arnoldi_level(arnoldi_level),
        _other_level(other_level),
        _handle(handle) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int first_matrix = _handle.first_index(member.league_rank());
    const int last_matrix  = _handle.last_index(member.league_rank());

    auto d = Kokkos::subview(_D, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto x = Kokkos::subview(_X, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto b = Kokkos::subview(_B, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);

    using Operator = KokkosBatched::CrsMatrix<ValuesViewType, IntView>;

    Operator A(d, _r, _c);

    if (UsePrec) {
      auto diag          = Kokkos::subview(_diag, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
      using PrecOperator = KokkosBatched::JacobiPrec<VectorViewType>;

      PrecOperator P(diag);
      P.setComputedInverse();

      KokkosBatched::TeamVectorGMRES<MemberType>::template invoke<Operator, VectorViewType, PrecOperator,
                                                                  KrylovHandleType>(member, A, b, x, P, _handle);
    } else {
      KokkosBatched::TeamVectorGMRES<MemberType>::template invoke<Operator, VectorViewType>(member, A, b, x, _handle);
    }
  }

  inline double run() {
    std::string name("KokkosBatched::Test::TeamVectorGMRES");
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion(name.c_str());

    int maximum_iteration = _handle.get_max_iteration();

    _handle.set_memory_strategy(1);

    _handle.tmp_view =
        typename KrylovHandleType::TemporaryViewType("", _X.extent(0), _X.extent(1) + maximum_iteration + 3);

    Kokkos::TeamPolicy<DeviceType> auto_policy(_handle.get_number_of_teams(), Kokkos::AUTO(), Kokkos::AUTO());
    Kokkos::TeamPolicy<DeviceType> tuned_policy(_handle.get_number_of_teams(), _team_size, _vector_length);
    Kokkos::TeamPolicy<DeviceType> policy;

    if (_team_size < 1)
      policy = auto_policy;
    else
      policy = tuned_policy;

    exec_space().fence();
    timer.reset();
    Kokkos::parallel_for(name.c_str(), policy, *this);
    exec_space().fence();
    double sec = timer.seconds();
    Kokkos::Profiling::popRegion();

    return sec;
  }
};