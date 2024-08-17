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
          typename KrylovHandleType>
struct Functor_TestBatchedTeamVectorCG_1 {
  const ValuesViewType _D;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;
  const int _N_team, _team_size, _vector_length;
  KrylovHandleType _handle;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVectorCG_1(const ValuesViewType &D, const IntView &r, const IntView &c,
                                    const VectorViewType &X, const VectorViewType &B, const int N_team,
                                    const int team_size, const int vector_length, KrylovHandleType &handle)
      : _D(D),
        _r(r),
        _c(c),
        _X(X),
        _B(B),
        _N_team(N_team),
        _team_size(team_size),
        _vector_length(vector_length),
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

    KokkosBatched::TeamVectorCG<MemberType>::template invoke<Operator, VectorViewType>(member, A, b, x, _handle);
  }

  inline double run() {
    std::string name("KokkosBatched::Test::TeamCG");
    Kokkos::Timer timer;
    Kokkos::Profiling::pushRegion(name.c_str());

    _handle.set_memory_strategy(1);

    _handle.tmp_view = typename KrylovHandleType::TemporaryViewType("", _X.extent(0), 4 * _X.extent(1));

    Kokkos::TeamPolicy<DeviceType> auto_policy(_handle.get_number_of_teams(), Kokkos::AUTO(), Kokkos::AUTO());
    Kokkos::TeamPolicy<DeviceType> tuned_policy(_handle.get_number_of_teams(), _team_size, _vector_length);
    Kokkos::TeamPolicy<DeviceType> policy;

    if (_team_size < 1)
      policy = auto_policy;
    else
      policy = tuned_policy;

    size_t bytes_0 = ValuesViewType::shmem_size(_N_team, 5);
    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_0));

    timer.reset();
    Kokkos::parallel_for(name.c_str(), policy, *this);
    exec_space().fence();
    double sec = timer.seconds();
    Kokkos::Profiling::popRegion();

    return sec;
  }
};