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

#include <fstream>

#define KOKKOSKERNELS_DEBUG_LEVEL 0

#include "Kokkos_Core.hpp"
#include "Kokkos_Timer.hpp"
#include "Kokkos_Random.hpp"
#include "Kokkos_UnorderedMap.hpp"
#include "Kokkos_Sort.hpp"

/// KokkosKernels headers
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Vector.hpp"
#include "KokkosKernels_IOUtils.hpp"

#include <Kokkos_ArithTraits.hpp>
#include <KokkosBatched_Util.hpp>
#include "examples_helper.hpp"
#include <KokkosBatched_Spmv.hpp>
#include <KokkosBatched_GMRES.hpp>
#include <KokkosBatched_CrsMatrix.hpp>
#include <KokkosBatched_Krylov_Handle.hpp>
#include <KokkosBatched_JacobiPrec.hpp>

typedef Kokkos::DefaultExecutionSpace exec_space;

template <typename DeviceType, typename ValuesViewType, typename IntView, typename VectorViewType,
          typename KrylovHandleType, bool UsePrec>
struct Functor_TestBatchedTeamVectorGMRES {
  const ValuesViewType _values;
  const ValuesViewType _diag;
  const IntView _r;
  const IntView _c;
  const VectorViewType _X;
  const VectorViewType _B;
  const int _team_size, _vector_length;
  KrylovHandleType _handle;

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVectorGMRES(const ValuesViewType &values, const IntView &r, const IntView &c,
                                     const VectorViewType &X, const VectorViewType &B, const int team_size,
                                     const int vector_length, KrylovHandleType &handle)
      : _values(values),
        _r(r),
        _c(c),
        _X(X),
        _B(B),
        _team_size(team_size),
        _vector_length(vector_length),
        _handle(handle) {}

  KOKKOS_INLINE_FUNCTION
  Functor_TestBatchedTeamVectorGMRES(const ValuesViewType &values, const ValuesViewType &diag, const IntView &r,
                                     const IntView &c, const VectorViewType &X, const VectorViewType &B,
                                     const int team_size, const int vector_length, KrylovHandleType &handle)
      : _values(values),
        _diag(diag),
        _r(r),
        _c(c),
        _X(X),
        _B(B),
        _team_size(team_size),
        _vector_length(vector_length),
        _handle(handle) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int first_matrix = _handle.first_index(member.league_rank());
    const int last_matrix  = _handle.last_index(member.league_rank());
    using TeamVectorCopy1D = KokkosBatched::TeamVectorCopy<MemberType, KokkosBatched::Trans::NoTranspose, 1>;

    auto d = Kokkos::subview(_values, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto x = Kokkos::subview(_X, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);
    auto b = Kokkos::subview(_B, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL);

    using ScratchPadIntViewType = Kokkos::View<typename IntView::non_const_value_type *, typename IntView::array_layout,
                                               typename IntView::execution_space::scratch_memory_space>;
    using ScratchPadValuesViewType =
        Kokkos::View<typename ValuesViewType::non_const_value_type **, typename ValuesViewType::array_layout,
                     typename ValuesViewType::execution_space::scratch_memory_space>;

    using Operator = KokkosBatched::CrsMatrix<ValuesViewType, ScratchPadIntViewType>;

    ScratchPadIntViewType tmp_1D_int(member.team_scratch(0), _r.extent(0) + _c.extent(0));

    auto r = Kokkos::subview(tmp_1D_int, Kokkos::make_pair(0, (int)_r.extent(0)));
    auto c = Kokkos::subview(tmp_1D_int, Kokkos::make_pair((int)_r.extent(0), (int)tmp_1D_int.extent(0)));

    TeamVectorCopy1D::invoke(member, _r, r);
    TeamVectorCopy1D::invoke(member, _c, c);
    Operator A(d, r, c);

    if (UsePrec) {
      ScratchPadValuesViewType diag(member.team_scratch(0), last_matrix - first_matrix, _diag.extent(1));
      using PrecOperator = KokkosBatched::JacobiPrec<ScratchPadValuesViewType>;

      KokkosBatched::TeamVectorCopy<MemberType>::invoke(
          member, Kokkos::subview(_diag, Kokkos::make_pair(first_matrix, last_matrix), Kokkos::ALL), diag);
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

    Kokkos::TeamPolicy<DeviceType> auto_policy(_handle.get_number_of_teams(), Kokkos::AUTO(), Kokkos::AUTO());
    Kokkos::TeamPolicy<DeviceType> tuned_policy(_handle.get_number_of_teams(), _team_size, _vector_length);
    Kokkos::TeamPolicy<DeviceType> policy;

    if (_team_size < 1)
      policy = auto_policy;
    else
      policy = tuned_policy;

    int maximum_iteration = _handle.get_max_iteration();

    using ScalarType = typename ValuesViewType::non_const_value_type;
    using Layout     = typename ValuesViewType::array_layout;
    using EXSP       = typename ValuesViewType::execution_space;

    using ViewType2D = Kokkos::View<ScalarType **, Layout, EXSP>;

    size_t bytes_1D      = ViewType2D::shmem_size(_handle.get_number_of_systems_per_team(), 1);
    size_t bytes_row_ptr = IntView::shmem_size(_r.extent(0));
    size_t bytes_col_idc = IntView::shmem_size(_c.extent(0));
    size_t bytes_2D_1    = ViewType2D::shmem_size(_handle.get_number_of_systems_per_team(), _X.extent(1));
    size_t bytes_2D_2    = ViewType2D::shmem_size(_handle.get_number_of_systems_per_team(), maximum_iteration + 1);

    size_t bytes_int  = bytes_row_ptr + bytes_col_idc;
    size_t bytes_diag = bytes_2D_1;
    size_t bytes_tmp  = 2 * bytes_2D_1 + 2 * bytes_1D + bytes_2D_2;

    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_tmp + bytes_diag + bytes_int));

    exec_space().fence();
    timer.reset();
    Kokkos::parallel_for(name.c_str(), policy, *this);
    exec_space().fence();
    double sec = timer.seconds();

    return sec;
  }
};

int main(int /*argc*/, char ** /*argv*/) {
  Kokkos::initialize();
  {
    using layout = Kokkos::LayoutLeft;

    using IntView          = Kokkos::View<int *, layout, exec_space>;
    using AMatrixValueView = Kokkos::View<double **, layout, exec_space>;
    using XYType           = Kokkos::View<double **, layout, exec_space>;

    std::string name_A = "mat.mm";
    std::string name_B = "rhs.mm";

    int N, Blk, nnz;

    Blk = 10;
    N   = 100;
    nnz = (Blk - 2) * 3 + 2 * 2;

    IntView rowOffsets("rowOffsets", Blk + 1);
    IntView colIndices("colIndices", nnz);
    AMatrixValueView values("values", N, nnz);
    AMatrixValueView diag("diag", N, Blk);
    XYType x("x", N, Blk);
    XYType y("y", N, Blk);

    printf("N = %d, Blk = %d, nnz = %d\n", N, Blk, nnz);

    create_tridiagonal_batched_matrices(nnz, Blk, N, rowOffsets, colIndices, values, x, y);

    // Replace y by ones:
    Kokkos::deep_copy(y, 1.);

    // Replace x by zeros:
    // Kokkos::deep_copy(x, 0.);

    getInvDiagFromCRS(values, rowOffsets, colIndices, diag);

    using ScalarType = typename AMatrixValueView::non_const_value_type;
    using Layout     = typename AMatrixValueView::array_layout;
    using EXSP       = typename AMatrixValueView::execution_space;

    using MagnitudeType = typename Kokkos::ArithTraits<ScalarType>::mag_type;

    using Norm2DViewType   = Kokkos::View<MagnitudeType **, Layout, EXSP>;
    using Scalar3DViewType = Kokkos::View<ScalarType ***, Layout, EXSP>;
    using IntViewType      = Kokkos::View<int *, Layout, EXSP>;

    using KrylovHandleType = KokkosBatched::KrylovHandle<Norm2DViewType, IntViewType, Scalar3DViewType>;

    const int N_team       = 2;
    const int n_iterations = 150;

    const int team_size      = -1;
    const int vector_length  = -1;
    const double tol         = 1e-8;
    const int ortho_strategy = 0;

    KrylovHandleType handle(N, N_team, n_iterations, true);
    handle.Arnoldi_view = Scalar3DViewType("", N, n_iterations, Blk + n_iterations + 3);

    handle.set_max_iteration(n_iterations);
    handle.set_tolerance(tol);
    handle.set_ortho_strategy(ortho_strategy);
    handle.set_scratch_pad_level(0);
    handle.set_compute_last_residual(true);

    double time =
        Functor_TestBatchedTeamVectorGMRES<exec_space, AMatrixValueView, IntView, XYType, KrylovHandleType, true>(
            values, diag, rowOffsets, colIndices, x, y, team_size, vector_length, handle)
            .run();

    printf("times = %f secondes\n", time);

    for (int i = 0; i < N; ++i) {
      if (handle.is_converged_host(i)) {
        std::cout << "System " << i << " converged in " << handle.get_iteration_host(i)
                  << " iterations, the initial absolute norm of the residual was " << handle.get_norm_host(i, 0)
                  << " and is now " << handle.get_last_norm_host(i) << std::endl;
      } else {
        std::cout << "System " << i << " did not converge in " << handle.get_max_iteration()
                  << " iterations, the initial absolute norm of the residual was " << handle.get_norm_host(i, 0)
                  << " and is now " << handle.get_last_norm_host(i) << std::endl;
      }
    }
    if (handle.is_converged_host())
      std::cout << "All the systems have converged." << std::endl;
    else
      std::cout << "There is at least one system that did not converge." << std::endl;
  }
  Kokkos::finalize();
}
