//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.4
//       Copyright (2021) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
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
#include <KokkosBatched_Trsv_Decl.hpp>
#include <KokkosBatched_Trsv_Serial_Impl.hpp>
#include <KokkosBatched_Trsv_Team_Impl.hpp>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_LU_Serial_Impl.hpp>
#include <KokkosBatched_LU_Team_Impl.hpp>
#include "KokkosBatched_Gesv.hpp"

typedef Kokkos::DefaultExecutionSpace exec_space;

template <typename DeviceType, typename AViewType, typename XYViewType>
struct Functor_TeamTestStaticPivoting {
  const AViewType _A;
  const XYViewType _X;
  const XYViewType _Y;

  KOKKOS_INLINE_FUNCTION
  Functor_TeamTestStaticPivoting(const AViewType &A, const XYViewType &X,
                                 const XYViewType &Y)
      : _A(A), _X(X), _Y(Y) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const {
    const int matrix_id = static_cast<int>(member.league_rank());

    auto A = Kokkos::subview(_A, matrix_id, Kokkos::ALL, Kokkos::ALL);
    auto X = Kokkos::subview(_X, matrix_id, Kokkos::ALL);
    auto Y = Kokkos::subview(_Y, matrix_id, Kokkos::ALL);
    member.team_barrier();
    KokkosBatched::TeamGesv<MemberType,
                            KokkosBatched::Gesv::StaticPivoting>::invoke(member,
                                                                         A, X,
                                                                         Y);
    member.team_barrier();
  }

  inline void run() {
    std::string name("KokkosBatched::Test::StaticPivoting");
    Kokkos::TeamPolicy<DeviceType> policy(_A.extent(0), Kokkos::AUTO(),
                                          Kokkos::AUTO());

    using MatrixViewType =
        Kokkos::View<typename AViewType::non_const_value_type **,
                     typename AViewType::array_layout,
                     typename AViewType::execution_space>;

    const int n    = _A.extent(1);
    size_t bytes_0 = MatrixViewType::shmem_size(n, n + 4);

    policy.set_scratch_size(0, Kokkos::PerTeam(bytes_0));

    Kokkos::parallel_for(name.c_str(), policy, *this);
  }
};

template <typename DeviceType, typename AViewType, typename XYViewType>
struct Functor_SerialTestStaticPivoting {
  const AViewType _A;
  const AViewType _tmp;
  const XYViewType _X;
  const XYViewType _Y;

  KOKKOS_INLINE_FUNCTION
  Functor_SerialTestStaticPivoting(const AViewType &A, const AViewType &tmp,
                                   const XYViewType &X, const XYViewType &Y)
      : _A(A), _tmp(tmp), _X(X), _Y(Y) {}

  KOKKOS_INLINE_FUNCTION void operator()(const int &matrix_id) const {
    auto A   = Kokkos::subview(_A, matrix_id, Kokkos::ALL, Kokkos::ALL);
    auto tmp = Kokkos::subview(_tmp, matrix_id, Kokkos::ALL, Kokkos::ALL);
    auto X   = Kokkos::subview(_X, matrix_id, Kokkos::ALL);
    auto Y   = Kokkos::subview(_Y, matrix_id, Kokkos::ALL);
    KokkosBatched::SerialGesv<KokkosBatched::Gesv::StaticPivoting>::invoke(
        A, X, Y, tmp);
  }

  inline void run() {
    std::string name("KokkosBatched::Test::StaticPivoting");

    const int N = _A.extent(0);
    Kokkos::parallel_for(name.c_str(), N, *this);
  }
};

int main(int /*argc*/, char ** /*argv[]*/) {
  Kokkos::initialize();
  {
    using layout = Kokkos::LayoutLeft;

    using AViewType  = Kokkos::View<double ***, layout, exec_space>;
    using XYViewType = Kokkos::View<double **, layout, exec_space>;

    int N = 1;
    int n = 10;

    AViewType A("A", N, n, n);
    AViewType tmp("tmp", N, n, n + 4);
    XYViewType X("X", N, n);
    XYViewType Y("Y", N, n);

    create_saddle_point_matrices(A, Y);

    // The matrices are modified by the GESV so we have to copy them if we want
    // to solve the same systems twice.
    AViewType A2("A2", N, n, n);
    XYViewType Y2("Y2", N, n);
    Kokkos::deep_copy(A2, A);
    Kokkos::deep_copy(Y2, Y);

    KokkosKernels::Impl::kk_write_3Dview_to_file(A, "A.txt");
    KokkosKernels::Impl::kk_write_2Dview_to_file(Y, "Y.txt");

    Functor_SerialTestStaticPivoting<exec_space, AViewType, XYViewType>(A, tmp,
                                                                        X, Y)
        .run();
    KokkosKernels::Impl::kk_write_2Dview_to_file(X, "X_serial.txt");
    Functor_TeamTestStaticPivoting<exec_space, AViewType, XYViewType>(A2, X, Y2)
        .run();
    KokkosKernels::Impl::kk_write_2Dview_to_file(X, "X_team.txt");
  }
  Kokkos::finalize();
}
