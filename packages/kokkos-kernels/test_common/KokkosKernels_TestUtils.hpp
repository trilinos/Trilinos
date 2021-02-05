/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
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
*/

#ifndef KOKKOSKERNELS_TEST_UTILS_HPP
#define KOKKOSKERNELS_TEST_UTILS_HPP

#include "KokkosKernels_Utils.hpp"
namespace Test {
  template<class ViewType, bool strided = std::is_same<typename ViewType::array_layout, Kokkos::LayoutStride>::value>
  struct multivector_layout_adapter;

  template<class ViewType>
  struct multivector_layout_adapter<ViewType,true> {
    typedef typename ViewType::value_type Scalar;
    typedef typename ViewType::device_type Device;
    typedef Kokkos::View<Scalar**[2], Kokkos::LayoutRight, Device> BaseTypeRight;
    typedef Kokkos::View<Scalar**, typename ViewType::array_layout, Device> BaseTypeDefault;
    typedef typename std::conditional<
                std::is_same<typename ViewType::array_layout,Kokkos::LayoutStride>::value,
                BaseTypeRight, BaseTypeDefault>::type BaseType;

    static ViewType view(const BaseType& v) {
      return Kokkos::subview(v,Kokkos::ALL,Kokkos::ALL,0);
    };
  };

  template<class ViewType>
  struct multivector_layout_adapter<ViewType,false> {
    typedef typename ViewType::value_type Scalar;
    typedef typename ViewType::device_type Device;
    typedef Kokkos::View<Scalar**[2], Kokkos::LayoutRight, Device> BaseTypeRight;
    typedef Kokkos::View<Scalar**, typename ViewType::array_layout, Device> BaseTypeDefault;
    typedef typename std::conditional<
                std::is_same<typename ViewType::array_layout,Kokkos::LayoutStride>::value,
                BaseTypeRight, BaseTypeDefault>::type BaseType;

    static ViewType view(const BaseType& v) {
      return Kokkos::subview(v,Kokkos::ALL,Kokkos::ALL);
    };
  };

  template<class Scalar1, class Scalar2, class Scalar3>
  void EXPECT_NEAR_KK(Scalar1 val1, Scalar2 val2, Scalar3 tol) {
    typedef Kokkos::Details::ArithTraits<Scalar1> AT1;
    typedef Kokkos::Details::ArithTraits<Scalar2> AT2;
    typedef Kokkos::Details::ArithTraits<Scalar3> AT3;
    EXPECT_NEAR(double(AT1::abs(val1)),double(AT2::abs(val2)),double(AT3::abs(tol)));
  }

  template<class ViewType1, class ViewType2, class Scalar>
  void EXPECT_NEAR_KK_1DVIEW(ViewType1 v1, ViewType2 v2, Scalar tol) {
    size_t v1_size = v1.extent(0);
    size_t v2_size = v2.extent(0);
    EXPECT_NEAR_KK(v1_size, v2_size, 0);


    typename ViewType1::HostMirror h_v1 = Kokkos::create_mirror_view(v1);
    typename ViewType2::HostMirror h_v2 = Kokkos::create_mirror_view(v2);

    KokkosKernels::Impl::safe_device_to_host_deep_copy (v1.extent(0), v1, h_v1);
    KokkosKernels::Impl::safe_device_to_host_deep_copy (v2.extent(0), v2, h_v2);

    for (size_t i = 0; i < v1_size; ++i){
      EXPECT_NEAR_KK(h_v1(i), h_v2(i), tol);
    }
  }

  #if defined(KOKKOS_HALF_T_IS_FLOAT)
  using halfScalarType = Kokkos::Experimental::half_t;
  #endif // KOKKOS_HALF_T_IS_FLOAT

  template<class ViewTypeA, class ViewTypeB, class ViewTypeC, class ExecutionSpace>
  struct SharedVanillaGEMM {
    bool A_t, B_t, A_c, B_c;
    int C_rows, C_cols, A_cols;
    ViewTypeA A;
    ViewTypeB B;
    ViewTypeC C;

    typedef typename ViewTypeA::value_type ScalarA;
    typedef typename ViewTypeB::value_type ScalarB;
    typedef typename ViewTypeC::value_type ScalarC;
    typedef Kokkos::Details::ArithTraits<ScalarC> APT;
    typedef typename APT::mag_type mag_type;
    ScalarA alpha;
    ScalarC beta;

    KOKKOS_INLINE_FUNCTION
    void operator() (const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team) const {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,C_rows), [&] (const int& i) {
        // Give each kokkos thread a vector of A
        auto a_vec = A_t ? Kokkos::subview(A, Kokkos::ALL(), i) : Kokkos::subview(A, i, Kokkos::ALL());

        // Have all vector lanes perform the dot product
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,C_cols), [&] (const int& j) {
          auto b_vec = B_t ? Kokkos::subview(B, j, Kokkos::ALL()) : Kokkos::subview(B, Kokkos::ALL(), j);
          ScalarC ab = ScalarC(0);
          for (int k = 0; k < A_cols; k++) {
            auto a = A_c ? APT::conj(a_vec(k)) : a_vec(k);
            auto b = B_c ? APT::conj(b_vec(k)) : b_vec(k);
            ab += a * b;
          }
          C(i,j) = beta * C(i,j) + alpha * ab;
        });
      });
    }
  };
  // C(i,:,:) = alpha * (A(i,:,:) * B(i,:,:)) + beta * C(i,:,:)
  template<class ViewTypeA, class ViewTypeB, class ViewTypeC, class ExecutionSpace>
  struct Functor_BatchedVanillaGEMM {
    bool A_t, B_t, A_c, B_c;
    ViewTypeA A;
    ViewTypeB B;
    ViewTypeC C;

    using ScalarA = typename ViewTypeA::value_type;
    using ScalarB = typename ViewTypeB::value_type;
    using ScalarC = typename ViewTypeC::value_type;
    ScalarA alpha;
    ScalarC beta;

    KOKKOS_INLINE_FUNCTION
    void operator()(const typename Kokkos::TeamPolicy<ExecutionSpace>::member_type& team) const {
      int i = team.league_rank();

      auto _A = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
      auto _B = Kokkos::subview(B, i, Kokkos::ALL(), Kokkos::ALL());
      auto _C = Kokkos::subview(C, i, Kokkos::ALL(), Kokkos::ALL());
      using SubviewTypeA = decltype(_A);
      using SubviewTypeB = decltype(_B);
      using SubviewTypeC = decltype(_C);
      struct SharedVanillaGEMM<SubviewTypeA,SubviewTypeB,SubviewTypeC,ExecutionSpace> vgemm;
      vgemm.A_t = A_t; vgemm.B_t = B_t;
      vgemm.A_c = A_c; vgemm.B_c = B_c;
      vgemm.C_rows = C.extent(1);
      vgemm.C_cols = C.extent(2);    
      vgemm.A_cols = A_t?A.extent(1):A.extent(2);
      vgemm.A = _A;
      vgemm.B = _B;
      vgemm.C = _C;
      vgemm.alpha = alpha;
      vgemm.beta = beta;
      vgemm(team);
    }

    inline
    void run() {
      Kokkos::parallel_for(
          "Test::VanillaGEMM",
          Kokkos::TeamPolicy<ExecutionSpace>(C.extent(0), Kokkos::AUTO, 16),
          *this);
    }
  };
}
#endif
