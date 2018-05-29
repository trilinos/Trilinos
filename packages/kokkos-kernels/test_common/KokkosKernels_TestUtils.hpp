/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
}

#endif
