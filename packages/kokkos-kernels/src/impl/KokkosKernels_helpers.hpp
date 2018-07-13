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
#ifndef KOKKOSKERNELS_HELPERS_HPP_
#define KOKKOSKERNELS_HELPERS_HPP_

namespace KokkosKernels {
namespace Impl {

// Unify Layout of a View to LayoutLeft if possible.
// Used to reduce number of code instantiations

template<class ViewType>
struct GetUnifiedLayout {
  typedef typename std::conditional<
        ( (ViewType::rank == 1) &&
          (std::is_same<typename ViewType::array_layout,Kokkos::LayoutRight>::value) ) ||
        ( (ViewType::rank == 0) )
       ,Kokkos::LayoutLeft,typename ViewType::array_layout>::type array_layout;
};

template<class T, class TX, bool do_const, bool isView = Kokkos::is_view<T>::value>
struct GetUnifiedScalarViewType {
  typedef typename TX::non_const_value_type type;
};

template<class T, class TX>
struct GetUnifiedScalarViewType<T,TX,false,true> {
  typedef Kokkos::View<typename T::non_const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<T>::array_layout,
                       typename T::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > type;
};

template<class T, class TX>
struct GetUnifiedScalarViewType<T,TX,true,true> {
  typedef Kokkos::View<typename T::const_value_type*,
                       typename KokkosKernels::Impl::GetUnifiedLayout<T>::array_layout,
                       typename T::device_type,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged> > type;
};

}
}
#endif
