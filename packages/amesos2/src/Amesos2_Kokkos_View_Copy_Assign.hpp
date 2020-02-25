// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/**
  \file   Amesos2_Kokkos_View_Copy_Assign.hpp
  \author
  \date   Fri Sept 13 6:00:00 2019

  \brief  Copy or assign views based on memory spaces
*/

#ifndef AMESOS2_KOKKOS_VIEW_COPY_ASSIGN_HPP
#define AMESOS2_KOKKOS_VIEW_COPY_ASSIGN_HPP

namespace Amesos2 {

// allocate dst size if necessary - 2 methods handle 1d and 2d
template<class dst_t, class src_t> // version for 1d view
typename std::enable_if<static_cast<int>(dst_t::Rank) == 1>::type
update_dst_size(dst_t & dst, const src_t & src) {
  if(dst.extent(0) != src.extent(0)) { // templated just for 1D
    dst = dst_t(Kokkos::ViewAllocateWithoutInitializing("dst"),
      src.extent(0));
  }
}

template<class dst_t, class src_t> // version for 2d view
typename std::enable_if<static_cast<int>(dst_t::Rank) == 2>::type
update_dst_size(dst_t & dst, const src_t & src) {  // templated just for 2d
  if(dst.extent(0) != src.extent(0) || dst.extent(1) != src.extent(1)) {
    dst = dst_t(Kokkos::ViewAllocateWithoutInitializing("dst"),
      src.extent(0), src.extent(1));
  }
}

// now handle type mismatch for same memory space - here types are same
template<class dst_t, class src_t> // version for same memory spaces
typename std::enable_if<std::is_same<typename dst_t::value_type,
  typename src_t::value_type>::value>::type
implement_copy_or_assign_same_mem_check_types(dst_t & dst, const src_t & src) {
  dst = src; // just assign the ptr - no need to copy
}

// now handle type mismatch for same memory space - now types are different
template<class dst_t, class src_t> // version for same memory spaces
typename std::enable_if<!std::is_same<typename dst_t::value_type,
  typename src_t::value_type>::value>::type
implement_copy_or_assign_same_mem_check_types(dst_t & dst, const src_t & src) {
  update_dst_size(dst, src); // allocates if necessary
  Kokkos::deep_copy(dst, src); // full copy
}

// implement_copy_or_assign has 2 versions for matched memory and
// mismatched memory. Right now we just check the memory space.
// a layout mismatch is going to compile fail so probably reflects an error
// in the initial setup.
template<class dst_t, class src_t> // version for same memory spaces
typename std::enable_if<std::is_same<typename dst_t::memory_space,
  typename src_t::memory_space>::value>::type
deep_copy_or_assign_view(dst_t & dst, const src_t & src) {
  implement_copy_or_assign_same_mem_check_types(dst, src);
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<std::is_same<typename dst_t::value_type,
  typename src_t::value_type>::value>::type
implement_copy_or_assign_diff_mem_check_types(dst_t & dst, const src_t & src) {
  update_dst_size(dst, src); // allocates if necessary
  Kokkos::deep_copy(dst, src); // full copy
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<static_cast<int>(dst_t::Rank) == 1>::type
implement_copy_or_assign_diff_mem_diff_types_check_dim(dst_t & dst, const src_t & src) {
  Kokkos::View<typename dst_t::value_type*, typename src_t::execution_space>
    intermediate(Kokkos::ViewAllocateWithoutInitializing("intermediate"), src.extent(0));
  Kokkos::deep_copy(intermediate, src); // to dst type
  Kokkos::deep_copy(dst, intermediate); // to dst mem
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<static_cast<int>(dst_t::Rank) == 2>::type
implement_copy_or_assign_diff_mem_diff_types_check_dim(dst_t & dst, const src_t & src) {
  Kokkos::View<typename dst_t::value_type**, Kokkos::LayoutLeft, typename src_t::execution_space>
    intermediate(Kokkos::ViewAllocateWithoutInitializing("intermediate"), src.extent(0), src.extent(1));
  Kokkos::deep_copy(intermediate, src); // to dst type
  Kokkos::deep_copy(dst, intermediate); // to dst mem
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<!std::is_same<typename dst_t::value_type,
  typename src_t::value_type>::value>::type
implement_copy_or_assign_diff_mem_check_types(dst_t & dst, const src_t & src) {
  update_dst_size(dst, src); // allocates if necessary
  // since mem space and types are different, we specify the order of operations
  // Kokkos::deep_copy won't do both since it would be a hidden deep_copy
  implement_copy_or_assign_diff_mem_diff_types_check_dim(dst, src);
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<!std::is_same<typename dst_t::memory_space,
  typename src_t::memory_space>::value>::type
deep_copy_or_assign_view(dst_t & dst, const src_t & src) {
  implement_copy_or_assign_diff_mem_check_types(dst, src); // full copy
}

// special utility method to hard cast src Kokkos::complex to std::complex
// to handle case where Kokkos::View needs to send data to Tpetra method
// in std::complex form. If the src is not Kokkos::complex this method is just
// going to call the regular deep_copy_or_assign_view above.
// MDM-TODO This should be improved but set it up to open discussion on how to
// handle some std::complex Kokkos::complex conversion issues.
// I would like to remove the reinterpret_cast but not sure how to do that and
// preserve the optimization pathways this file sets up.
template<class dst_t, class src_t> // Kokkos::complex<real_t> version
typename std::enable_if<
  std::is_same<typename src_t::value_type, typename Kokkos::complex<double>>::value
  >::type
deep_copy_or_assign_view_make_src_std_complex(dst_t & dst, src_t & src) {
  typedef const Kokkos::View<std::complex<double>**, Kokkos::LayoutLeft, typename src_t::execution_space> std_complex_t;
  // MDM-TODO - make this safer
  std_complex_t * std_complex_src = reinterpret_cast<std_complex_t*>(&src);
  deep_copy_or_assign_view(dst, *std_complex_src);
}

template<class dst_t, class src_t> // Kokkos::complex<real_t> version
typename std::enable_if<
  std::is_same<typename src_t::value_type, typename Kokkos::complex<float>>::value
  >::type
deep_copy_or_assign_view_make_src_std_complex(dst_t & dst, src_t & src) {
  typedef const Kokkos::View<std::complex<float>**, Kokkos::LayoutLeft, typename src_t::execution_space> std_complex_t;
  // MDM-TODO - make this safer
  std_complex_t * std_complex_src = reinterpret_cast<std_complex_t*>(&src);
  deep_copy_or_assign_view(dst, *std_complex_src);
}

template<class dst_t, class src_t> // Kokkos::complex<real_t> version
typename std::enable_if<
  !std::is_same<typename src_t::value_type, typename Kokkos::complex<double>>::value &&
  !std::is_same<typename src_t::value_type, typename Kokkos::complex<float>>::value
  >::type
deep_copy_or_assign_view_make_src_std_complex(dst_t & dst, src_t & src) {
  deep_copy_or_assign_view(dst, src);
}

} // end namespace Amesos2

#endif  // AMESOS2_KOKKOS_VIEW_COPY_ASSIGN_HPP
