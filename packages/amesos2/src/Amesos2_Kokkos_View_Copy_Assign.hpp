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
typename std::enable_if<std::is_same<typename dst_t::scalar_array_type,
  typename dst_t::value_type*>::value>::type
update_dst_size(dst_t & dst, const src_t & src) {
  if(dst.extent(0) != src.extent(0)) { // templated just for 1D
    dst = dst_t(Kokkos::ViewAllocateWithoutInitializing("dst"),
      src.extent(0));
  }
}

template<class dst_t, class src_t> // version for 2d view
typename std::enable_if<std::is_same<typename dst_t::scalar_array_type,
  typename dst_t::value_type**>::value>::type
update_dst_size(dst_t & dst, const src_t & src) {  // templated just for 2d
  if(dst.extent(0) != src.extent(0) || dst.extent(1) != src.extent(1)) {
    dst = dst_t(Kokkos::ViewAllocateWithoutInitializing("dst"),
      src.extent(0), src.extent(1));
  }
}

// implement_copy_or_assign has 2 versions for matched memory and
// mismatched memory. Right now we just check the memory space.
// a layout mismatch is going to compile fail so probably reflects an error
// in the initial setup.
template<class dst_t, class src_t> // version for same memory spaces
typename std::enable_if<std::is_same<typename dst_t::memory_space,
  typename src_t::memory_space>::value>::type
implement_copy_or_assign(dst_t & dst, const src_t & src) {
  std::cout << "Assign View: " << src_t::memory_space::name() <<
    " to " << dst_t::memory_space::name() << std::endl;
  dst = src; // just assign the ptr - no need to copy
}

template<class dst_t, class src_t> // version for different memory spaces
typename std::enable_if<!std::is_same<typename dst_t::memory_space,
  typename src_t::memory_space>::value>::type
implement_copy_or_assign(dst_t & dst, const src_t & src) {
  std::cout << "Deep Copy View: " << src_t::memory_space::name() <<
    " to " << dst_t::memory_space::name() << std::endl;
  update_dst_size(dst, src); // allocates if necessary
  Kokkos::deep_copy(dst, src); // full copy
}

// MDM-TODO: Remove this test system when a general adapter is implemented which
// can deliver any requested memory space. This test system was a quick way to
// have the existing Tpetra adapter deliver different memory types. I'm using
// a static to avoid sending the param through the pipe line and then remove
// it later. A general adapter will completely replace this.
class TestChooseMemorySpace {
  public:
    static std::string & src() {
      static std::string src_memory_space_name = "Default";
      return src_memory_space_name;
    }
};

// run test mode using the test_memory_t as the memory space of the src
template<class dst_t, class src_t, class test_memory_t>
void implement_test(dst_t & dst, const src_t & src) {
  // create empty view in the requested test memory space
  Kokkos::View<typename dst_t::scalar_array_type,
                       typename src_t::array_layout,
                       test_memory_t> test_src;
  // now first copy will allocate test_src and fill it
  implement_copy_or_assign(test_src, src);
  // now we run the normal copy procedure but using the test as the src
  implement_copy_or_assign(dst, test_src);
}

// this is the testing mode which will use the param string to determine
// a new src memory space and copy the data to that before continuing with
// the normal copy manager procedure.
template<class dst_t, class src_t>
void test_deep_copy_or_assign_view(dst_t & dst, const src_t & src) {
  const std::string & src_memory_space_name = TestChooseMemorySpace::src();
  if(src_memory_space_name == "Host") {
    implement_test<dst_t, src_t, Kokkos::HostSpace>(dst, src);
  }
  else if(src_memory_space_name == "CudaUVM") {
#ifdef KOKKOS_ENABLE_CUDA
    implement_test<dst_t, src_t, Kokkos::CudaUVMSpace>(dst, src);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true,
      std::runtime_error, "Cannot test CudaUVM since Cuda is not enabled.");
#endif
  }
  else if(src_memory_space_name == "Cuda") {
#ifdef KOKKOS_ENABLE_CUDA
    implement_test<dst_t, src_t, Kokkos::CudaSpace>(dst, src);
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true,
      std::runtime_error, "Cannot test Cuda since Cuda is not enabled.");
#endif
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
      "deep_copy_or_assign_ptr cannot understand argument "
      "test-src-memory-space-name set to "
        << src_memory_space_name << ".");
  }
}

// This method is currently used by the adapters.
// It will call above methods to copy or assign the view.
// If we are testing src memory spaces, this will convert the view to that
// memory space first, then feed it to the implement_copy_or_assign method.
template<class dst_t, class src_t>
void deep_copy_or_assign_view(dst_t & dst, const src_t & src) {
  // Default means we are not testing - do normal thing and copy or assign
  if(TestChooseMemorySpace::src() == "Default") {
    implement_copy_or_assign<dst_t, src_t>(dst, src);
  }
  else {
    // otherwise test - this will copy to the requested memory space, then
    // call the normal implement_copy_or_assign so it simulates different
    // source memory spaces.
    test_deep_copy_or_assign_view(dst, src);
  }
}

} // end namespace Amesos2

#endif  // AMESOS2_KOKKOS_VIEW_COPY_ASSIGN_HPP
