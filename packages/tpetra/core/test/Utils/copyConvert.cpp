/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_copyConvert.hpp"

namespace { // (anonymous)

  //
  // UNIT TESTS
  //

  template<class OutputViewType, class InputViewType>
  void testCopyConvert (bool& success, std::ostream& out)
  {
    typedef typename OutputViewType::non_const_value_type output_value_type;
    typedef typename InputViewType::non_const_value_type input_value_type;
    typedef typename OutputViewType::execution_space output_execution_space;
    typedef typename InputViewType::execution_space input_execution_space;

    const int size = 42;

    typename InputViewType::non_const_type x ("x", size);
    auto x_h = Kokkos::create_mirror_view (x);
    input_execution_space::fence (); // for UVM's sake
    for (int i = 0; i < size; ++i) {
      x_h(i) = static_cast<input_value_type> (i+1); // no entries are zero
    }
    Kokkos::deep_copy (x, x_h);
    input_execution_space::fence (); // for UVM's sake

    typename OutputViewType::non_const_type y ("y", size);
    output_execution_space::fence (); // for UVM's sake
    Tpetra::Details::copyConvert (y, x);
    output_execution_space::fence (); // for UVM's sake

    auto y_h = Kokkos::create_mirror_view (y);
    output_execution_space::fence (); // for UVM's sake
    Kokkos::deep_copy (y_h, y);
    output_execution_space::fence (); // for UVM's sake

    for (int i = 0; i < size; ++i) {
      const auto curEnt = y_h(i);
      TEST_EQUALITY( curEnt, static_cast<output_value_type> (i+1) );
    }
  }

  TEUCHOS_UNIT_TEST( TpetraUtils, copy_convert )
  {
    using Tpetra::Details::copyConvert;
    typedef Tpetra::Map<> map_type;
    typedef map_type::device_type default_device_type;

    auto comm = Tpetra::TestingUtilities::getDefaultComm ();
    // Create a Map just to ensure that Kokkos gets initialized and
    // finalized correctly.
    const map_type map (comm->getSize (), 1, 0, comm);

    // Test the case where the Views have the same value types and
    // memory spaces.
    {
      typedef int output_value_type;
      typedef int input_value_type;
      typedef default_device_type output_device_type;
      typedef Kokkos::View<output_value_type*, output_device_type> output_view_type;
      typedef default_device_type input_device_type;
      typedef Kokkos::View<input_value_type*, input_device_type> input_view_type;
      testCopyConvert<output_view_type, input_view_type> (success, out);
    }

    // Test the case where the Views have the same value types, but
    // different memory spaces, such that the output View's execution
    // space cannot access the source View's memory space.  This is
    // only possible in a CUDA build.
#ifdef KOKKOS_ENABLE_CUDA
    {
      typedef int output_value_type;
      typedef int input_value_type;
      // CudaUVMSpace is technically accessible from host, so use CudaSpace explicitly.
      typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> output_device_type;
      typedef Kokkos::View<output_value_type*, output_device_type> output_view_type;
      typedef typename output_view_type::HostMirror::execution_space host_execution_space;
      typedef Kokkos::Device<host_execution_space, Kokkos::HostSpace> input_device_type;
      typedef Kokkos::View<input_value_type*, input_device_type> input_view_type;
      testCopyConvert<output_view_type, input_view_type> (success, out);
    }
#endif // KOKKOS_ENABLE_CUDA

    // Test the case where the Views have different value types, but
    // the same memory spaces.
    {
      // Input and output Views need to have distinct types that
      // nevertheless can be converted one to the other (ignoring
      // overflow).
      typedef long output_value_type;
      typedef short input_value_type;
      typedef default_device_type output_device_type;
      typedef Kokkos::View<output_value_type*, output_device_type> output_view_type;
      typedef default_device_type input_device_type;
      typedef Kokkos::View<input_value_type*, input_device_type> input_view_type;
      testCopyConvert<output_view_type, input_view_type> (success, out);
    }

    // Test the case where the Views have different value types, and
    // different memory spaces, such that the output View's execution
    // space cannot access the source View's memory space.  This is
    // only possible in a CUDA build.
#ifdef KOKKOS_ENABLE_CUDA
    {
      // Input and output Views need to have distinct types that
      // nevertheless can be converted one to the other (ignoring
      // overflow).
      typedef long output_value_type;
      typedef short input_value_type;
      // CudaUVMSpace is technically accessible from host, so use CudaSpace explicitly.
      typedef Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace> output_device_type;
      typedef Kokkos::View<output_value_type*, output_device_type> output_view_type;
      typedef typename output_view_type::HostMirror::execution_space host_execution_space;
      typedef Kokkos::Device<host_execution_space, Kokkos::HostSpace> input_device_type;
      typedef Kokkos::View<input_value_type*, input_device_type> input_view_type;
      testCopyConvert<output_view_type, input_view_type> (success, out);
    }
#endif // KOKKOS_ENABLE_CUDA
  }

} // namespace (anonymous)


