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

  template<class OutputValueType,
           class OutputArrayLayout,
           class OutputDeviceType,
           class InputValueType,
           class InputArrayLayout,
           class InputDeviceType>
  void
  testRank1CopyConvert (bool& success, Teuchos::FancyOStream& out)
  {
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    using output_view_type =
      Kokkos::View<OutputValueType*, OutputArrayLayout, OutputDeviceType>;
    using input_view_type =
      Kokkos::View<InputValueType*, InputArrayLayout, InputDeviceType>;
    using output_value_type =
      typename output_view_type::non_const_value_type;
    using input_value_type =
      typename input_view_type::non_const_value_type;

    const int size = 5;

    input_view_type x (view_alloc ("x", WithoutInitializing), size);
    output_view_type y (view_alloc ("y", WithoutInitializing), size);
    Kokkos::fence (); // for UVM's sake

    auto x_h = Kokkos::create_mirror_view (x);
    for (int i = 0; i < size; ++i) {
      x_h(i) = static_cast<input_value_type> (i+1); // no entries are zero
    }
    Kokkos::deep_copy (x, x_h);

    Tpetra::Details::copyConvert (y, x);
    Kokkos::fence (); // for UVM's sake

    auto y_h = Kokkos::create_mirror_view (y);
    Kokkos::deep_copy (y_h, y);

    bool equal = true;
    for (int i = 0; i < size; ++i) {
      const output_value_type expected_val =
        static_cast<output_value_type> (i+1);
      const output_value_type actual_val = y_h(i);
      if (actual_val != expected_val) {
        equal = false;
        break;
      }
    }

    if (! equal) {
      for (int i = 0; i < size; ++i) {
        const output_value_type expected_val =
          static_cast<output_value_type> (i+1);
        const output_value_type actual_val = y_h(i);        
        TEST_EQUALITY( actual_val, expected_val );
      }
    }
  }

  template<class OutputValueType,
           class OutputArrayLayout,
           class OutputDeviceType,
           class InputValueType,
           class InputArrayLayout,
           class InputDeviceType>
  void
  testRank2CopyConvert (bool& success, Teuchos::FancyOStream& out)
  {
    using Kokkos::view_alloc;
    using Kokkos::WithoutInitializing;
    using output_view_type =
      Kokkos::View<OutputValueType**, OutputArrayLayout, OutputDeviceType>;
    using input_view_type =
      Kokkos::View<InputValueType**, InputArrayLayout, InputDeviceType>;
    using output_value_type =
      typename output_view_type::non_const_value_type;
    using input_value_type =
      typename input_view_type::non_const_value_type;

    const int numRows = 7;
    const int numCols = 11;

    input_view_type x (view_alloc ("x", WithoutInitializing),
                       numRows, numCols);
    output_view_type y (view_alloc ("y", WithoutInitializing),
                        numRows, numCols);
    Kokkos::fence (); // for UVM's sake

    auto x_h = Kokkos::create_mirror_view (x);
    for (int i = 0; i < numRows; ++i) {
      for (int j = 0; j < numCols; ++j) {
        const input_value_type val =
          static_cast<input_value_type> ((i+1) + numCols*(j+1));
        x_h(i,j) = val;
      }
    }
    Kokkos::deep_copy (x, x_h);

    Tpetra::Details::copyConvert (y, x);
    Kokkos::fence (); // for UVM's sake

    auto y_h = Kokkos::create_mirror_view (y);
    Kokkos::deep_copy (y_h, y);

    bool equal = true;
    for (int i = 0; i < numRows; ++i) {
      for (int j = 0; j < numCols; ++j) {
        const output_value_type expected_val =
          static_cast<output_value_type> ((i+1) + numCols*(j+1));
        const output_value_type actual_val = y_h(i,j);
        if (expected_val != actual_val) {
          equal = false;
          break;
        }
      }
      if (! equal) {
        break;
      }
    }

    if (! equal) {
      for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
          const output_value_type expected_val =
            static_cast<output_value_type> ((i+1) + numCols*(j+1));
          const output_value_type actual_val = y_h(i,j);
          TEST_EQUALITY( expected_val, actual_val );
        }
      }
    }
  }

  template<class OutputValueType,
           class OutputArrayLayout,
           class OutputDeviceType,
           class InputValueType,
           class InputArrayLayout,
           class InputDeviceType>
  void
  testCopyConvert (bool& success,
                   Teuchos::FancyOStream& out,
                   const char outputValueName[],
                   const char inputValueName[])
  {
    out << "OutputValueType: " << outputValueName << std::endl
        << "InputValueType: " << inputValueName << std::endl;
    Teuchos::OSTab tab1 (out);
    
    testRank1CopyConvert<OutputValueType, OutputArrayLayout, OutputDeviceType,
                         InputValueType, InputArrayLayout, InputDeviceType>
      (success, out);
    if (! success) {
      out << "Returning early" << std::endl;
      return;
    }
    testRank2CopyConvert<OutputValueType, OutputArrayLayout, OutputDeviceType,
                         InputValueType, InputArrayLayout, InputDeviceType>
      (success, out);
  }

  template<class OutputArrayLayout,
           class OutputDeviceType,
           class InputArrayLayout,
           class InputDeviceType>
  void
  test_for_all_value_types (bool& success,
                            Teuchos::FancyOStream& out,
                            const char outputArrayLayoutName[],
                            const char inputArrayLayoutName[])
  {
    out << "OutputArrayLayout: " << outputArrayLayoutName << std::endl
        << "InputArrayLayout: " << inputArrayLayoutName << std::endl;
    Teuchos::OSTab tab1 (out);

    {
      using output_value_type = int;
      using input_value_type = int;
      testCopyConvert<output_value_type, OutputArrayLayout, OutputDeviceType,
                      input_value_type, InputArrayLayout, InputDeviceType>
        (success, out, "int", "int");
      if (! success) {
        out << "Returning early" << std::endl;
        return;
      }
    }
    {
      using output_value_type = long long;
      using input_value_type = int;
      testCopyConvert<output_value_type, OutputArrayLayout, OutputDeviceType,
                      input_value_type, InputArrayLayout, InputDeviceType>
        (success, out, "long long", "int");
      if (! success) {
        out << "Returning early" << std::endl;
        return;
      }
    }
    {
      using output_value_type = int;
      using input_value_type = long long;
      testCopyConvert<output_value_type, OutputArrayLayout, OutputDeviceType,
                      input_value_type, InputArrayLayout, InputDeviceType>
        (success, out, "int", "long long");
      if (! success) {
        out << "Returning early" << std::endl;
        return;
      }
    }
    {
      using output_value_type = Kokkos::complex<float>;
      using input_value_type = double;
      testCopyConvert<output_value_type, OutputArrayLayout, OutputDeviceType,
                      input_value_type, InputArrayLayout, InputDeviceType>
        (success, out, "Kokkos::complex<float>", "double");
      if (! success) {
        out << "Returning early" << std::endl;
        return;
      }
    }
    {
      using output_value_type = double;
      using input_value_type = Kokkos::complex<float>;
      testCopyConvert<output_value_type, OutputArrayLayout, OutputDeviceType,
                      input_value_type, InputArrayLayout, InputDeviceType>
        (success, out, "double", "Kokkos::complex<float>");
      if (! success) {
        out << "Returning early" << std::endl;
        return;
      }
    }
  }

  template<class OutputDeviceType,
           class InputDeviceType>
  void
  test_for_all_value_types_and_layouts
  (bool& success, Teuchos::FancyOStream& out,
   const char outputDeviceName[],
   const char inputDeviceName[])
  {
    out << "OutputDeviceType: " << outputDeviceName << std::endl
        << "InputDeviceType: " << inputDeviceName << std::endl;
    Teuchos::OSTab tab1 (out);

    {
      using output_array_layout = Kokkos::LayoutLeft;
      using input_array_layout = Kokkos::LayoutLeft;
      test_for_all_value_types<output_array_layout,
                               OutputDeviceType,
                               input_array_layout,
                               InputDeviceType> (success, out,
                                                 "Left", "Left");
    }
    {
      using output_array_layout = Kokkos::LayoutRight;
      using input_array_layout = Kokkos::LayoutLeft;
      test_for_all_value_types<output_array_layout,
                               OutputDeviceType,
                               input_array_layout,
                               InputDeviceType> (success, out,
                                                 "Right", "Left");
    }
    {
      using output_array_layout = Kokkos::LayoutLeft;
      using input_array_layout = Kokkos::LayoutRight;
      test_for_all_value_types<output_array_layout,
                               OutputDeviceType,
                               input_array_layout,
                               InputDeviceType> (success, out,
                                                 "Left", "Right");
    }
    {
      using output_array_layout = Kokkos::LayoutRight;
      using input_array_layout = Kokkos::LayoutRight;
      test_for_all_value_types<output_array_layout,
                               OutputDeviceType,
                               input_array_layout,
                               InputDeviceType> (success, out,
                                                 "Right", "Right");
    }
  }

  void
  test_for_all_value_types_and_layouts_and_devices
  (bool& success, Teuchos::FancyOStream& out)
  {
    {
      using output_device_type =
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "HostSpace", "HostSpace");
    }

#ifdef KOKKOS_ENABLE_CUDA
    // CudaSpace
    {
      using output_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "CudaSpace", "CudaSpace");
    }
    {
      using output_device_type =
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "HostSpace", "CudaSpace");
    }
    {
      using output_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "CudaSpace", "HostSpace");
    }

    // Mixed CudaSpace and CudaUVMSpace
    {
      using output_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "CudaSpace", "CudaUVMSpace");
    }
    {
      using output_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "CudaUVMSpace", "CudaSpace");
    }

    // CudaUVMSpace
    {
      using output_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "CudaUVMSpace", "CudaUVMSpace");
    }
    {
      using output_device_type =
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "HostSpace", "CudaUVMSpace");
    }
    {
      using output_device_type =
        Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      using input_device_type =
        Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;
      test_for_all_value_types_and_layouts<output_device_type,
                                           input_device_type>
        (success, out, "CudaUVMSpace", "HostSpace");
    }
#endif // KOKKOS_ENABLE_CUDA
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( TpetraUtils, CopyConvert )
  {
    test_for_all_value_types_and_layouts_and_devices (success, out);
  }

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  // Initialize MPI (if enabled) before initializing Kokkos.  This
  // lets MPI control things like pinning processes to sockets.
  Teuchos::GlobalMPISession mpiSession (&argc, &argv);
  Kokkos::initialize (argc, argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  Kokkos::finalize ();
  return errCode;
}
