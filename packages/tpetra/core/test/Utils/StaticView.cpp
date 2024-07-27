// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_Details_StaticView.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <algorithm>
#include <sstream>
#include <type_traits>

namespace { // (anonymous)

  //
  // UNIT TESTS
  //

  template<class ViewType1, class ViewType2>
  bool
  view1dSame (ViewType1 x, ViewType2 y)
  {
    using execution_space = typename ViewType1::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, size_t>;

    if (x.extent (0) != y.extent (0)) {
      return false;
    }
    const size_t size = static_cast<size_t> (x.extent (0));
    int allSame = 0;
    Kokkos::parallel_reduce
      ("view1dSame", range_type (0, size),
       KOKKOS_LAMBDA (const size_t k, int& curResult) {
        const int equal = (x(k) == y(k)) ? 1 : 0;
        curResult = curResult && equal;
      }, Kokkos::LAnd<int> (allSame));
    return allSame == 1;
  }

  template<class ViewType1, class ViewType2>
  bool
  view2dSame (ViewType1 x, ViewType2 y)
  {
    using execution_space = typename ViewType1::execution_space;
    using range_type = Kokkos::RangePolicy<execution_space, size_t>;

    if (x.extent (0) != y.extent (0)) {
      return false;
    }
    if (x.extent (1) != y.extent (1)) {
      return false;
    }

    const size_t nrow = static_cast<size_t> (x.extent (0));
    const size_t ncol = static_cast<size_t> (x.extent (1));
    int allSame = 0;
    Kokkos::LAnd<int> reducer (allSame);
    Kokkos::parallel_reduce
      ("view2dSame", range_type (0, nrow),
       KOKKOS_LAMBDA (const size_t row, int& curResult) {
        int rowEqual = 1;
        for (size_t col = 0; col < ncol; ++col) {
          if (x(row,col) != y(row,col)) {
            rowEqual = 0;
            break;
          }
        }
        curResult = curResult && rowEqual;
      }, Kokkos::LAnd<int> (allSame));
    return allSame == 1;
  }

  template<class ValueType>
  KOKKOS_INLINE_FUNCTION ValueType toValue (const size_t k)
  {
    using mag_type = typename Kokkos::ArithTraits<ValueType>::mag_type;
    return static_cast<ValueType> (static_cast<mag_type> (k));
  }

  template<class ValueType, class DeviceType>
  void
  view1dIota (const Kokkos::View<ValueType*, DeviceType>& x,
              const ValueType startValue)
  {
    using exec_space = typename DeviceType::execution_space;
    using range_type = Kokkos::RangePolicy<exec_space>;

    Kokkos::parallel_for
      ("view1dIota", range_type (0, x.extent (0)),
       KOKKOS_LAMBDA (const size_t k) {
         x(k) = startValue + toValue<ValueType> (k);
       });
  }

  template<class ViewType>
  void
  view2dIota (const ViewType& x,
              const typename ViewType::non_const_value_type& startValue)
  {
    using value_type = typename ViewType::non_const_value_type;
    using exec_space = typename ViewType::execution_space;
    using range_type = Kokkos::RangePolicy<exec_space, size_t>;

    const size_t nrow = static_cast<size_t> (x.extent (0));
    const size_t ncol = static_cast<size_t> (x.extent (1));

    Kokkos::parallel_for
      ("view2dIota", range_type (0, nrow),
       KOKKOS_LAMBDA (const size_t row) {
         for (size_t col = 0; col < ncol; ++col) {
           x(row,col) = startValue +
             toValue<value_type> (ncol) * toValue<value_type> (row) +
             toValue<value_type> (row);
         }
       });
  }

  template<class ValueType, class DeviceType>
  void
  testStatic1dView (bool& success, Teuchos::FancyOStream& out,
                    const std::string& valueName,
                    const std::string& deviceName)
  {
    using view_type = Kokkos::View<ValueType*, DeviceType>;

    out << "Test ValueType=" << valueName
        << ", DeviceType=" << deviceName << std::endl;
    Teuchos::OSTab tab1 (out);

    const std::vector<size_t> sizes {{ 0, 1, 0, 5, 0, 17, 31, 0, 65, 11 }};
    const size_t maxSize = *std::max_element (sizes.begin (), sizes.end ());
    const ValueType initVal = static_cast<ValueType> (1.0);

    view_type referenceView
      (Kokkos::view_alloc ("Reference View", Kokkos::WithoutInitializing),
       maxSize);
    view1dIota (referenceView, initVal);

    for (size_t size : sizes) {
      view_type view =
        Tpetra::Details::getStatic1dView<ValueType, DeviceType> (size);
      TEST_ASSERT( size_t (view.extent (0)) == size );
      if (size != 0) {
        const ValueType* rawPtr = view.data ();
        TEST_ASSERT( rawPtr != nullptr );
        TEST_ASSERT( reinterpret_cast<size_t> (rawPtr) % sizeof (ValueType) == 0 );
        view_type refView =
          Kokkos::subview (referenceView,
                           std::pair<size_t, size_t> (0, size));
        TEST_ASSERT( refView.data () != rawPtr );

        view1dIota (view, initVal);
        TEST_ASSERT( view1dSame (view, refView) );
      }
      TEST_NOTHROW( out << "View label: " << view.label () );
    }
  }

  template<class DeviceType>
  void
  testStatic1dViewValues (bool& success, Teuchos::FancyOStream& out,
                          const std::string& deviceName)
  {
    testStatic1dView<double, DeviceType>
      (success, out, "double", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic1dView<float, DeviceType>
      (success, out, "float", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic1dView<Kokkos::complex<double>, DeviceType>
      (success, out, "complex<double>", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic1dView<Kokkos::complex<float>, DeviceType>
      (success, out, "complex<float>", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic1dView<size_t, DeviceType>
      (success, out, "size_t", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic1dView<int, DeviceType>
      (success, out, "int", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic1dView<long long, DeviceType>
      (success, out, "long long", deviceName);
  }

  void
  testStatic1dViewValuesAndDevices (bool& success, Teuchos::FancyOStream& out)
  {
#ifdef KOKKOS_ENABLE_CUDA
    {
      using device_type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      testStatic1dViewValues<device_type>
        (success, out, "(Cuda,CudaUVMSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
    {
      using device_type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      testStatic1dViewValues<device_type>
        (success, out, "(Cuda,CudaSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
    {
      using device_type = Kokkos::Device<Kokkos::HostSpace::execution_space,
                                         Kokkos::CudaHostPinnedSpace>;
      testStatic1dViewValues<device_type>
        (success, out, "(Host,CudaHostPinnedSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_OPENMP
    {
      using device_type = Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>;
      testStatic1dViewValues<device_type>
        (success, out, "(OpenMP,HostSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_THREADS
    {
      using device_type = Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>;
      testStatic1dViewValues<device_type>
        (success, out, "(Threads,HostSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_THREADS

#ifdef KOKKOS_ENABLE_SERIAL
    {
      using device_type = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
      testStatic1dViewValues<device_type>
        (success, out, "(Serial,HostSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_SERIAL
  }

  TEUCHOS_UNIT_TEST( StaticView, Get1dView )
  {
    out << "Test getStatic1dView" << std::endl;
    Teuchos::OSTab tab1 (out);
    testStatic1dViewValuesAndDevices (success, out);
  }

  template<class ValueType, class DeviceType>
  void
  testStatic2dView (bool& success, Teuchos::FancyOStream& out,
                    const std::string& valueName,
                    const std::string& deviceName)
  {
    using view_type = Kokkos::View<ValueType**, Kokkos::LayoutLeft, DeviceType>;

    out << "Test ValueType=" << valueName
        << ", DeviceType=" << deviceName << std::endl;
    Teuchos::OSTab tab1 (out);

    const std::vector<size_t> row_sizes {{ 0, 1, 0, 5, 0, 17, 31, 0, 65, 11 }};
    const std::vector<size_t> col_sizes {{ 0, 1, 2, 3, 5 }};

    const size_t max_row_size =
      *std::max_element (row_sizes.begin (), row_sizes.end ());
    const size_t max_col_size =
      *std::max_element (col_sizes.begin (), col_sizes.end ());

    const ValueType initVal = static_cast<ValueType> (1.0);

    view_type referenceView
      (Kokkos::view_alloc ("Reference View", Kokkos::WithoutInitializing),
       max_row_size, max_col_size);

    for (size_t row_size : row_sizes) {
      for (size_t col_size : col_sizes) {
        view_type view =
          Tpetra::Details::getStatic2dView<ValueType, DeviceType>
          (row_size, col_size);
        TEST_ASSERT( size_t (view.extent (0)) == row_size );
        TEST_ASSERT( size_t (view.extent (1)) == col_size );

        if (row_size != 0 && col_size != 0) {
          const ValueType* rawPtr = view.data ();
          TEST_ASSERT( rawPtr != nullptr );
          TEST_ASSERT( reinterpret_cast<size_t> (rawPtr) % sizeof (ValueType) == 0 );

          auto refView =
            Kokkos::subview (referenceView,
                             std::pair<size_t, size_t> (0, row_size),
                             std::pair<size_t, size_t> (0, col_size));
          TEST_ASSERT( refView.data () != rawPtr );

          if (! success) {
            out << "Stopping test early" << std::endl;
            return;
          }
          view2dIota (view, initVal);
          view2dIota (refView, initVal);
          TEST_ASSERT( view2dSame (view, refView) );
        }
        TEST_NOTHROW( out << "View label: " << view.label () );
      }
    }
  }

  template<class DeviceType>
  void
  testStatic2dViewValues (bool& success, Teuchos::FancyOStream& out,
                          const std::string& deviceName)
  {
    testStatic2dView<double, DeviceType>
      (success, out, "double", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dView<float, DeviceType>
      (success, out, "float", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dView<Kokkos::complex<double>, DeviceType>
      (success, out, "complex<double>", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dView<Kokkos::complex<float>, DeviceType>
      (success, out, "complex<float>", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dView<size_t, DeviceType>
      (success, out, "size_t", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dView<int, DeviceType>
      (success, out, "int", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dView<long long, DeviceType>
      (success, out, "long long", deviceName);
  }

  void
  testStatic2dViewValuesAndDevices (bool& success, Teuchos::FancyOStream& out)
  {
#ifdef KOKKOS_ENABLE_CUDA
    {
      using device_type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      testStatic2dViewValues<device_type>
        (success, out, "(Cuda,CudaUVMSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
    {
      using device_type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      testStatic2dViewValues<device_type>
        (success, out, "(Cuda,CudaSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
    {
      using device_type = Kokkos::Device<Kokkos::HostSpace::execution_space,
                                         Kokkos::CudaHostPinnedSpace>;
      testStatic2dViewValues<device_type>
        (success, out, "(Host,CudaHostPinnedSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_OPENMP
    {
      using device_type = Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>;
      testStatic2dViewValues<device_type>
        (success, out, "(OpenMP,HostSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_THREADS
    {
      using device_type = Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>;
      testStatic2dViewValues<device_type>
        (success, out, "(Threads,HostSpace)");
    }
#endif // KOKKOS_ENABLE_THREADS

#ifdef KOKKOS_ENABLE_SERIAL
    {
      using device_type = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
      testStatic2dViewValues<device_type>
        (success, out, "(Serial,HostSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_SERIAL
  }

  TEUCHOS_UNIT_TEST( StaticView, Get2dView )
  {
    out << "Test getStatic2dView" << std::endl;
    Teuchos::OSTab tab1 (out);
    TEST_NOTHROW( testStatic2dViewValuesAndDevices (success, out) );
  }

  template<class ValueType, class DeviceType>
  void
  testStatic2dDualView (bool& success, Teuchos::FancyOStream& out,
                        const std::string& valueName,
                        const std::string& deviceName)
  {
    using dual_view_type =
      Kokkos::DualView<ValueType**, Kokkos::LayoutLeft, DeviceType>;

    out << "Test ValueType=" << valueName
        << ", DeviceType=" << deviceName << std::endl;
    Teuchos::OSTab tab1 (out);

    const std::vector<size_t> row_sizes {{ 0, 1, 0, 5, 0, 17, 31, 0, 65, 11 }};
    const std::vector<size_t> col_sizes {{ 0, 1, 2, 3, 5 }};

    const size_t max_row_size =
      *std::max_element (row_sizes.begin (), row_sizes.end ());
    const size_t max_col_size =
      *std::max_element (col_sizes.begin (), col_sizes.end ());

    const ValueType initVal = static_cast<ValueType> (1.0);
    dual_view_type referenceView ("Reference View", max_row_size, max_col_size);

    for (size_t row_size : row_sizes) {
      for (size_t col_size : col_sizes) {
        dual_view_type view =
          Tpetra::Details::getStatic2dDualView<ValueType, DeviceType>
          (row_size, col_size);
        TEST_ASSERT( size_t (view.extent (0)) == row_size );
        TEST_ASSERT( size_t (view.extent (1)) == col_size );
        TEST_ASSERT( size_t (view.d_view.extent (0)) == row_size );
        TEST_ASSERT( size_t (view.d_view.extent (1)) == col_size );
        TEST_ASSERT( size_t (view.h_view.extent (0)) == row_size );
        TEST_ASSERT( size_t (view.h_view.extent (1)) == col_size );
        TEST_ASSERT( ! (view.need_sync_device () && view.need_sync_host ()) );

        if (row_size != 0 && col_size != 0) {
          const ValueType* rawPtr_d = view.d_view.data ();
          TEST_ASSERT( rawPtr_d != nullptr );
          TEST_ASSERT( reinterpret_cast<size_t> (rawPtr_d) %
                       sizeof (ValueType) == 0 );

          const ValueType* rawPtr_h = view.h_view.data ();
          TEST_ASSERT( rawPtr_h != nullptr );
          TEST_ASSERT( reinterpret_cast<size_t> (rawPtr_h) %
                       sizeof (ValueType) == 0 );

          // Result of getStatic2dDualView must behave as if the host
          // View were created by Kokkos::create_mirror_view.
          if (std::is_same<
                typename dual_view_type::t_dev::memory_space,
                typename dual_view_type::t_host::memory_space>::value) {
            TEST_ASSERT( rawPtr_d == rawPtr_h );
          }

          auto refView =
            Kokkos::subview (referenceView,
                             std::pair<size_t, size_t> (0, row_size),
                             std::pair<size_t, size_t> (0, col_size));
          TEST_ASSERT( refView.d_view.data () != rawPtr_d );
          TEST_ASSERT( refView.h_view.data () != rawPtr_h );

          refView.clear_sync_state ();
          view.clear_sync_state ();

          // mfh 26 Mar 2019: Work-around for Kokkos bug:
          // https://github.com/kokkos/kokkos/issues/2051
#ifdef KOKKOS_ENABLE_CUDA
          constexpr bool ok_to_test_sync =
            ! std::is_same<typename DeviceType::memory_space,
                Kokkos::CudaSpace>::value;
#else
          constexpr bool ok_to_test_sync = true;
#endif // KOKKOS_ENABLE_CUDA

          if (ok_to_test_sync) {
            refView.modify_host ();
            view.modify_host ();
            view2dIota (view.h_view, initVal);
            view2dIota (refView.h_view, initVal);
            TEST_ASSERT( view2dSame (view.h_view, refView.h_view) );

            refView.sync_device ();
            view.sync_device ();

            refView.modify_device ();
            view.modify_device ();
            view2dIota (view.d_view, initVal);
            view2dIota (refView.d_view, initVal);
            TEST_ASSERT( view2dSame (view.d_view, refView.d_view) );
          }
          else {
            refView.modify_device ();
            view.modify_device ();
            view2dIota (view.d_view, initVal);
            view2dIota (refView.d_view, initVal);
            TEST_ASSERT( view2dSame (view.d_view, refView.d_view) );
          }
        }

        TEST_NOTHROW( out << "Device View label: " << view.d_view.label () );
        TEST_NOTHROW( out << "Host View label: " << view.h_view.label () );
      }
    }
  }

  template<class DeviceType>
  void
  testStatic2dDualViewValues (bool& success, Teuchos::FancyOStream& out,
                              const std::string& deviceName)
  {
    testStatic2dDualView<double, DeviceType>
      (success, out, "double", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dDualView<float, DeviceType>
      (success, out, "float", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dDualView<Kokkos::complex<double>, DeviceType>
      (success, out, "complex<double>", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dDualView<Kokkos::complex<float>, DeviceType>
      (success, out, "complex<float>", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dDualView<size_t, DeviceType>
      (success, out, "size_t", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dDualView<int, DeviceType>
      (success, out, "int", deviceName);
    if (! success) {
      out << "Stopping test early" << std::endl;
      return;
    }

    testStatic2dDualView<long long, DeviceType>
      (success, out, "long long", deviceName);
  }

  void
  testStatic2dDualViewValuesAndDevices (bool& success, Teuchos::FancyOStream& out)
  {
#ifdef KOKKOS_ENABLE_CUDA
    {
      using device_type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>;
      testStatic2dDualViewValues<device_type>
        (success, out, "(Cuda,CudaUVMSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
    {
      using device_type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;
      testStatic2dDualViewValues<device_type>
        (success, out, "(Cuda,CudaSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
    {
      using device_type = Kokkos::Device<Kokkos::HostSpace::execution_space,
                                         Kokkos::CudaHostPinnedSpace>;
      testStatic2dDualViewValues<device_type>
        (success, out, "(Host,CudaHostPinnedSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_OPENMP
    {
      using device_type = Kokkos::Device<Kokkos::OpenMP, Kokkos::HostSpace>;
      testStatic2dDualViewValues<device_type>
        (success, out, "(OpenMP,HostSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_THREADS
    {
      using device_type = Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>;
      testStatic2dDualViewValues<device_type>
        (success, out, "(Threads,HostSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_THREADS

#ifdef KOKKOS_ENABLE_SERIAL
    {
      using device_type = Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>;
      testStatic2dDualViewValues<device_type>
        (success, out, "(Serial,HostSpace)");
      if (! success) {
        out << "Stopping test early" << std::endl;
        return;
      }
    }
#endif // KOKKOS_ENABLE_SERIAL
  }

  TEUCHOS_UNIT_TEST( StaticView, Get2dDualView )
  {
    out << "Test getStatic2dDualView" << std::endl;
    Teuchos::OSTab tab1 (out);
    testStatic2dDualViewValuesAndDevices (success, out);
  }

} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  // Initialize MPI (if enabled) before initializing Kokkos.  This
  // lets MPI control things like pinning processes to sockets.
  Teuchos::GlobalMPISession mpiScope (&argc, &argv);
  Kokkos::ScopeGuard kokkosScope (argc, argv);
  const int errCode =
    Teuchos::UnitTestRepository::runUnitTestsFromMain (argc, argv);
  return errCode;
}
