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

#include "Teuchos_UnitTestHarness.hpp"
#include "Kokkos_Sparse_impl_copyIntegers.hpp"
#include "Kokkos_ArithTraits.hpp"
#include <limits>
#include <type_traits>
#include <string>
#include <typeinfo>
#include <utility> // std::pair
#include <vector>

// CUDA 7.5 with GCC 4.8.4 doesn't always accept functors in an
// anonymous namespace.  If I name the namespace, it builds just fine.
namespace KokkosSparseTest {

template<class ViewType>
class Iota {
public:
  Iota (const ViewType& x, const bool negative = false) :
    x_ (x),
    negative_ (negative)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const typename ViewType::size_type& i) const
  {
    typedef typename ViewType::non_const_value_type val_type;
    typedef typename Kokkos::Details::ArithTraits<val_type>::mag_type mag_type;

    const val_type i_val = static_cast<val_type> (static_cast<mag_type> (i));
    x_(i) = negative_ ? -i_val : i_val;
  }

private:
  ViewType x_;
  bool negative_;
};

// Fill x with [0, 1, 2, ..., x.dimension_0() - 1],
// or the negative of that [0, -1, -2, ...] if negative is true.
template<class ViewType>
void
iota (const ViewType& x,
      const bool negative = false)
{
  static_assert (std::is_same<typename ViewType::value_type,
                   typename ViewType::non_const_value_type>::value,
                 "ViewType must be a nonconst View.");
  Kokkos::RangePolicy<typename ViewType::execution_space,
    typename ViewType::size_type> range (0, x.dimension_0 ());
  Kokkos::parallel_for (range, Iota<ViewType> (x, negative));
}

template<class XViewType,
         class YViewType>
class ViewsEqual {
public:
  typedef int value_type;
  ViewsEqual (const XViewType& x, const YViewType& y) : x_ (x), y_ (y) {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const typename XViewType::size_type& i, int& dst) const
  {
    if (x_(i) != y_(i)) {
      dst = 0;
    }
  }

  KOKKOS_INLINE_FUNCTION void init (int& dst) const {
    dst = 1;
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile int& dst, const volatile int& src) const
  {
    dst = dst && src;
  }

private:
  XViewType x_;
  YViewType y_;
};

template<class XViewType,
         class YViewType>
bool viewsEqual (const XViewType& x, const YViewType& y)
{
  Kokkos::RangePolicy<typename XViewType::execution_space,
    typename XViewType::size_type> range (0, x.dimension_0 ());
  int result = 1;
  Kokkos::parallel_reduce (range, ViewsEqual<XViewType, YViewType> (x, y), result);
  return result != 0;
}

} // namespace KokkosSparseTest

namespace { // (anonymous)

using KokkosSparseTest::iota;
using KokkosSparseTest::viewsEqual;
using std::endl;

template<class OutputIntegerType,
         class OutputDeviceType,
         class InputIntegerType,
         class InputDeviceType>
void
testCopyIntegersOneCase (bool& success,
                         Teuchos::FancyOStream& out,
                         const std::string& outputIntegerType,
                         const std::string& outputDeviceType,
                         const std::string& inputIntegerType,
                         const std::string& inputDeviceType)
{
  using KokkosSparse::Impl::copyIntegers;

  out << "OutputDeviceType: " << outputDeviceType
      << ", InputDeviceType: " << inputDeviceType
      << ", OutputIntegerType: " << outputIntegerType
      << ", InputIntegerType: " << inputIntegerType
      << endl;
  Teuchos::OSTab tab0 (out);
  // out << "Min OutputIntegerType: "
  //     << std::numeric_limits<OutputIntegerType>::min () << endl
  //     << "Min InputIntegerType: "
  //     << std::numeric_limits<InputIntegerType>::min () << endl;
  // out << "Max OutputIntegerType: "
  //     << std::numeric_limits<OutputIntegerType>::max () << endl
  //     << "Max InputIntegerType: "
  //     << std::numeric_limits<InputIntegerType>::max () << endl;
  out << "Max OutputIntegerType, cast to InputIntegerType: "
      << static_cast<InputIntegerType> (std::numeric_limits<OutputIntegerType>::max ())
      << endl;

  const size_t len = 10;

  // Both out and in are allocated, and out is nonconst.
  Kokkos::View<OutputIntegerType*, OutputDeviceType> outView ("outView", len);
  Kokkos::View<InputIntegerType*, InputDeviceType> inView ("inView", len);
  Kokkos::View<OutputIntegerType*, OutputDeviceType>
    outView_expected ("outView_expected", len);

  {
    out << "Test the no-overflow case: values in "
      "[0, 1, 2, ..., " << (len-1) << "]" << endl;
    Teuchos::OSTab tab1 (out);

    iota (inView); // fill out with [0, 1, 2, ..., len-1]
    iota (outView_expected);
    TEST_NOTHROW( copyIntegers (outView, inView) );
    TEST_ASSERT( viewsEqual (outView, outView_expected) );
  }

  // If possible given the two integer types, test the overflow case.
  constexpr bool inputBiggerThanOutput =
    sizeof (InputIntegerType) > sizeof (OutputIntegerType) ||
    (sizeof (InputIntegerType) == sizeof (OutputIntegerType) &&
     ! std::numeric_limits<InputIntegerType>::is_signed &&
     std::numeric_limits<OutputIntegerType>::is_signed);

  if (inputBiggerThanOutput) {
    out << "Test the overflow case: max input value is "
        << std::numeric_limits<InputIntegerType>::max ()
        << ", but max allowed output value is "
        << std::numeric_limits<OutputIntegerType>::max ()
        << endl;
    Teuchos::OSTab tab1 (out);

    iota (inView); // same as before, except...
    { // ... last entry gets a number too big to fit in InputIntegerType.
      auto inView_host = Kokkos::create_mirror_view (inView);
      inView_host[inView.dimension_0 () - 1] =
        std::numeric_limits<InputIntegerType>::max ();
      Kokkos::deep_copy (inView, inView_host);
    }
    // The function should catch overflow.
    TEST_THROW( copyIntegers (outView, inView), std::runtime_error );
  }

  if (std::is_signed<InputIntegerType>::value &&
      ! std::is_signed<OutputIntegerType>::value) {
    out << "Test that casts of negative signed numbers to unsigned also throw"
        << endl;
    Teuchos::OSTab tab1 (out);

    iota (inView, true);
    TEST_THROW( copyIntegers (outView, inView), std::runtime_error );
  }
}

template<class OutputIntegerType,
         class InputIntegerType>
void
testCopyIntegersAllDevices (bool& success,
                            Teuchos::FancyOStream& out,
                            const std::string& outputIntegerType,
                            const std::string& inputIntegerType,
                            const bool testOnlyDefaultDevice)
{
  std::string deviceType;
  std::string otherDeviceType;


#ifdef KOKKOS_HAVE_CUDA
  // Test Cuda and its host execution space, since both are initialized.
  typedef Kokkos::Cuda execution_space;
  typedef Kokkos::View<int*, execution_space>::HostMirror::execution_space other_execution_space;

  deviceType = "Device<Cuda>";
  otherDeviceType = "HostMirror";

#else // NOT KOKKOS_HAVE_CUDA

  // If we don't have Cuda, see if at least two different execution
  // spaces exist, one of which is Kokkos::Serial, and the other of
  // which is either Kokkos::OpenMP or Kokkos::Threads (not both!).

#  if defined(KOKKOS_HAVE_SERIAL)
  // If Serial exists, use it as the first execution space.  Try to
  // pick something different (if it exists) as the second space.
  typedef Kokkos::Serial execution_space;
  deviceType = "Device<Serial>";

#    if defined(KOKKOS_HAVE_OPENMP)
  typedef Kokkos::OpenMP other_execution_space;
  otherDeviceType = "Device<OpenMP>";
#    elif defined(KOKKOS_HAVE_PTHREAD)
  typedef Kokkos::Threads other_execution_space;
  otherDeviceType = "Device<Threads>";
#    else // ! KOKKOS_HAVE_OPENMP && ! KOKKOS_HAVE_PTHREAD
  // Serial is the only space that exists out of the four known spaces
  // (Cuda, OpenMP, Serial, Threads).  Other spaces may exist; one of
  // those could be the default, or Serial could be the default.  Just
  // use the default space; it's OK if it's the same as Serial.
  typedef Kokkos::DefaultExecutionSpace other_execution_space;
  otherDeviceType = "Default";

#    endif // KOKKOS_HAVE_OPENMP, KOKKOS_HAVE_PTHREAD
#  else // NOT defined(KOKKOS_HAVE_SERIAL)

  // Serial does NOT exist.  Key here is to make sure that we don't
  // use OpenMP and Threads at the same time.  It's OK if
  // execution_space and other_execution_space are the same.

#    if defined(KOKKOS_HAVE_OPENMP)
  typedef Kokkos::OpenMP execution_space;
  deviceType = "Device<OpenMP>";
  typedef Kokkos::OpenMP other_execution_space;
  otherDeviceType = "Device<OpenMP>";

#    elif defined(KOKKOS_HAVE_PTHREAD)

  typedef Kokkos::Threads execution_space;
  deviceType = "Device<Threads>";
  typedef Kokkos::Threads other_execution_space;
  otherDeviceType = "Device<Threads>";

#    else // ! KOKKOS_HAVE_OPENMP && ! KOKKOS_HAVE_PTHREAD

  // Hm, none of the four spaces (Cuda, OpenMP, Serial, Threads)
  // exist.  Hopefully some other space exists.
  typedef Kokkos::DefaultExecutionSpace execution_space;
  deviceType = "Default";
  typedef Kokkos::DefaultExecutionSpace other_execution_space;
  otherDeviceType = "Default";

#    endif // KOKKOS_HAVE_OPENMP, KOKKOS_HAVE_PTHREAD
#  endif // defined(KOKKOS_HAVE_SERIAL)
#endif // KOKKOS_HAVE_CUDA

  TEUCHOS_TEST_FOR_EXCEPTION
    (deviceType == "" || otherDeviceType == "", std::logic_error,
     "Failed to select two execution spaces for the test!");

  if (testOnlyDefaultDevice) {
    typedef Kokkos::Device<Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space> device_type;
    deviceType = "Default";
    testCopyIntegersOneCase<OutputIntegerType, device_type,
      InputIntegerType, device_type> (success, out, outputIntegerType,
                                      deviceType, inputIntegerType,
                                      deviceType);
    if (! std::is_same<OutputIntegerType, InputIntegerType>::value) {
      testCopyIntegersOneCase<InputIntegerType, device_type,
        OutputIntegerType, device_type> (success, out, inputIntegerType,
                                         deviceType, outputIntegerType,
                                         deviceType);
    }
  }
  else {
    typedef Kokkos::Device<execution_space, execution_space::memory_space> device_type;
    typedef Kokkos::Device<other_execution_space, other_execution_space::memory_space> other_device_type;

    testCopyIntegersOneCase<OutputIntegerType, device_type,
      InputIntegerType, other_device_type> (success, out, outputIntegerType,
                                            deviceType, inputIntegerType,
                                            otherDeviceType);
    if (! std::is_same<OutputIntegerType, InputIntegerType>::value) {
      testCopyIntegersOneCase<InputIntegerType, device_type,
        OutputIntegerType, other_device_type> (success, out, inputIntegerType,
                                               deviceType, outputIntegerType,
                                               otherDeviceType);
    }

    if (! std::is_same<device_type, other_device_type>::value) {
      testCopyIntegersOneCase<OutputIntegerType, other_device_type,
        InputIntegerType, device_type> (success, out, outputIntegerType,
                                        otherDeviceType, inputIntegerType,
                                        deviceType);
      if (! std::is_same<OutputIntegerType, InputIntegerType>::value) {
        testCopyIntegersOneCase<InputIntegerType, other_device_type,
          OutputIntegerType, device_type> (success, out, inputIntegerType,
                                           otherDeviceType, outputIntegerType,
                                           deviceType);

      }
    }
  }
}

void
testCopyIntegers (bool& success,
                  Teuchos::FancyOStream& out,
                  const bool testOnlyDefaultDevice)
{
  // output, input
  std::vector<std::pair<std::string, std::string> > failedTypePairs;

  {
    typedef short output_integer_type;
    typedef short input_integer_type;

    const std::string outputIntegerType = "short";
    const std::string inputIntegerType = "short";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef unsigned short output_integer_type;
    typedef short input_integer_type;

    const std::string outputIntegerType = "unsigned short";
    const std::string inputIntegerType = "short";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef short output_integer_type;
    typedef unsigned short input_integer_type;

    const std::string outputIntegerType = "short";
    const std::string inputIntegerType = "unsigned short";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef unsigned short output_integer_type;
    typedef unsigned short input_integer_type;

    const std::string outputIntegerType = "unsigned short";
    const std::string inputIntegerType = "unsigned short";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef int output_integer_type;
    typedef int input_integer_type;

    const std::string outputIntegerType = "int";
    const std::string inputIntegerType = "int";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef int output_integer_type;
    typedef unsigned int input_integer_type;

    const std::string outputIntegerType = "int";
    const std::string inputIntegerType = "unsigned int";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef unsigned int output_integer_type;
    typedef unsigned int input_integer_type;

    const std::string outputIntegerType = "unsigned int";
    const std::string inputIntegerType = "unsigned int";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef int output_integer_type;
    typedef long long input_integer_type;

    const std::string outputIntegerType = "int";
    const std::string inputIntegerType = "long long";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef int output_integer_type;
    typedef unsigned long long input_integer_type;

    const std::string outputIntegerType = "int";
    const std::string inputIntegerType = "unsigned long long";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef unsigned int output_integer_type;
    typedef long long input_integer_type;

    const std::string outputIntegerType = "unsigned int";
    const std::string inputIntegerType = "long long";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef unsigned int output_integer_type;
    typedef unsigned long long input_integer_type;

    const std::string outputIntegerType = "unsigned int";
    const std::string inputIntegerType = "unsigned long long";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef long long output_integer_type;
    typedef long long input_integer_type;

    const std::string outputIntegerType = "long long";
    const std::string inputIntegerType = "long long";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef long long output_integer_type;
    typedef unsigned long long input_integer_type;

    const std::string outputIntegerType = "long long";
    const std::string inputIntegerType = "unsigned long long";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  {
    typedef unsigned long long output_integer_type;
    typedef unsigned long long input_integer_type;

    const std::string outputIntegerType = "unsigned long long";
    const std::string inputIntegerType = "unsigned long long";
    bool lclSuccess = true;
    testCopyIntegersAllDevices<output_integer_type,
      input_integer_type> (lclSuccess, out, outputIntegerType, inputIntegerType, testOnlyDefaultDevice);
    if (! lclSuccess) {
      failedTypePairs.push_back (std::make_pair (outputIntegerType,
                                                 inputIntegerType));
    }
    success = success && lclSuccess;
  }

  if (failedTypePairs.size () != 0) {
    out << "The following (output type, input type) pairs FAILED:" << endl;
    Teuchos::OSTab tab1 (out);
    for (size_t k = 0; k < failedTypePairs.size (); ++k) {
      out << "(" << failedTypePairs[k].first << ","
          << failedTypePairs[k].second << ")" << endl;
    }
  }
}

} // namespace (anonymous)

int
main (int argc, char* argv[])
{
  using std::endl;

  Teuchos::RCP<Teuchos::FancyOStream> outPtr =
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cout));
  Teuchos::FancyOStream& out = *outPtr;

  out << "Call Kokkos::initialize" << endl;
  Kokkos::initialize (argc, argv);

  bool success = true;
  out << "Run test" << endl;
  const bool testOnlyDefaultDevice = false;
  testCopyIntegers (success, out, testOnlyDefaultDevice);

  out << "Call Kokkos::finalize" << endl;
  Kokkos::finalize ();

  if (success) {
    out << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  }
  else {
    out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}


//////////////////////////////////////////////////////////////////////



