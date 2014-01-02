//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER

#include <iostream>

#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_ArithTraitsTest.hpp"
#include "Kokkos_Serial.hpp"

#ifdef KOKKOS_HAVE_OPENMP
#  include "Kokkos_OpenMP.hpp"
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
#  include "Kokkos_Threads.hpp"
#endif // KOKKOS_HAVE_PTHREAD


template<class DeviceType>
bool runHostTests (std::ostream& out)
{
  bool success = true;

  //
  // Built-in char(acter) types
  //

  success = success && testArithTraitsOnHost<char, DeviceType> (out);
  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  success = success && testArithTraitsOnHost<signed char, DeviceType> (out);
  success = success && testArithTraitsOnHost<unsigned char, DeviceType> (out);

  //
  // Built-in integer types
  //

  success = success && testArithTraitsOnHost<short, DeviceType> (out);
  success = success && testArithTraitsOnHost<unsigned short, DeviceType> (out);
  success = success && testArithTraitsOnHost<int8_t, DeviceType> (out);
  success = success && testArithTraitsOnHost<uint8_t, DeviceType> (out);
  success = success && testArithTraitsOnHost<int16_t, DeviceType> (out);
  success = success && testArithTraitsOnHost<uint16_t, DeviceType> (out);
  success = success && testArithTraitsOnHost<int32_t, DeviceType> (out);
  success = success && testArithTraitsOnHost<uint32_t, DeviceType> (out);
  success = success && testArithTraitsOnHost<int, DeviceType> (out);
  success = success && testArithTraitsOnHost<unsigned int, DeviceType> (out);
  success = success && testArithTraitsOnHost<int64_t, DeviceType> (out);
  success = success && testArithTraitsOnHost<uint64_t, DeviceType> (out);
  success = success && testArithTraitsOnHost<long, DeviceType> (out);
  success = success && testArithTraitsOnHost<unsigned long, DeviceType> (out);
  success = success && testArithTraitsOnHost<long long, DeviceType> (out);
  success = success && testArithTraitsOnHost<unsigned long long, DeviceType> (out);

  //
  // Built-in real and complex floating-point types
  //

  success = success && testArithTraitsOnHost<float, DeviceType> (out);
  success = success && testArithTraitsOnHost<double, DeviceType> (out);
  success = success && testArithTraitsOnHost<long double, DeviceType> (out);
  success = success && testArithTraitsOnHost<std::complex<float>, DeviceType> (out);
  success = success && testArithTraitsOnHost<std::complex<double>, DeviceType> (out);
  success = success && testArithTraitsOnHost<std::complex<long double>, DeviceType> (out);

  return success;
}


template<class DeviceType>
bool runDeviceTests (std::ostream& out)
{
  bool success = true;

  //
  // Built-in char(acter) types
  //

  success = success && testArithTraitsOnDevice<char, DeviceType> (out);
  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  success = success && testArithTraitsOnDevice<signed char, DeviceType> (out);
  success = success && testArithTraitsOnDevice<unsigned char, DeviceType> (out);

  //
  // Built-in integer types
  //

  success = success && testArithTraitsOnDevice<short, DeviceType> (out);
  success = success && testArithTraitsOnDevice<unsigned short, DeviceType> (out);
  success = success && testArithTraitsOnDevice<int8_t, DeviceType> (out);
  success = success && testArithTraitsOnDevice<uint8_t, DeviceType> (out);
  success = success && testArithTraitsOnDevice<int16_t, DeviceType> (out);
  success = success && testArithTraitsOnDevice<uint16_t, DeviceType> (out);
  success = success && testArithTraitsOnDevice<int32_t, DeviceType> (out);
  success = success && testArithTraitsOnDevice<uint32_t, DeviceType> (out);
  success = success && testArithTraitsOnDevice<int, DeviceType> (out);
  success = success && testArithTraitsOnDevice<unsigned int, DeviceType> (out);
  success = success && testArithTraitsOnDevice<int64_t, DeviceType> (out);
  success = success && testArithTraitsOnDevice<uint64_t, DeviceType> (out);
  success = success && testArithTraitsOnDevice<long, DeviceType> (out);
  success = success && testArithTraitsOnDevice<unsigned long, DeviceType> (out);
  success = success && testArithTraitsOnDevice<long long, DeviceType> (out);
  success = success && testArithTraitsOnDevice<unsigned long long, DeviceType> (out);

  //
  // Built-in real floating-point types
  //

  success = success && testArithTraitsOnDevice<float, DeviceType> (out);
  success = success && testArithTraitsOnDevice<double, DeviceType> (out);

  return success;
}


int main () {
  using std::cout;
  using std::endl;

  bool success = true;

  bool serialSuccess = true;
  serialSuccess = serialSuccess && runHostTests<Kokkos::Serial> (cout);
  serialSuccess = serialSuccess && runDeviceTests<Kokkos::Serial> (cout);
  success = success && serialSuccess;
  if (serialSuccess) {
    cout << endl << "Kokkos::Serial host and device: TEST PASSED" << endl;
  } else {
    cout << endl << "Kokkos::Serial host and device: TEST FAILED" << endl;
  }

#ifdef KOKKOS_HAVE_OPENMP
  bool openmpSuccess = true;
  openmpSuccess = openmpSuccess && runHostTests<Kokkos::OpenMP> (cout);
  openmpSuccess = openmpSuccess && runDeviceTests<Kokkos::OpenMP> (cout);
  success = success && openmpSuccess;
  if (openmpSuccess) {
    cout << endl << "Kokkos::OpenMP host and device: TEST PASSED" << endl;
  } else {
    cout << endl << "Kokkos::OpenMP host and device: TEST FAILED" << endl;
  }
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
  bool threadsSuccess = true;
  threadsSuccess = threadsSuccess && runHostTests<Kokkos::Threads> (cout);
  threadsSuccess = threadsSuccess && runDeviceTests<Kokkos::Threads> (cout);
  success = success && threadsSuccess;
  if (threadsSuccess) {
    cout << endl << "Kokkos::Threads host and device: TEST PASSED" << endl;
  } else {
    cout << endl << "Kokkos::Threads host and device: TEST FAILED" << endl;
  }
#endif // KOKKOS_HAVE_PTHREAD

  if (success) {
    cout << endl << "End Result: TEST PASSED" << endl;
  } else {
    cout << endl << "End Result: TEST FAILED" << endl;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
