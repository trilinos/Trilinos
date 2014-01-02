#include <iostream>

#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_ArithTraitsTest.hpp"
#include "Kokkos_Cuda.hpp"

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

  success = success && runHostTests<Kokkos::Cuda> (cout);
  success = success && runDeviceTests<Kokkos::Cuda> (cout);

  if (success) {
    cout << endl << "End Result: TEST PASSED" << endl;
  } else {
    cout << endl << "End Result: TEST FAILED" << endl;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
