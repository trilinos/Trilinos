#include <iostream>

#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_ArithTraitsTest.hpp"

using Kokkos::Details::ArithTraits;

const int N = 16;
const int blocksize = 16;

__global__ void hello (char *a, int *b) {
  //if (threadIdx.x < std::numeric_limits<int>::max ()) {
  if (threadIdx.x > ArithTraits<int>::min () &&
      threadIdx.x < ArithTraits<int>::max () &&
      threadIdx.x >= ArithTraits<int>::abs (0)) {
    a[threadIdx.x] += b[threadIdx.x];
  }
  //}
}

int main () {
  using std::cout;
  using std::endl;

  char a[N] = "Hello \0\0\0\0\0\0";
  int b[N] = {15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  char *ad;
  int *bd;
  const int csize = N*sizeof(char);
  const int isize = N*sizeof(int);

  //cout << "1 << 1 = " << (1 << 1) << endl;
  //cout << "1 << 5 = " << (1 << 5) << endl;
  //cout << "highestBit<int>(42) = " << highestBit (42) << endl;

  cout << "pow<char> (2, 3) = " << (int) ArithTraits<char>::pow (2, 3) << endl;
  cout << "pow<signed char> (2, 3) = " << (int) ArithTraits<signed char>::pow (2, 3) << endl;
  cout << "pow<signed char> (-2, 3) = " << (int) ArithTraits<signed char>::pow (-2, 3) << endl;
  cout << "pow<signed char> (2, -3) = " << (int) ArithTraits<signed char>::pow (2, -3) << endl;

  cout << "First: " << a << endl;

  cudaMalloc ((void**)&ad, csize);
  cudaMalloc ((void**)&bd, isize);
  cudaMemcpy (ad, a, csize, cudaMemcpyHostToDevice);
  cudaMemcpy (bd, b, isize, cudaMemcpyHostToDevice);

  dim3 dimBlock (blocksize, 1);
  dim3 dimGrid (1, 1);
  hello<<<dimGrid, dimBlock>>> (ad, bd);
  cudaMemcpy (a, ad, csize, cudaMemcpyDeviceToHost);
  cudaFree (ad);
  cudaFree (bd);

  cout << "Second: " << a << endl;

  // When printing char, remember to cast to int, because else the
  // number gets translated into an ASCII character.  (Small ASCII
  // values are control characters and are not easily printable.)

  bool success = true;
  bool localSuccess = true;

  localSuccess = ArithTraitsTester<char>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "char passed" << endl;
  } else {
    cout << "char FAILED" << endl;
  }

  // Interestingly enough, char and int8_t are different types, but
  // signed char and int8_t are the same (on my system).
  localSuccess = ArithTraitsTester<signed char>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "signed char passed" << endl;
  } else {
    cout << "signed char FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<unsigned char>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "unsigned char passed" << endl;
  } else {
    cout << "unsigned char FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<short>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "short passed" << endl;
  } else {
    cout << "short FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<unsigned short>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "unsigned short passed" << endl;
  } else {
    cout << "unsigned short FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<int>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "int passed" << endl;
  } else {
    cout << "int FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<unsigned int>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "unsigned int passed" << endl;
  } else {
    cout << "unsigned int FAILED" << endl;
  }

  // localSuccess = ArithTraitsTester<long>::test (cout);
  // success = success && localSuccess;
  // if (localSuccess) {
  //   cout << "long passed" << endl;
  // } else {
  //   cout << "long FAILED" << endl;
  // }

  // localSuccess = ArithTraitsTester<unsigned long>::test (cout);
  // success = success && localSuccess;
  // if (localSuccess) {
  //   cout << "unsigned long passed" << endl;
  // } else {
  //   cout << "unsigned long FAILED" << endl;
  // }

  localSuccess = ArithTraitsTester<int8_t>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "int8_t passed" << endl;
  } else {
    cout << "int8_t FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<uint8_t>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "uint8_t passed" << endl;
  } else {
    cout << "uint8_t FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<int16_t>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "int16_t passed" << endl;
  } else {
    cout << "int16_t FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<uint16_t>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "uint16_t passed" << endl;
  } else {
    cout << "uint16_t FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<int32_t>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "int32_t passed" << endl;
  } else {
    cout << "int32_t FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<uint32_t>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "uint32_t passed" << endl;
  } else {
    cout << "uint32_t FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<long>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "long passed" << endl;
  } else {
    cout << "long FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<unsigned long>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "unsigned long passed" << endl;
  } else {
    cout << "unsigned long FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<int64_t>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "int64_t passed" << endl;
  } else {
    cout << "int64_t FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<uint64_t>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "uint64_t passed" << endl;
  } else {
    cout << "uint64_t FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<long long>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "long long passed" << endl;
  } else {
    cout << "long long FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<unsigned long long>::test (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "unsigned long long passed" << endl;
  } else {
    cout << "unsigned long long FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<float>::test (cout) &&
    ArithTraitsTester<float>::testFloatingPoint (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "float passed" << endl;
  } else {
    cout << "float FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<double>::test (cout) &&
    ArithTraitsTester<double>::testFloatingPoint (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "double passed" << endl;
  } else {
    cout << "double FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<long double>::test (cout) &&
    ArithTraitsTester<long double>::testFloatingPoint (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "long double passed" << endl;
  } else {
    cout << "long double FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<std::complex<float> >::test (cout) &&
    ArithTraitsTester<std::complex<float> >::testFloatingPoint (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "std::complex<float> passed" << endl;
  } else {
    cout << "std::complex<float> FAILED" << endl;
  }

  localSuccess = ArithTraitsTester<std::complex<double> >::test (cout) &&
    ArithTraitsTester<std::complex<double> >::testFloatingPoint (cout);
  success = success && localSuccess;
  if (localSuccess) {
    cout << "std::complex<double> passed" << endl;
  } else {
    cout << "std::complex<double> FAILED" << endl;
  }

  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
