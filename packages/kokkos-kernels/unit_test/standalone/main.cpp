#include<Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA
#include<Test_Cuda.hpp>
#endif
#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL
#include<Test_Serial.hpp>
#endif
#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP
#include<Test_OpenMP.hpp>
#endif
#ifdef KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS
#include<Test_Threads.hpp>
#endif

// Include your testfile here:
//#include<Test_Blas1_sum.hpp>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc,argv);
  {
    ::testing::InitGoogleTest( &argc, argv );
    int result = RUN_ALL_TESTS();
  }
  Kokkos::finalize();
}
