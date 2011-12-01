#include <gtest/gtest.h>
#include <Kokkos_Cuda.hpp>

namespace Test {

extern void test_device_cuda_init();

class cuda : public ::testing::Test {
  protected:
    static void SetUpTestCase() {
      Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
    }
    static void TearDownTestCase() {
      Kokkos::Cuda::finalize();
    }
};

extern void test_cuda_hexgrad(int exp_beg, int exp_end);
extern void test_cuda_gramschmidt(int exp_beg, int exp_end);

TEST_F( cuda, hexgrad ) {
  EXPECT_NO_THROW(test_cuda_hexgrad( 10, 20 ));
}

TEST_F( cuda, gramschmidt ) {
  EXPECT_NO_THROW(test_cuda_gramschmidt( 10, 20 ));
}

} // namespace Test

