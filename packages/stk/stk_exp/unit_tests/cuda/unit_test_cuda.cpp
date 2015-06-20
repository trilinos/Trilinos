#include <gtest/gtest.h>

#include <stk_util/stk_config.h>

// restrict this file to only build if cuda is enabled
#ifdef HAVE_STK_CUDA

#include <check_cuda.h>
#include <cuda_memory.h>
#include <cuda_saxpy.h>

#include <iostream>

namespace {

TEST( stk_exp_cuda, check_cuda )
{
  int expected_return_code = 0;
  EXPECT_EQ(expected_return_code, check_cuda());
}

TEST( stk_exp_cuda, cuda_saxpy )
{
  const int num_times = 100;
  const size_t N = 1000000;
  const double alpha = 0.5;

  float* x = cuda_alloc_float(N, 2.0);
  float* y = cuda_alloc_float(N, 1.0);

  int err = cuda_saxpy(num_times, N, x, alpha, y);
  EXPECT_EQ(0, err);

  std::vector<float> host_y;
  copy_cuda_memory_to_host(N, y, host_y);

  double expected_value = 101.0;
  EXPECT_EQ(expected_value, host_y[0]);
}

} // namespace

#endif

