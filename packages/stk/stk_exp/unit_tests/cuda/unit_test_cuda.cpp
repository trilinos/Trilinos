#include <gtest/gtest.h>

#include <stk_util/stk_config.h>

// restrict this file to only build if cuda is enabled
#ifdef HAVE_STK_CUDA

#include <check_cuda.h>
#include <cuda_saxpy.h>

#include <iostream>

namespace {

TEST( stk_cuda, check_cuda )
{
  int expected_return_code = 0;
  EXPECT_EQ(expected_return_code, check_cuda());
}

TEST( stk_cuda, cuda_saxpy )
{
  const int num_times = 100;
  const size_t N = 1000000;
  const double alpha = 0.5;
  std::vector<float> x(N, 2.0), y(N, 1.0);

  int error = cuda_saxpy(num_times, N, &x[0], alpha, &y[0]);

  double expected_value = 101.0;
  EXPECT_EQ(expected_value, y[0]);
}

} // namespace

#endif

