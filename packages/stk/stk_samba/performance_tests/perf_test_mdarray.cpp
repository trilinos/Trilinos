#define ENABLE_MDARRAY_PERF_TEST false

#if ENABLE_MDARRAY_PERF_TEST

#include <gtest/gtest.h>

#include <samba/utility/mdarray.hpp>
#include <performance_tests/cpu_time.hpp>

#include <numeric>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {

//  | a b c |
//  | d e f |
//  | g h i |
//
//  det = aei + bfg + cdh - ceg -bdi -afh;

struct result
{
  double det;
  double fill_time;
  double compute_time;
};

template <size_t NumNodes>
result run_cstyle_mdarray(int seed)
{
  typedef double array[NumNodes][3][3];

  array a;

  std::srand(seed);

  result r;
  r.det = 0;

  double start_time = cpu_time();
  //setup random data
  for (size_t node=0; node<NumNodes; ++node) {
    for (size_t x=0; x<3u; ++x) {
    for (size_t y=0; y<3u; ++y) {
      a[node][x][y] = (std::rand() % 100) / 100.0;
    }}
  }

  r.fill_time = cpu_time() - start_time;

  start_time = cpu_time();

  for (size_t node=0; node<NumNodes; ++node)
  {
    r.det += a[node][0][0] * a[node][1][1] * a[node][2][2]
           + a[node][0][1] * a[node][1][2] * a[node][2][0]
           + a[node][0][2] * a[node][1][0] * a[node][2][1]
           - a[node][0][2] * a[node][1][1] * a[node][2][0]
           - a[node][0][1] * a[node][1][0] * a[node][2][2]
           - a[node][0][0] * a[node][1][2] * a[node][2][1];
  }

  r.compute_time = cpu_time() - start_time;

  return r;
}

template <size_t NumNodes>
result run_mdarray(int seed)
{
  typedef samba::mdarray<double[NumNodes][3][3]> array;

  array a;

  result r;
  r.det = 0;

  std::srand(seed);

  double start_time = cpu_time();
  //setup random data
  for (size_t node=0; node<NumNodes; ++node) {
    for (size_t x=0; x<3u; ++x) {
    for (size_t y=0; y<3u; ++y) {
      a[node][x][y] = (std::rand() % 100) / 100.0;
    }}
  }

  r.fill_time = cpu_time() - start_time;

  start_time = cpu_time();

  for (size_t node=0; node<NumNodes; ++node)
  {
    r.det += a[node][0][0] * a[node][1][1] * a[node][2][2]
           + a[node][0][1] * a[node][1][2] * a[node][2][0]
           + a[node][0][2] * a[node][1][0] * a[node][2][1]
           - a[node][0][2] * a[node][1][1] * a[node][2][0]
           - a[node][0][1] * a[node][1][0] * a[node][2][2]
           - a[node][0][0] * a[node][1][2] * a[node][2][1];
  }
  r.compute_time = cpu_time() - start_time;

  return r;
}


} // namespace

TEST(samba, mdarray_perf)
{
  std::cout << std::setprecision(3);

  const size_t num_iterations = 1000;

  result volatile c_result = {0,0,0};
  result volatile m_result = {0,0,0};

  const size_t num_nodes = 50000;


  {
    std::cout << "c-style mdarray: " << std::endl;

    for (size_t i=0; i<num_iterations; ++i) {
      result tmp = run_cstyle_mdarray<num_nodes>(i);
      c_result.det += tmp.det;
      c_result.fill_time += tmp.fill_time;
      c_result.compute_time += tmp.compute_time;
    }

    std::cout << "     fill: " << c_result.fill_time << " seconds" << std::endl;
    std::cout << "  compute: " << c_result.compute_time << " seconds" << std::endl;
    std::cout << "    total: " << c_result.fill_time + c_result.compute_time << " seconds" << std::endl << std::endl;

  }

  {
    std::cout << "mdarray: " << std::endl;

    for (size_t i=0; i<num_iterations; ++i) {
      result tmp = run_cstyle_mdarray<num_nodes>(i);
      m_result.det += tmp.det;
      m_result.fill_time += tmp.fill_time;
      m_result.compute_time += tmp.compute_time;
    }

    std::cout << "     fill: " << m_result.fill_time << " seconds" << std::endl;
    std::cout << "  compute: " << m_result.compute_time << " seconds" << std::endl;
    std::cout << "    total: " << m_result.fill_time + m_result.compute_time << " seconds" << std::endl << std::endl;
  }

  EXPECT_TRUE(c_result.det == m_result.det);
}

#endif

