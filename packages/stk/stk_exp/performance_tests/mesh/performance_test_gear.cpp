#ifndef __IBMCPP__
#include <gtest/gtest.h>
#include <sierra/mesh/fixture/gear.hpp>
#include <sierra/mesh/details/selected_buckets.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <iostream>

using namespace sierra::mesh;

TEST( gear, basic)
{
  double total_start_time = stk::cpu_time();
  {
    double start_time = stk::cpu_time();
    sierra::mesh::fixture::gear gear;
    double end_time = stk::cpu_time() - start_time;

    std::cout << "gear setup time: " << end_time << std::endl;

    start_time = stk::cpu_time();
    gear.generate();
    end_time = stk::cpu_time() - start_time;

    std::cout << "Generate time: " << end_time << std::endl;
  }
  double total_time = stk::cpu_time() - total_start_time;
  std::cout << "Total Time: " << total_time << std::endl;
}
#endif
