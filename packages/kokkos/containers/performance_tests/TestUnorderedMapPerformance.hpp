#ifndef KOKKOS_TEST_UNORDERED_MAP_PERFORMANCE_HPP
#define KOKKOS_TEST_UNORDERED_MAP_PERFORMANCE_HPP

#include <impl/Kokkos_Timer.hpp>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>


namespace Perf {

template <typename Device, bool Near>
struct UnorderedMapTest
{
  typedef Device device_type;
  typedef Kokkos::UnorderedMap<uint32_t, uint32_t, device_type> map_type;
  typedef typename map_type::histogram_type histogram_type;

  uint32_t capacity;
  uint32_t inserts;
  uint32_t collisions;
  double   seconds;
  map_type map;
  histogram_type histogram;

  UnorderedMapTest( uint32_t arg_capacity, uint32_t arg_inserts, uint32_t arg_collisions)
    : capacity(arg_capacity)
    , inserts(arg_inserts)
    , collisions(arg_collisions)
    , seconds(0)
    , map(capacity)
    , histogram(map.get_histogram())
  {
    Kokkos::Impl::Timer wall_clock ;
    wall_clock.reset();

    Kokkos::parallel_for(inserts, *this);
    Device::fence();

    seconds = wall_clock.seconds();

    histogram.calculate();
    Device::fence();
  }

  void print(std::ostream & metrics_out, std::ostream & length_out, std::ostream & distance_out, std::ostream & block_distance_out)
  {
    metrics_out << capacity << " , ";
    metrics_out << inserts/collisions << " , ";
    metrics_out << (100.0 * inserts/collisions) / capacity << " , ";
    metrics_out << inserts << " , ";
    metrics_out << (map.has_failed_inserts() ? "true" : "false") << " , ";
    metrics_out << collisions << " , ";
    metrics_out << 1e9*(seconds/inserts) << " , ";
    metrics_out << seconds << std::endl;

    length_out << capacity << " , ";
    length_out << ((100.0 *inserts/collisions) / capacity) << " , ";
    length_out << collisions << " , ";
    histogram.print_length(length_out);

    distance_out << capacity << " , ";
    distance_out << ((100.0 *inserts/collisions) / capacity) << " , ";
    distance_out << collisions << " , ";
    histogram.print_distance(distance_out);

    block_distance_out << capacity << " , ";
    block_distance_out << ((100.0 *inserts/collisions) / capacity) << " , ";
    block_distance_out << collisions << " , ";
    histogram.print_block_distance(block_distance_out);
  }


  KOKKOS_INLINE_FUNCTION
  void operator()(uint32_t i) const
  {
    if (Near) {
      map.insert(i/collisions, i);
    }
    else {
      map.insert(i%(inserts/collisions), i);
    }
  }

};

//#define KOKKOS_COLLECT_UNORDERED_MAP_METRICS

template <typename Device, bool Near>
void run_performance_tests(std::string const & base_file_name)
{
#if defined(KOKKOS_COLLECT_UNORDERED_MAP_METRICS)
  std::string metrics_file_name = base_file_name + std::string("-metrics.csv");
  std::string length_file_name = base_file_name  + std::string("-length.csv");
  std::string distance_file_name = base_file_name + std::string("-distance.csv");
  std::string block_distance_file_name = base_file_name + std::string("-block_distance.csv");

  std::ofstream metrics_out( metrics_file_name.c_str(), std::ofstream::out );
  std::ofstream length_out( length_file_name.c_str(), std::ofstream::out );
  std::ofstream distance_out( distance_file_name.c_str(), std::ofstream::out );
  std::ofstream block_distance_out( block_distance_file_name.c_str(), std::ofstream::out );

  // set up file headers
  metrics_out << "Capacity , Unique , Percent Full , Attempted Inserts , Failed Inserts , Collision Ratio , Nanoseconds/Inserts, Seconds" << std::endl;
  length_out << "Capacity , Percent Full , ";
  distance_out << "Capacity , Percent Full , ";
  block_distance_out << "Capacity , Percent Full , ";

  for (int i=0; i<100; ++i) {
    length_out << i << " , ";
    distance_out << i << " , ";
    block_distance_out << i << " , ";
  }

  length_out << "\b\b\b   " << std::endl;
  distance_out << "\b\b\b   " << std::endl;
  block_distance_out << "\b\b\b   " << std::endl;

  Kokkos::Impl::Timer wall_clock ;
  for (uint32_t collisions = 1;  collisions <= 64u; collisions = collisions << 1) {
    wall_clock.reset();
    std::cout << "Collisions: " << collisions << std::endl;
    for (uint32_t i = 1; i <= 12; ++i) {
      std::cout << "  percent full (" << std::setprecision(3) << std::fixed << (100.0*i)/8 << ") " << std::flush;
      for (uint32_t capacity = 1<<12; capacity < 1<<24; capacity = capacity << 1) {
        uint32_t inserts = i*(capacity/8u);
        UnorderedMapTest<Device, Near> test(capacity, inserts*collisions, collisions);
        test.print(metrics_out, length_out, distance_out, block_distance_out);
        std::cout << capacity << ".." << std::flush;
      }
      std::cout << std::endl;
    }
    std::cout << "  " << wall_clock.seconds() << " secs" << std::endl;
  }
  metrics_out.close();
  length_out.close();
  distance_out.close();
  block_distance_out.close();
#else
  (void)base_file_name;
  std::cout << "skipping test" << std::endl;
#endif
}


} // namespace Perf

#endif //KOKKOS_TEST_UNORDERED_MAP_PERFORMANCE_HPP
