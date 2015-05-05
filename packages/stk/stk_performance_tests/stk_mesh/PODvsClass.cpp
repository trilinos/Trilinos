#include <iostream>
#include <vector>
#include <iomanip>
#include <stdint.h>
#include <type_traits>
#include <sys/time.h>
#include <sys/resource.h>
#include <gtest/gtest.h>

struct TimingData
{
  double defaultInitTime;
  double pushBackInitTime;
  double copyTime;
};

double
cpu_time()
{
  struct rusage my_rusage;

  ::getrusage(RUSAGE_SELF, &my_rusage);

  double seconds = my_rusage.ru_utime.tv_sec + my_rusage.ru_stime.tv_sec;
  double micro_seconds = my_rusage.ru_utime.tv_usec + my_rusage.ru_stime.tv_usec;

  return seconds + micro_seconds*1.0e-6;
}

struct PodThingy
{
  enum stuff_t {
    InvalidThingy = 0ull,
    MinThingy = 1ull,
    MaxThingy = ~0ull
  };

  uint64_t m_value;

  PodThingy operator=(stuff_t val) { m_value = val; return *this; }
};

struct NonPodThingy
{
  enum stuff_t {
    InvalidThingy = 0ull,
    MinThingy = 1ull,
    MaxThingy = ~0ull
  };

  uint64_t m_value;

  NonPodThingy() : m_value(InvalidThingy) {};
  NonPodThingy(uint64_t value) : m_value(value) {};
  NonPodThingy operator=(stuff_t val) { m_value = val; return *this; }
};

template<class THINGY>
TimingData time_thingy(size_t num_thingies, size_t num_repeat)
{
  TimingData timingData;
  {
    double time_start = cpu_time();
    for (size_t i = 0; i < num_repeat; ++i) {
      std::vector<THINGY> thingies(num_thingies);
    }
    double time_stop = cpu_time();

    timingData.defaultInitTime = time_stop - time_start;
  }

  {
    double time_start = cpu_time();
    for (size_t i = 0; i < num_repeat; ++i) {
      std::vector<THINGY> thingies;
      thingies.reserve(num_thingies);
      for (size_t i = 0; i < num_thingies; ++i) {
        thingies.push_back(THINGY());
      }
    }
    double time_stop = cpu_time();

    timingData.pushBackInitTime = time_stop - time_start;

    std::vector<THINGY> thingies;
    thingies.reserve(num_thingies);
    for (size_t i = 0; i < num_thingies; ++i) {
      thingies.push_back(THINGY());
    }

    time_start = cpu_time();
    for (size_t i = 0; i < num_repeat; ++i) {
      std::vector<THINGY> other_thingies(thingies);
    }
    time_stop = cpu_time();

    timingData.copyTime = time_stop - time_start;
  }

  return timingData;
}

void print_timings(const std::string & typeName, const TimingData & newTiming, const TimingData & referenceTiming)
{
  std::cout << std::setprecision(3) << std::fixed;
  std::cout << "--------------------------------------------------------" << std::endl;
  std::cout << typeName << ":      default init = " << newTiming.defaultInitTime << " s";
  if (&newTiming == &referenceTiming) {
    std::cout << std::endl;
  } else {
    std::cout << "  (" << std::showpos << (newTiming.defaultInitTime - referenceTiming.defaultInitTime) * 100. / referenceTiming.defaultInitTime
              << "%)" << std::noshowpos << std::endl;
  }

  std::cout << typeName << ":    push_back init = " << newTiming.pushBackInitTime << " s";
  if (&newTiming == &referenceTiming) {
    std::cout << std::endl;
  } else {
    std::cout << "  (" << std::showpos << (newTiming.pushBackInitTime - referenceTiming.pushBackInitTime) * 100. / referenceTiming.pushBackInitTime
              << "%)" << std::noshowpos << std::endl;
  }

  std::cout << typeName << ": copy construction = " << newTiming.copyTime << " s";
  if (&newTiming == &referenceTiming) {
    std::cout << std::endl;
  } else {
    std::cout << "  (" << std::showpos << (newTiming.copyTime - referenceTiming.copyTime) * 100. / referenceTiming.copyTime
              << "%)" << std::noshowpos << std::endl;
  }
}


TEST(POD, performanceOfPODvsClass)
{
  static_assert(std::is_pod<PodThingy>::value, "PodThingy should be POD");
  static_assert(!std::is_pod<NonPodThingy>::value, "NonPodThingy should not be POD");

  const size_t NUM_THINGIES = 1000000;
  const size_t NUM_REPEAT = 500;

  time_thingy<uint64_t>(NUM_THINGIES, NUM_REPEAT);
  const TimingData timingData_uint64 = time_thingy<uint64_t>(NUM_THINGIES, NUM_REPEAT);
  print_timings("uint64_t", timingData_uint64, timingData_uint64);

  time_thingy<PodThingy>(NUM_THINGIES, NUM_REPEAT);
  const TimingData timingData_PodThingy = time_thingy<PodThingy>(NUM_THINGIES, NUM_REPEAT);
  print_timings("PodThingy", timingData_PodThingy, timingData_uint64);

  time_thingy<NonPodThingy>(NUM_THINGIES, NUM_REPEAT);
  const TimingData timingData_NonPodThingy = time_thingy<NonPodThingy>(NUM_THINGIES, NUM_REPEAT);
  print_timings("NonPodThingy", timingData_NonPodThingy, timingData_uint64);
}

