#ifndef STK_NGP_TEST_GLOBAL_REPORTER_HPP
#define STK_NGP_TEST_GLOBAL_REPORTER_HPP

#include "NgpTestDeviceMacros.hpp"

namespace ngp_testing {

template <typename DeviceType>
class Reporter;

using DeviceReporter = Reporter<Kokkos::DefaultExecutionSpace::device_type>;
using HostReporter = Reporter<Kokkos::DefaultHostExecutionSpace::device_type>;

void initialize_reporters();
void finalize_reporters();

NGP_TEST_FUNCTION HostReporter*  get_host_reporter();
NGP_TEST_FUNCTION DeviceReporter* get_device_reporter();
DeviceReporter* get_device_reporter_on_host();

}

#include "GlobalReporter.cpp"

#endif
