#ifndef _GlobalReporter_cpp_
#define _GlobalReporter_cpp_

#include "GlobalReporter.hpp"
#include "Reporter.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

namespace ngp_testing {
namespace global {
static constexpr int numReports = 5;
inline
HostReporter*& getHostReporter()
{
  static HostReporter* hostReporter = nullptr;
  return hostReporter;
}

inline
DeviceReporter*& getDeviceReporterOnHost()
{
  static DeviceReporter* deviceReporterOnHost = nullptr;
  return deviceReporterOnHost;
}

#ifdef KOKKOS_ENABLE_SYCL
namespace{
  sycl::ext::oneapi::experimental::device_global<DeviceReporter*> deviceReporterOnDevice;
}
#endif

NGP_TEST_INLINE DeviceReporter*& getDeviceReporterOnDevice()
{
#ifndef KOKKOS_ENABLE_SYCL
  KOKKOS_IF_ON_DEVICE((
    __device__ static DeviceReporter* deviceReporterOnDevice = nullptr;
    return deviceReporterOnDevice;
  ))
#else
  KOKKOS_IF_ON_DEVICE((
    return deviceReporterOnDevice;
  ))
#endif
  KOKKOS_IF_ON_HOST((
    static DeviceReporter* deviceReporterOnDevice = nullptr;
    return deviceReporterOnDevice;
  ))
}

inline
DeviceReporter*& getDeviceReporterAddress()
{
  static DeviceReporter* deviceReporterAddress = nullptr;
  return deviceReporterAddress;
}
}

inline
void copy_to_device(const DeviceReporter& reporter, DeviceReporter* const addr) {
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int){
    global::getDeviceReporterOnDevice() = addr;
    new (global::getDeviceReporterOnDevice()) DeviceReporter(reporter);
  });
  Kokkos::fence();
}

inline
void initialize_reporters() {
  global::getHostReporter() = new HostReporter(global::numReports);
  global::getDeviceReporterOnHost() = new DeviceReporter(global::numReports);
  global::getDeviceReporterAddress() = static_cast<DeviceReporter*>(Kokkos::kokkos_malloc(sizeof(DeviceReporter)));
  copy_to_device(dynamic_cast<DeviceReporter&>(*global::getDeviceReporterOnHost()), global::getDeviceReporterAddress());
}

inline
void finalize_reporters() {
  delete global::getHostReporter();
  global::getHostReporter() = nullptr;
  delete global::getDeviceReporterOnHost();
  global::getDeviceReporterOnHost() = nullptr;
  Kokkos::kokkos_free(global::getDeviceReporterAddress());
}

NGP_TEST_INLINE HostReporter* get_host_reporter() {
  return global::getHostReporter();
}

NGP_TEST_INLINE DeviceReporter* get_device_reporter() {
  return global::getDeviceReporterOnDevice();
}

inline
DeviceReporter* get_device_reporter_on_host() {
  return global::getDeviceReporterOnHost();
}

}

#endif

