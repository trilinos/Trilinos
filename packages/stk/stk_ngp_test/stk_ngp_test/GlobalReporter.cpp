#ifndef _GlobalReporter_cpp_
#define _GlobalReporter_cpp_

#include "GlobalReporter.hpp"
#include "Reporter.hpp"
#include "stk_util/ngp/NgpSpaces.hpp"

namespace ngp_testing {
namespace global {
static constexpr int numReports = 5;
inline
ReporterBase*& getHostReporter()
{
  static ReporterBase* hostReporter = nullptr;
  return hostReporter;
}

inline
ReporterBase*& getDeviceReporterOnHost()
{
  static ReporterBase* deviceReporterOnHost = nullptr;
  return deviceReporterOnHost;
}

#ifdef KOKKOS_ENABLE_SYCL
namespace{
sycl::ext::oneapi::experimental::device_global<
    ReporterBase*,
    decltype(sycl::ext::oneapi::experimental::properties(
        sycl::ext::oneapi::experimental::device_image_scope))>
  deviceReporterOnDevice;
}
#endif

NGP_TEST_INLINE ReporterBase*& getDeviceReporterOnDevice()
{
#ifndef KOKKOS_ENABLE_SYCL
  KOKKOS_IF_ON_DEVICE((
    __device__ static ReporterBase* deviceReporterOnDevice = nullptr;
    return deviceReporterOnDevice;
  ))
#else
  KOKKOS_IF_ON_DEVICE((
    return deviceReporterOnDevice;
  ))
#endif
  KOKKOS_IF_ON_HOST((
    static ReporterBase* deviceReporterOnDevice = nullptr;
    return deviceReporterOnDevice;
  ))
}

inline
ReporterBase*& getDeviceReporterAddress()
{
  static ReporterBase* deviceReporterAddress = nullptr;
  return deviceReporterAddress;
}
}

using DeviceReporter = Reporter<Kokkos::DefaultExecutionSpace::device_type>;
using HostReporter = Reporter<Kokkos::DefaultHostExecutionSpace::device_type>;

inline
void copy_to_device(const DeviceReporter& reporter,
                    ReporterBase* const addr) {
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

NGP_TEST_INLINE ReporterBase* get_reporter() {
  KOKKOS_IF_ON_DEVICE((
    return global::getDeviceReporterOnDevice();
  ))
  KOKKOS_IF_ON_HOST((
    return global::getHostReporter();
  ))
}

inline
ReporterBase* get_device_reporter() {
  return global::getDeviceReporterOnHost();
}

}

#endif

