#ifndef _GlobalReporter_cpp_
#define _GlobalReporter_cpp_

#include "GlobalReporter.hpp"
#include "Reporter.hpp"

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

NGP_TEST_INLINE ReporterBase*& getDeviceReporterOnDevice()
{
  #ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HIP_GPU
    __device__ static ReporterBase* deviceReporterOnDevice = nullptr;
  #else
    static ReporterBase* deviceReporterOnDevice = nullptr;
  #endif
  return deviceReporterOnDevice;
}

inline
ReporterBase*& getDeviceReporterAddress()
{
  static ReporterBase* deviceReporterAddress = nullptr;
  return deviceReporterAddress;
}
}

#ifdef KOKKOS_ENABLE_CUDA
using DeviceReporter = Reporter<Kokkos::CudaHostPinnedSpace>;
#else
using DeviceReporter = Reporter<Kokkos::DefaultExecutionSpace::device_type>;
#endif

#ifdef KOKKOS_ENABLE_OPENMP
using HostReporter = Reporter<Kokkos::OpenMP::device_type>;
#else
using HostReporter = Reporter<Kokkos::Serial::device_type>;
#endif

inline
void copy_to_device(const DeviceReporter& reporter,
                    ReporterBase* const addr) {
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,1), KOKKOS_LAMBDA(const int){
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
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  return global::getDeviceReporterOnDevice();
#else
  return global::getHostReporter();
#endif
}

inline
ReporterBase* get_device_reporter() {
  return global::getDeviceReporterOnHost();
}

}

#endif

