#include "GlobalReporter.hpp"
#include "Reporter.hpp"

namespace ngp_testing {
namespace global {
static constexpr int numReports = 5;
ReporterBase* hostReporter = nullptr;

ReporterBase* deviceReporterOnHost = nullptr;
NGP_TEST_DEVICE_ONLY ReporterBase* deviceReporterOnDevice = nullptr;
ReporterBase* deviceReporterAddress = nullptr;
}

using DeviceReporter = Reporter<Kokkos::DefaultExecutionSpace::device_type>;

#ifdef KOKKOS_ENABLE_OPENMP
using HostReporter = Reporter<Kokkos::OpenMP::device_type>;
#else
using HostReporter = Reporter<Kokkos::Serial::device_type>;
#endif

void copy_to_device(const DeviceReporter& reporter,
                    ReporterBase* const addr) {
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0,1), KOKKOS_LAMBDA(const int){
    global::deviceReporterOnDevice = addr;
    new (global::deviceReporterOnDevice) DeviceReporter(reporter);
  });
  Kokkos::fence();
}

void initialize_reporters() {
  global::hostReporter = new HostReporter(global::numReports);
  global::deviceReporterOnHost = new DeviceReporter(global::numReports);
  global::deviceReporterAddress = static_cast<DeviceReporter*>(Kokkos::kokkos_malloc(sizeof(DeviceReporter)));
  copy_to_device(dynamic_cast<DeviceReporter&>(*global::deviceReporterOnHost), global::deviceReporterAddress);
}

void finalize_reporters() {
  delete global::hostReporter;
  delete global::deviceReporterOnHost;
  Kokkos::kokkos_free(global::deviceReporterAddress);
}

NGP_TEST_FUNCTION ReporterBase* get_reporter() {
#ifdef __CUDA_ARCH__
  return global::deviceReporterOnDevice;
#else
  return global::hostReporter;
#endif
}

ReporterBase* get_device_reporter() {
  return global::deviceReporterOnHost;
}

}
