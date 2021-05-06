#include "Kokkos_Core.hpp"

#include "Tacho.hpp"
#include "Tacho_Driver.hpp"
#include "Tacho_Driver_Impl.hpp"

namespace Tacho {
#if defined(KOKKOS_ENABLE_CUDA)
  using eti_value_type = Kokkos::complex<float>;
  using eti_device_type = typename UseThisDevice<Kokkos::Cuda>::type; 
  template struct Driver<eti_value_type,eti_device_type>;
#endif
}
