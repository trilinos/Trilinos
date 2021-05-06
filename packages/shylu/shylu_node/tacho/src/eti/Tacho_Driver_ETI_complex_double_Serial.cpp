#include "Kokkos_Core.hpp"

#include "Tacho.hpp"
#include "Tacho_Driver.hpp"
#include "Tacho_Driver_Impl.hpp"

namespace Tacho {
#if defined(KOKKOS_ENABLE_SERIAL)
  using eti_value_type = Kokkos::complex<double>;
  using eti_device_type = typename UseThisDevice<Kokkos::Serial>::type; 
  template struct Driver<eti_value_type,eti_device_type>;
#endif
}
