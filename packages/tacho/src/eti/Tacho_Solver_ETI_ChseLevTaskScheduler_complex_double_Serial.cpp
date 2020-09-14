#include "Kokkos_Core.hpp"

#include "Tacho_Solver.hpp"
#include "Tacho_Solver_Impl.hpp"

namespace Tacho {
#if defined(KOKKOS_ENABLE_SERIAL)
  using eti_scheduler_type = Kokkos::ChaseLevTaskScheduler<Kokkos::Serial>;
  template struct Solver<Kokkos::complex<double>,eti_scheduler_type>;
#endif
}
