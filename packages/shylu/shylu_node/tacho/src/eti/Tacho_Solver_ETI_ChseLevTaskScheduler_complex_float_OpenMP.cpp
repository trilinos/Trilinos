#include "Kokkos_Core.hpp"

#include "Tacho_Solver.hpp"
#include "Tacho_Solver_Impl.hpp"

namespace Tacho {
#if defined(KOKKOS_ENABLE_OPENMP)
  using eti_scheduler_type = Kokkos::ChaseLevTaskScheduler<Kokkos::OpenMP>;
  template struct Solver<Kokkos::complex<float>,eti_scheduler_type>;
#endif
}
