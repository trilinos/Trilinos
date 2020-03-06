#include "Kokkos_Core.hpp"

#include "Tacho.hpp"
#include "Tacho_Solver.hpp"
#include "Tacho_Solver_Impl.hpp"

namespace Tacho {
#if defined(KOKKOS_ENABLE_CUDA)
  using eti_scheduler_type = Kokkos::TaskScheduler<Kokkos::Cuda>;
  template struct Solver<float,eti_scheduler_type>;
#endif
}
