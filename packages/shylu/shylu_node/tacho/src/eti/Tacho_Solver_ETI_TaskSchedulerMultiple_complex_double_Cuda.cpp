#include "Kokkos_Core.hpp"

#include "Tacho_Solver.hpp"
#include "Tacho_Solver_Impl.hpp"

namespace Tacho {
#if defined(KOKKOS_ENABLE_CUDA)
  using eti_scheduler_type = Kokkos::TaskSchedulerMultiple<Kokkos::Cuda>; 
  template struct Solver<Kokkos::complex<double>,eti_scheduler_type>;
#endif
}
