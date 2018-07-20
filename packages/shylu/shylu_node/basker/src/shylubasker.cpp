
#include "shylubasker_types.hpp"
#include "shylubasker_decl.hpp"
//#include "shylubasker_def.hpp"

#include <Kokkos_Core.hpp>

namespace BaskerNS{

  template class Basker<int, float, Kokkos::OpenMP>;
  template class Basker<long int, float, Kokkos::OpenMP>;
  template class Basker<long int, double, Kokkos::OpenMP>;
  template class Basker<int, double, Kokkos::OpenMP>;

  //template class Basker<unsigned long, double, Kokkos::OpenMP>;
  //template class Basker<long long, double, Kokkos::OpenMP>;
 
}
