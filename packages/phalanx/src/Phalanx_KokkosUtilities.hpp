#ifndef PHALANX_KOKKOS_UTILITIES_HPP
#define PHALANX_KOKKOS_UTILITIES_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDField.hpp"

namespace PHX {
  
  void InitializeKokkosDevice(const int& num_threads = 1);

  void FinalizeKokkosDevice();

  struct KokkosDeviceSession {

    KokkosDeviceSession();

    ~KokkosDeviceSession();
    
  };

}


#endif
