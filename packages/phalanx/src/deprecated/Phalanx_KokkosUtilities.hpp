#ifndef PHALANX_KOKKOS_UTILITIES_HPP
#define PHALANX_KOKKOS_UTILITIES_HPP

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_MDField.hpp"

namespace PHX {

  void InitializeKokkosDevice(int num_threads = -1);

  void InitializeKokkosDevice(int&  narg, char* arg[]);

  void FinalizeKokkosDevice();

  struct KokkosDeviceSession {

    KokkosDeviceSession();

    ~KokkosDeviceSession();
    
  };

}


#endif
