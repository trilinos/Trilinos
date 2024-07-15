// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
