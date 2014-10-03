// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

// Tests
#include "TestAssembly.hpp"

// Devices
#include "Kokkos_Core.hpp"

template <typename Storage>
void mainHost(const Teuchos::RCP<const Teuchos::Comm<int> >& comm ,
              const int use_print ,
              const int use_trials ,
              const int use_atomic ,
              const int use_nodes[] ,
              const bool check ,
              Kokkos::DeviceConfig dev_config) {
#ifdef __MIC__
  const int entry_min = 8;
  const int entry_max = 48;
  const int entry_step = 8;
#else
  const int entry_min = 4;
  const int entry_max = 32;
  const int entry_step = 4;
  // const int entry_min = 16;
  // const int entry_max = 16;
  // const int entry_step = 16;
#endif

  performance_test_driver<Storage,entry_min,entry_max,entry_step>(
    comm, use_print, use_trials, use_atomic, use_nodes, check, dev_config);
}

#ifdef KOKKOS_HAVE_PTHREAD
template void mainHost< Stokhos::StaticFixedStorage<int,double,1,Kokkos::Threads> >(const Teuchos::RCP<const Teuchos::Comm<int> >& comm , const int use_print , const int use_trials , const int use_atomic , const int use_nodes[] , const bool check , Kokkos::DeviceConfig dev_config);
#endif

#ifdef KOKKOS_HAVE_OPENMP
template void mainHost< Stokhos::StaticFixedStorage<int,double,1,Kokkos::OpenMP> >(const Teuchos::RCP<const Teuchos::Comm<int> >& comm , const int use_print , const int use_trials , const int use_atomic , const int use_nodes[] , const bool check , Kokkos::DeviceConfig dev_config);
#endif
