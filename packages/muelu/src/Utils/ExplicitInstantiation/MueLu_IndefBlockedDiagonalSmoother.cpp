// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MUELU_EXPERIMENTAL

#include "MueLu_ExplicitInstantiation.hpp"

#include "MueLu_IndefBlockedDiagonalSmoother_def.hpp"

#ifdef HAVE_MUELU_TPETRA
  #include <TpetraCore_config.h>
  #include "TpetraCore_ETIHelperMacros.h"

  #define MUELU_LOCAL_INSTANT(S,LO,GO,N) \
          template class MueLu::IndefBlockedDiagonalSmoother<S,LO,GO,N>;

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR(MUELU_LOCAL_INSTANT)
#endif

#ifdef HAVE_MUELU_EPETRA
  #ifdef HAVE_MUELU_TPETRA
    #ifndef HAVE_MUELU_TPETRA_INST_INT_INT
      #ifdef HAVE_TPETRA_INST_CUDA_DEFAULT
        // in case of Cuda being the default Node we need the default node type on the host
        typedef Kokkos::View<int>::HostMirror::execution_space default_host_execution_space;
        typedef Kokkos::Compat::KokkosDeviceWrapperNode<host_execution_space, Kokkos::HostSpace> default_node_type;
        template class MueLu::IndefBlockedDiagonalSmoother<double,int,int, default_node_type >;
      #elif HAVE_TPETRA_INST_OPENMP_DEFAULT
        template class MueLu::IndefBlockedDiagonalSmoother<double,int,int, Kokkos::Compat::KokkosOpenMPWrapperNode >;
      #elif HAVE_TPETRA_INST_SERIAL_DEFAULT
        template class MueLu::IndefBlockedDiagonalSmoother<double,int,int, Kokkos::Compat::KokkosSerialWrapperNode >;
      #elif HAVE_TPETRA_INST_PTHREAD_DEFAULT
        template class MueLu::IndefBlockedDiagonalSmoother<double,int,int, Kokkos::Compat::KokkosThreadsWrapperNode >;
      #else
        // TODO: there should be at least one default node active!! maybe we have to tweak MueLu CMakeLists.txt?
        template class MueLu::IndefBlockedDiagonalSmoother<double,int,int, Kokkos::Compat::KokkosSerialWrapperNode >;
      #endif
    #endif
  #else
    #ifndef HAVE_MUELU_TPETRA_INST_INT_INT
      // Epetra only case.  We need the host memory space and the default host execution space.
      // We have to specify the host here, because the default execution space could be CUDA.
      // If you don't want to use OpenMP, you could try to use Kokkos::Serial, as long as it is enabled.
      #ifdef KOKKOS_HAVE_SERIAL
        typedef Kokkos::Serial default_host_execution_space;
      #else
        typedef Kokkos::View<int>::HostMirror::execution_space default_host_execution_space;
      #endif // KOKKOS_HAVE_SERIAL
      typedef Kokkos::Compat::KokkosDeviceWrapperNode<host_execution_space, Kokkos::HostSpace> default_node_type;
      template class MueLu::IndefBlockedDiagonalSmoother<double, int, int, default_node_type>;
    #endif
  #endif // HAVE_MUELU_TPETRA
#endif // end ifdef HAVE_MUELU_EPETRA

#endif
