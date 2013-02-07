//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef KOKKOS_CONFIGDEFS_H
#define KOKKOS_CONFIGDEFS_H

#include <KokkosClassic_config.h>
// mem management
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

// TODO (mfh 07 Feb 2013) We want to move this into CMake.
#if defined(__GNUC__)
#  define KOKKOSCLASSIC_RESTRICT __restrict
#elif defined(__INTEL_COMPILER)
#  define KOKKOSCLASSIC_RESTRICT restrict
#else
#  define KOKKOSCLASSIC_RESTRICT 
#endif // __GNUC__

//! Namespace for Kokkos classes and methods
namespace Kokkos {
  // import Teuchos memory management classes into Kokkos
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::arcp;
  using Teuchos::arcpFromArrayView;
  using Teuchos::rcpFromRef;
  using Teuchos::null;
  using Teuchos::arcp_reinterpret_cast;
  using Teuchos::Array;
  using Teuchos::ParameterList;
}

//! Namespace for Kokkos example classes
namespace KokkosExamples {
}


namespace Kokkos {

//! Sweep direction for Gauss-Seidel or Symmetric Over-Relaxation (SOR).
enum ESweepDirection {
  Forward = 0,
  Backward
};

} // namespace Kokkos

#endif /* KOKKOS_CONFIGDEFS_H */
