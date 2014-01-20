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

#include "Kokkos_DefaultNode.hpp"

namespace KokkosClassic {

  Teuchos::RCP<DefaultNode::DefaultNodeType> DefaultNode::getDefaultNode()
  {
    return Details::getNode<DefaultNodeType> ();
  }

#ifdef HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT
  namespace Details {

  // Full specializations for the new Kokkos wrapper Node types.

#ifdef KOKKOS_HAVE_CUDA
  template<>
  Teuchos::RCP< ::Kokkos::Compat::KokkosCudaWrapperNode>
  getNode< ::Kokkos::Compat::KokkosCudaWrapperNode> () {
    Teuchos::ParameterList defaultParams;
    // This Node type knows how to handle multiple Node instances,
    // so we don't have to keep a static instance around.
    return Teuchos::rcp (new ::Kokkos::Compat::KokkosCudaWrapperNode (defaultParams));
  }
#endif // KOKKOS_HAVE_CUDA

#ifdef KOKKOS_HAVE_OPENMP
  template<>
  Teuchos::RCP< ::Kokkos::Compat::KokkosOpenMPWrapperNode>
  getNode< ::Kokkos::Compat::KokkosOpenMPWrapperNode> () {
    Teuchos::ParameterList defaultParams;
    // This Node type knows how to handle multiple Node instances,
    // so we don't have to keep a static instance around.
    return Teuchos::rcp (new ::Kokkos::Compat::KokkosOpenMPWrapperNode (defaultParams));
  }
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
  template<>
  Teuchos::RCP< ::Kokkos::Compat::KokkosThreadsWrapperNode>
  getNode< ::Kokkos::Compat::KokkosThreadsWrapperNode> () {
    Teuchos::ParameterList defaultParams;
    // This Node type knows how to handle multiple Node instances,
    // so we don't have to keep a static instance around.
    return Teuchos::rcp (new ::Kokkos::Compat::KokkosThreadsWrapperNode (defaultParams));
  }
#endif // KOKKOS_HAVE_PTHREAD

  } // namespace Details
#endif // HAVE_KOKKOSCLASSIC_KOKKOSCOMPAT

} // namespace KokkosClassic
