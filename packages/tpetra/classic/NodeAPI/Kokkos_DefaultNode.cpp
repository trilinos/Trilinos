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

namespace Details {

template<class NodeType>
Teuchos::RCP<NodeType>
getNode (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  static Teuchos::RCP<NodeType> theNode;
  if (theNode.is_null ()) {
    if (params.is_null ()) {
      Teuchos::ParameterList defaultParams;
      theNode = Teuchos::rcp (new NodeType (defaultParams));
    } else {
      theNode = Teuchos::rcp (new NodeType (*params));
    }
  }
  return theNode;
}

#ifdef KOKKOS_ENABLE_CUDA
template Teuchos::RCP< ::Kokkos::Compat::KokkosCudaWrapperNode> getNode< ::Kokkos::Compat::KokkosCudaWrapperNode> (const Teuchos::RCP<Teuchos::ParameterList>& );
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_OPENMP
template Teuchos::RCP< ::Kokkos::Compat::KokkosOpenMPWrapperNode> getNode< ::Kokkos::Compat::KokkosOpenMPWrapperNode> (const Teuchos::RCP<Teuchos::ParameterList>& );
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_SERIAL
template Teuchos::RCP< ::Kokkos::Compat::KokkosSerialWrapperNode> getNode< ::Kokkos::Compat::KokkosSerialWrapperNode> (const Teuchos::RCP<Teuchos::ParameterList>& );
#endif // KOKKOS_ENABLE_SERIAL

#if defined(KOKKOS_ENABLE_PTHREAD) || defined(KOKKOS_ENABLE_THREADS)
template Teuchos::RCP< ::Kokkos::Compat::KokkosThreadsWrapperNode> getNode< ::Kokkos::Compat::KokkosThreadsWrapperNode> (const Teuchos::RCP<Teuchos::ParameterList>& );
#endif // KOKKOS_ENABLE_PTHREAD || KOKKOS_ENABLE_THREADS

} // namespace Details

Teuchos::RCP<DefaultNode::DefaultNodeType>
DefaultNode::getDefaultNode()
{
  return Details::getNode<DefaultNodeType> ();
}

} // namespace KokkosClassic
