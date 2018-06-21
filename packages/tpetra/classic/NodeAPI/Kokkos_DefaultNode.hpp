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

#ifndef KOKKOS_DEFAULTNODE_HPP
#define KOKKOS_DEFAULTNODE_HPP

#include "Kokkos_ConfigDefs.hpp"
#include "KokkosClassic_DefaultNode_config.h"
#include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Teuchos_RCP.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
namespace Teuchos {
  // Dear users: This is just a forward declaration.
  // Please skip over it.
  class ParameterList;
} // namespace Teuchos
#endif // DOXYGEN_SHOULD_SKIP_THIS

namespace KokkosClassic {

namespace Details {
  /// \fn getNode
  /// \brief Create and return a Kokkos Node instance.
  /// \tparam NodeType The Kokkos Node type.
  ///
  /// \warning This function is <i>not</i> safe to be called by
  ///   multiple threads simultaneously.  The first call to this
  ///   function must be serialized.  Also, RCP is not currently
  ///   thread safe.
  ///
  /// \param params [in/out] On input: Any parameters that the Kokkos
  ///   Node accepts.  On output, the list may be modified to include
  ///   missing parameters and their default values.  If params is
  ///   null, default parameters will be used.
  ///
  /// Every Kokkos Node's constructor takes a Teuchos::ParameterList.
  /// We presume that for every Kokkos Node, if that list of
  /// parameters is empty, then the Node will use default parameters.
  /// This is true for all the Node types implemented in Kokkos.
  template<class NodeType>
  Teuchos::RCP<NodeType>
  getNode (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

#ifdef KOKKOS_ENABLE_CUDA
  extern template Teuchos::RCP< ::Kokkos::Compat::KokkosCudaWrapperNode> getNode< ::Kokkos::Compat::KokkosCudaWrapperNode> (const Teuchos::RCP<Teuchos::ParameterList>& );
#endif // KOKKOS_ENABLE_CUDA

#ifdef KOKKOS_ENABLE_OPENMP
  extern template Teuchos::RCP< ::Kokkos::Compat::KokkosOpenMPWrapperNode> getNode< ::Kokkos::Compat::KokkosOpenMPWrapperNode> (const Teuchos::RCP<Teuchos::ParameterList>& );
#endif // KOKKOS_ENABLE_OPENMP

#ifdef KOKKOS_ENABLE_SERIAL
  extern template Teuchos::RCP< ::Kokkos::Compat::KokkosSerialWrapperNode> getNode< ::Kokkos::Compat::KokkosSerialWrapperNode> (const Teuchos::RCP<Teuchos::ParameterList>& );
#endif // KOKKOS_ENABLE_SERIAL

#if defined(KOKKOS_ENABLE_PTHREAD) || defined(KOKKOS_ENABLE_THREADS)
  extern template Teuchos::RCP< ::Kokkos::Compat::KokkosThreadsWrapperNode> getNode< ::Kokkos::Compat::KokkosThreadsWrapperNode> (const Teuchos::RCP<Teuchos::ParameterList>& );
#endif // KOKKOS_ENABLE_PTHREAD || KOKKOS_ENABLE_THREADS

} // namespace Details

  /// \brief Specify Tpetra's default Node type.
  ///
  /// Tpetra::Map uses this class to get Tpetra's default Node type.
  /// <i>This is an implementation detail of Tpetra</i>.  If you want
  /// to know the default Node type, just ask Tpetra::Map, like this:
  /// \code
  /// typedef Tpetra::Map<>::node_type default_node_type;
  /// \endcode
  class DefaultNode {
  public:
#if defined(HAVE_TPETRA_DEFAULTNODE_CUDAWRAPPERNODE)
    typedef ::Kokkos::Compat::KokkosCudaWrapperNode DefaultNodeType;
#elif defined(HAVE_TPETRA_DEFAULTNODE_OPENMPWRAPPERNODE)
    typedef ::Kokkos::Compat::KokkosOpenMPWrapperNode DefaultNodeType;
#elif defined(HAVE_TPETRA_DEFAULTNODE_THREADSWRAPPERNODE)
    typedef ::Kokkos::Compat::KokkosThreadsWrapperNode DefaultNodeType;
#elif defined(HAVE_TPETRA_DEFAULTNODE_SERIALWRAPPERNODE)
    typedef ::Kokkos::Compat::KokkosSerialWrapperNode DefaultNodeType;
#else
#    error "No default Kokkos Node type specified.  Please set the CMake option Tpetra_DefaultNode to a valid Node type."
#endif

    //! \brief Return a pointer to the default node.
    static Teuchos::RCP<DefaultNodeType> getDefaultNode();
  };

} // namespace KokkosClassic

#endif // KOKKOS_DEFAULTNODE_HPP

