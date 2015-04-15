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
#include "Kokkos_NodeAPIConfigDefs.hpp"
#include "KokkosClassic_DefaultNode_config.h"
#include "Kokkos_BufferMacros.hpp"

#ifdef HAVE_TPETRACLASSIC_SERIAL
#  include "Kokkos_SerialNode.hpp"
#endif // HAVE_TPETRACLASSIC_SERIAL
#ifdef HAVE_TPETRACLASSIC_TBB
#  include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_TPETRACLASSIC_THREADPOOL
#  include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_TPETRACLASSIC_OPENMP
#  include "Kokkos_OpenMPNode.hpp"
#endif
#ifdef HAVE_TPETRACLASSIC_TEUCHOSKOKKOSCOMPAT
#  include "KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

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
  getNode (const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
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

} // namespace Details

  /** \brief Class to specify %Kokkos default node type and instantiate the default node.
      \ingroup kokkos_node_api
    */
  class DefaultNode {
  public:
#if defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TPINODE)
    typedef TPINode DefaultNodeType;
#elif defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_TBBNODE)
    typedef TBBNode DefaultNodeType;
#elif defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPNODE)
    typedef OpenMPNode DefaultNodeType;
#elif defined(HAVE_TPETRACLASSIC_TEUCHOSKOKKOSCOMPAT)
#  if defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_CUDAWRAPPERNODE)
    typedef ::Kokkos::Compat::KokkosCudaWrapperNode DefaultNodeType;
#  elif defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_OPENMPWRAPPERNODE)
    typedef ::Kokkos::Compat::KokkosOpenMPWrapperNode DefaultNodeType;
#  elif defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_THREADSWRAPPERNODE)
    typedef ::Kokkos::Compat::KokkosThreadsWrapperNode DefaultNodeType;
#  elif defined(HAVE_KOKKOSCLASSIC_DEFAULTNODE_SERIALWRAPPERNODE)
    typedef ::Kokkos::Compat::KokkosSerialWrapperNode DefaultNodeType;
#  elif defined(HAVE_TPETRACLASSIC_SERIAL)
    typedef ::KokkosClassic::DoNotUse::SerialNode DefaultNodeType;
#  else
#    error "No default Kokkos Node type specified.  Please set the CMake option KokkosClassic_DefaultNode to a valid Node type."
#  endif // defined(HAVE_TPETRACLASSIC_TEUCHOSKOKKOSCOMPAT)
#elif defined(HAVE_TPETRACLASSIC_SERIAL)
    //! Typedef specifying the default node type.
    typedef ::KokkosClassic::DoNotUse::SerialNode DefaultNodeType;
#else
#  error "No default Kokkos Node type specified.  Please set the CMake option KokkosClassic_DefaultNode to a valid Node type."
#endif

      //! \brief Return a pointer to the default node.
      static RCP<DefaultNodeType> getDefaultNode();
  };

} // namespace KokkosClassic

#endif // KOKKOS_DEFAULTNODE_HPP

