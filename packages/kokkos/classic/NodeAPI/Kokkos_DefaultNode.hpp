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

#ifndef KOKKOS_DEFAULT_NODE_HPP_
#define KOKKOS_DEFAULT_NODE_HPP_

#include "Kokkos_ConfigDefs.hpp"
#include "KokkosClassic_DefaultNode_config.h"
#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#include "Kokkos_OpenMPNode.hpp"
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
#include "Kokkos_ThrustGPUNode.hpp"
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace Kokkos {

namespace Details {
  /// \fn getNode
  /// \brief Create a Kokkos Node instance with default parameters.
  /// \tparam NodeType The Kokkos Node type.
  ///
  /// \warning This function is <i>not</i> safe to be called by
  ///   multiple threads simultaneously.  The first call to this
  ///   function must be serialized.  Also, RCP is not currently
  ///   thread safe.
  ///
  /// Every Kokkos Node's constructor takes a Teuchos::ParameterList.
  /// We presume that for every Kokkos Node, if that list of
  /// parameters is empty, then the Node will use default parameters.
  /// This is true for all the Node types implemented in Kokkos.
  template<class NodeType>
  Teuchos::RCP<NodeType> getNode () {
    static Teuchos::RCP<NodeType> theNode;
    if (theNode.is_null ()) {
      Teuchos::ParameterList defaultParams;
      theNode = Teuchos::rcp (new NodeType (defaultParams));
    }
    return theNode;
  }
} // namespace Details

/// \class DefaultNode
/// \brief Class that specifies Kokkos' default Node type, and
///   creates and holds an instance of it for you to use.
/// \ingroup kokkos_node_api
class DefaultNode {
public:
  //! The default Node type.
  typedef SerialNode DefaultNodeType;

  //! Return a pointer to the default Node.
  static RCP<DefaultNodeType> getDefaultNode();
};

} // namespace Kokkos

#endif
