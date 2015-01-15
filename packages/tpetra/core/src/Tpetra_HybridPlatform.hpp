// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_HYBRIDPLATFORM_HPP
#define TPETRA_HYBRIDPLATFORM_HPP

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Kokkos_DefaultNode.hpp>

// This macro is only for use by Tpetra developers.
// It should only be invoked in the Tpetra namespace,
// outside of the HybridPlatform class declaration.
#define TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(N) \
  template <> bool HybridPlatform::isNodeSupported<N> ();

namespace Tpetra {

/// \brief A platform class for hybrid nodes.
/// \warning This class is DEPRECATED and will be removed in the next
///   major release (12.0) of Trilinos.
class TPETRA_DEPRECATED HybridPlatform : public Teuchos::Describable {
private:
  typedef KokkosClassic::DefaultNode::DefaultNodeType default_node_type;

public:
  //! @name Constructor/Destructor Methods
  //@{

  //! Constructor
  HybridPlatform (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                  Teuchos::ParameterList& pl);

  //! Destructor
  ~HybridPlatform ();

  //@}
  //! @name Class Query, Creation and Accessor Methods
  //@{

  //! Comm Instance
  Teuchos::RCP<const Teuchos::Comm<int> > getComm () const;

  //! List of supported nodes and their valid parameters.
  static Teuchos::RCP<Teuchos::ParameterList> listSupportedNodes ();

  //! Whether HybridPlatform supports the given \c Node type.
  template <class Node>
  static bool isNodeSupported ();

  /// \brief Run user code with the runtime-selected Node type.
  ///
  /// This method assumes that UserCode is a class with a template
  /// parameter Node, which has a class ("static") method run():
  /// \code
  /// template<class Node>
  /// class UserCode {
  /// public:
  ///   static void
  ///   run (Teuchos::ParameterList& plist,
  ///        Teuchos::RCP<const Teuchos::Comm<int> > comm,
  ///        Teuchos::RCP<Node> node);
  /// };
  /// \endcode
  /// Note that this method depends on the "template parameter that
  /// takes a template parameter" feature of C++11.  Your compiler
  /// may or may not support this feature.  If it does, you may have
  /// to use a special compiler flag to enable the feature.
  template <template <class Node> class UserCode>
  void runUserCode ();

  /// \brief Run user code with the runtime-selected Node type.
  ///
  /// This method, unlike the version of runUserCode that takes no
  /// arguments above, assumes that UserCode is a class with an
  /// <i>instance</i> (not class) method run():
  /// \code
  /// class UserCode {
  /// public:
  ///   template<class Node>
  ///   void
  ///   run (Teuchos::ParameterList& plist,
  ///        Teuchos::RCP<const Teuchos::Comm<int> > comm,
  ///        Teuchos::RCP<Node> node);
  /// };
  /// \endcode
  template <class UserCode>
  void runUserCode (UserCode &code);

  //@}

protected:
  void createNode ();

private:
  HybridPlatform(const HybridPlatform &platform); // not supported
  const Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  Teuchos::ParameterList instList_;
  bool nodeCreated_;

  // DEFAULTNODE is always the default Node type, whatever that
  // happens to be.  We include this just in case SERIALNODE is not
  // enabled, so that the syntax won't break (e.g., so that the enum
  // won't be empty).
  enum NodeType {
    DEFAULTNODE
#ifdef HAVE_KOKKOSCLASSIC_SERIAL
    , SERIALNODE
#endif // HAVE_KOKKOSCLASSIC_SERIAL
#ifdef HAVE_KOKKOSCLASSIC_TBB
    , TBBNODE
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
    , TPINODE
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
    , OMPNODE
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
    , THRUSTGPUNODE
#endif
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#ifdef KOKKOS_HAVE_SERIAL
    , SERIAL_WRAPPER_NODE
#endif // KOKKOS_HAVE_SERIAL
#ifdef KOKKOS_HAVE_OPENMP
    , OPENMP_WRAPPER_NODE
#endif // KOKKOS_HAVE_OPENMP
#ifdef KOKKOS_HAVE_PTHREAD
    , THREADS_WRAPPER_NODE
#endif // KOKKOS_HAVE_PTHREAD
#ifdef KOKKOS_HAVE_CUDA
    , CUDA_WRAPPER_NODE
#endif // KOKKOS_HAVE_CUDA
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
  } nodeType_;

  //! Instance of the default Node type.
  Teuchos::RCP<default_node_type> defaultNode_;

#ifdef HAVE_KOKKOSCLASSIC_SERIAL
  Teuchos::RCP<KokkosClassic::SerialNode> serialNode_;
#endif // HAVE_KOKKOSCLASSIC_SERIAL
#ifdef HAVE_KOKKOSCLASSIC_TBB
  Teuchos::RCP<KokkosClassic::TBBNode> tbbNode_;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
  Teuchos::RCP<KokkosClassic::TPINode> tpiNode_;
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
  Teuchos::RCP<KokkosClassic::OpenMPNode> ompNode_;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
  Teuchos::RCP<KokkosClassic::ThrustGPUNode> thrustNode_;
#endif
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#ifdef KOKKOS_HAVE_SERIAL
  Teuchos::RCP<Kokkos::Compat::KokkosSerialWrapperNode> serialWrapperNode_;
#endif // KOKKOS_HAVE_SERIAL
#ifdef KOKKOS_HAVE_OPENMP
  Teuchos::RCP<Kokkos::Compat::KokkosOpenMPWrapperNode> ompWrapperNode_;
#endif // KOKKOS_HAVE_OPENMP
#ifdef KOKKOS_HAVE_PTHREAD
  Teuchos::RCP<Kokkos::Compat::KokkosThreadsWrapperNode> threadsWrapperNode_;
#endif // KOKKOS_HAVE_PTHREAD
#ifdef KOKKOS_HAVE_CUDA
  Teuchos::RCP<Kokkos::Compat::KokkosCudaWrapperNode> cudaWrapperNode_;
#endif // KOKKOS_HAVE_CUDA
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
};


template <class Node>
bool HybridPlatform::isNodeSupported ()
{
  return false;
}


template <class UserCode>
void
HybridPlatform::runUserCode (UserCode& codeobj)
{
  createNode ();
  switch (nodeType_) {
  case DEFAULTNODE:
    codeobj.template run<default_node_type> (instList_, comm_, defaultNode_);
    break;
#ifdef HAVE_KOKKOSCLASSIC_SERIAL
  case SERIALNODE:
    codeobj.template run<KokkosClassic::SerialNode> (instList_, comm_, serialNode_);
    break;
#endif // HAVE_KOKKOSCLASSIC_SERIAL
#ifdef HAVE_KOKKOSCLASSIC_TBB
  case TBBNODE:
    codeobj.template run<KokkosClassic::TBBNode> (instList_, comm_, tbbNode_);
    break;
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
  case OMPNODE:
    codeobj.template run<KokkosClassic::OpenMPNode> (instList_, comm_, ompNode_);
    break;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
  case TPINODE:
    codeobj.template run<KokkosClassic::TPINode> (instList_, comm_, tpiNode_);
    break;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
  case THRUSTGPUNODE:
    codeobj.template run<KokkosClassic::ThrustGPUNode> (instList_, comm_, thrustNode_);
    break;
#endif
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#ifdef KOKKOS_HAVE_SERIAL
  case SERIAL_WRAPPER_NODE:
    codeobj.template run<Kokkos::Compat::KokkosSerialWrapperNode> (instList_, comm_,
                                                                   serialWrapperNode_);
    break;
#endif // KOKKOS_HAVE_SERIAL
#ifdef KOKKOS_HAVE_OPENMP
  case OPENMP_WRAPPER_NODE:
    codeobj.template run<Kokkos::Compat::KokkosOpenMPWrapperNode> (instList_, comm_,
                                                                   ompWrapperNode_);
    break;
#endif // KOKKOS_HAVE_OPENMP
#ifdef KOKKOS_HAVE_PTHREAD
  case THREADS_WRAPPER_NODE:
    codeobj.template run<Kokkos::Compat::KokkosThreadsWrapperNode> (instList_, comm_,
                                                                    threadsWrapperNode_);
    break;
#endif // KOKKOS_HAVE_PTHREAD
#ifdef KOKKOS_HAVE_CUDA
  case CUDA_WRAPPER_NODE:
    codeobj.template run<Kokkos::Compat::KokkosCudaWrapperNode> (instList_, comm_,
                                                                 cudaWrapperNode_);
    break;
#endif // KOKKOS_HAVE_CUDA
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, Teuchos::typeName(*this) << "::runUserCode: "
      "Invalid node type." << std::endl);
  } // end of switch
}

template <template<class Node> class UserCode>
void
HybridPlatform::runUserCode ()
{
  createNode ();
  switch (nodeType_) {
  case DEFAULTNODE:
    UserCode<default_node_type>::run (instList_, comm_, defaultNode_);
    break;
#ifdef HAVE_KOKKOSCLASSIC_SERIAL
  case SERIALNODE:
    UserCode<KokkosClassic::SerialNode>::run (instList_, comm_, serialNode_);
    break;
#endif // HAVE_KOKKOSCLASSIC_SERIAL
#ifdef HAVE_KOKKOSCLASSIC_TBB
  case TBBNODE:
    UserCode<KokkosClassic::TBBNode>::run (instList_, comm_, tbbNode_);
    break;
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
  case OMPNODE:
    UserCode<KokkosClassic::OpenMPNode>::run (instList_, comm_, ompNode_);
    break;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
  case TPINODE:
    UserCode<KokkosClassic::TPINode>::run (instList_, comm_, tpiNode_);
    break;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
  case THRUSTGPUNODE:
    UserCode<KokkosClassic::ThrustGPUNode>::run (instList_, comm_, thrustNode_);
    break;
#endif
#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#ifdef KOKKOS_HAVE_SERIAL
  case SERIAL_WRAPPER_NODE:
    UserCode<Kokkos::Compat::KokkosSerialWrapperNode>::run (instList_, comm_,
                                                            serialWrapperNode_);
    break;
#endif // KOKKOS_HAVE_SERIAL
#ifdef KOKKOS_HAVE_OPENMP
  case OPENMP_WRAPPER_NODE:
    UserCode<Kokkos::Compat::KokkosOpenMPWrapperNode>::run (instList_, comm_,
                                                            ompWrapperNode_);
    break;
#endif // KOKKOS_HAVE_OPENMP
#ifdef KOKKOS_HAVE_PTHREAD
  case THREADS_WRAPPER_NODE:
    UserCode<Kokkos::Compat::KokkosThreadsWrapperNode>::run (instList_, comm_,
                                                             threadsWrapperNode_);
    break;
#endif // KOKKOS_HAVE_PTHREAD
#ifdef KOKKOS_HAVE_CUDA
  case CUDA_WRAPPER_NODE:
    UserCode<Kokkos::Compat::KokkosCudaWrapperNode>::run (instList_, comm_,
                                                          cudaWrapperNode_);
    break;
#endif // KOKKOS_HAVE_CUDA
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, Teuchos::typeName(*this) << "::runUserCode: "
      "Invalid node type." << std::endl);
  } // end of switch
}

#ifdef HAVE_KOKKOSCLASSIC_SERIAL
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(KokkosClassic::SerialNode)
#endif // HAVE_KOKKOSCLASSIC_SERIAL

#ifdef HAVE_KOKKOSCLASSIC_TBB
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(KokkosClassic::TBBNode)
#endif

#ifdef HAVE_KOKKOSCLASSIC_OPENMP
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(KokkosClassic::OpenMPNode)
#endif

#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(KokkosClassic::TPINode)
#endif

#ifdef HAVE_KOKKOSCLASSIC_THRUST
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(KokkosClassic::ThrustGPUNode)
#endif

#ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#ifdef KOKKOS_HAVE_SERIAL
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::Compat::KokkosSerialWrapperNode)
#endif // KOKKOS_HAVE_SERIAL

#ifdef KOKKOS_HAVE_OPENMP
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::Compat::KokkosOpenMPWrapperNode)
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::Compat::KokkosThreadsWrapperNode)
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::Compat::KokkosCudaWrapperNode)
#endif // KOKKOS_HAVE_CUDA
#endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT

// Make sure that the default Node type is always supported.  We can
// only do this if it's not any of the Node types listed above,
// because the default Node type is just a typedef, and we don't want
// to instantiate support for it twice.
#if ! defined(HAVE_KOKKOSCLASSIC_SERIAL) && ! defined(HAVE_KOKKOSCLASSIC_TBB) && ! defined(HAVE_KOKKOSCLASSIC_OPENMP) && ! defined(HAVE_KOKKOSCLASSIC_THREADPOOL) && ! defined(HAVE_KOKKOSCLASSIC_THRUST)
#  ifdef HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#    if ! defined (KOKKOS_HAVE_SERIAL) && ! defined(KOKKOS_HAVE_OPENMP) && ! defined(KOKKOS_HAVE_PTHREAD) && ! defined(KOKKOS_HAVE_CUDA)
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(KokkosClassic::DefaultNode::DefaultNodeType)
#    endif
#  else
TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(KokkosClassic::DefaultNode::DefaultNodeType)
#  endif // HAVE_TPETRACORE_TEUCHOSKOKKOSCOMPAT
#endif

} // namespace Tpetra

#endif // TPETRA_HYBRIDPLATFORM_HPP
