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

#include <Kokkos_SerialNode.hpp>
#ifdef HAVE_KOKKOSCLASSIC_TBB
#include <Kokkos_TBBNode.hpp>
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
#include <Kokkos_TPINode.hpp>
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
#include <Kokkos_OpenMPNode.hpp>
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
#include <Kokkos_ThrustGPUNode.hpp>
#endif

// This macro is only for use by Tpetra developers.
// It should only be invoked in the Tpetra namespace,
// outside of the HybridPlatform class declaration.
#define TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(N) \
  template <> bool HybridPlatform::isNodeSupported<N> ();

namespace Tpetra {

  //! A platform class for hybrid nodes.
  class HybridPlatform : public Teuchos::Describable {
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
    Teuchos::RCP<Kokkos::SerialNode>    serialNode_;
    bool nodeCreated_;
#ifdef HAVE_KOKKOSCLASSIC_TBB
    Teuchos::RCP<Kokkos::TBBNode>       tbbNode_;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
    Teuchos::RCP<Kokkos::TPINode>       tpiNode_;
#endif
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
    Teuchos::RCP<Kokkos::OpenMPNode>    ompNode_;
#endif
#ifdef HAVE_KOKKOSCLASSIC_THRUST
    Teuchos::RCP<Kokkos::ThrustGPUNode> thrustNode_;
#endif

    enum NodeType {
      SERIALNODE
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
    } nodeType_;
  };

  template <class Node>
  bool HybridPlatform::isNodeSupported ()
  {
    return false;
  }
  
  template <class UserCode>
  void HybridPlatform::runUserCode (UserCode& codeobj) {
    createNode();
    switch (nodeType_) {
      case SERIALNODE:
        codeobj.template run<Kokkos::SerialNode>(instList_,comm_, serialNode_);
        break;
#ifdef HAVE_KOKKOSCLASSIC_TBB
      case TBBNODE:
        codeobj.template run<Kokkos::TBBNode>(instList_,comm_, tbbNode_);
        break;
#endif        
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
      case OMPNODE:
        codeobj.template run<Kokkos::OpenMPNode>(instList_,comm_, ompNode_);
        break;
#endif        
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
      case TPINODE:
        codeobj.template run<Kokkos::TPINode>(instList_,comm_, tpiNode_);
        break;
#endif        
#ifdef HAVE_KOKKOSCLASSIC_THRUST
      case THRUSTGPUNODE:
        codeobj.template run<Kokkos::ThrustGPUNode>(instList_,comm_, thrustNode_);
        break;
#endif        
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
            Teuchos::typeName(*this) << "::runUserCode(): Invalid node type." << std::endl);
    } // end of switch
  }

  template <template<class Node> class UserCode>
  void HybridPlatform::runUserCode() {
    createNode();
    switch (nodeType_) {
      case SERIALNODE:
        UserCode<Kokkos::SerialNode>::run(instList_,comm_, serialNode_);
        break;
#ifdef HAVE_KOKKOSCLASSIC_TBB
      case TBBNODE:
        UserCode<Kokkos::TBBNode>::run(instList_,comm_, tbbNode_);
        break;
#endif        
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
      case OMPNODE:
        UserCode<Kokkos::OpenMPNode>::run(instList_,comm_, ompNode_);
        break;
#endif        
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
      case TPINODE:
        UserCode<Kokkos::TPINode>::run(instList_,comm_, tpiNode_);
        break;
#endif        
#ifdef HAVE_KOKKOSCLASSIC_THRUST
      case THRUSTGPUNODE:
        UserCode<Kokkos::ThrustGPUNode>::run(instList_,comm_, thrustNode_);
        break;
#endif        
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, 
            Teuchos::typeName(*this) << "::runUserCode(): Invalid node type." << std::endl);
    } // end of switch
  }

  TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::SerialNode)
#ifdef HAVE_KOKKOSCLASSIC_TBB
  TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::TBBNode)
#endif        
#ifdef HAVE_KOKKOSCLASSIC_OPENMP
  TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::OpenMPNode)
#endif        
#ifdef HAVE_KOKKOSCLASSIC_THREADPOOL
  TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::TPINode)
#endif        
#ifdef HAVE_KOKKOSCLASSIC_THRUST
  TPETRA_HYBRIDPLATFORM_ADD_NODE_SUPPORT_DECL(Kokkos::ThrustGPUNode)
#endif

} // namespace Tpetra

#endif // TPETRA_HYBRIDPLATFORM_HPP
