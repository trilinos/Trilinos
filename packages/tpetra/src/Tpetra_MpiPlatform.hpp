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

#ifndef TPETRA_MPIPLATFORM_HPP
#define TPETRA_MPIPLATFORM_HPP

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Describable.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Tpetra {

namespace Details {
  //! Create a Kokkos Node instance with default parameters.
  template<class NodeType>
  Teuchos::RCP<NodeType>
  getNode() {
    Teuchos::ParameterList defaultParams;
    return Teuchos::rcp (new NodeType (defaultParams));
  }
} // namespace Details

  /// \brief Implementation of the Platform concept for MPI-based platforms.
  ///
  /// MpiPlatform is an implementation of Tpetra's Platform concept.
  /// Classes implementing the Platform concept are templated on the
  /// Kokkos Node type.  They have at least the following public
  /// interface:
  /// \code
  /// // This is not a real class; it just illustrates the concept.
  /// template<class Node>
  /// class Platform {
  /// public:
  ///   typedef Node NodeType;
  ///   explicit Platform (const RCP<Node>& node);
  ///   RCP<const Comm<int> > getComm() const;
  ///   RCP<Node> getNode() const;
  /// };
  /// \endcode
  /// MpiPlatform also has a constructor that accepts an MPI
  /// communicator, over which the application using the platform
  /// should perform communication.  The default communicator is
  /// MPI_COMM_WORLD.  MpiPlatform is only available if Trilinos was
  /// built with MPI.
  template <class Node>
  class MpiPlatform : public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{ 

    //! Kokkos Node type; the template parameter of this class.
    typedef Node NodeType;

    //@}
    //! @name Constructors and destructor
    //@{ 

    /// \brief Constructor that accepts a Kokkos Node.
    ///
    /// This version of the constructor uses MPI_COMM_WORLD as the
    /// communicator.  It is declared "explicit" to forbid silent
    /// conversions via assignment from the Node instance to an
    /// MpiPlatform.
    ///
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters.
    explicit MpiPlatform (const RCP<Node> &node) :
      comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD))),
      node_ (node.is_null () ? Details::getNode<Node> () : node)
    {}

    /// \brief Constructor that accepts a Kokkos Node and a wrapped MPI communicator.
    ///
    /// This version of the constructor accepts an arbitrary MPI
    /// communicator.  It requires that you first wrap the MPI
    /// communicator in a Teuchos::OpaqueWrapper.  This is helpful if
    /// you want to "free" the communicator automatically after use.
    ///
    /// If you just have a raw \c MPI_Comm and you want to be
    /// responsible for freeing it after use, use the constructor
    /// version that takes a raw \c MPI_Comm.  Otherwise, see the
    /// documentation of Teuchos::OpaqueWrapper to learn how to wrap
    /// an \c MPI_Comm and how to set the wrapper to free the
    /// communicator automatically after use.
    ///
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters.
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
    MpiPlatform (const RCP<Node> &node, 
                 const RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm)
      : comm_ (Teuchos::createMpiComm<int> (rawMpiComm)),
	node_ (node.is_null () ? Details::getNode<Node> () : node)
    {}

    /// \brief Constructor that accepts a Kokkos Node and a raw MPI communicator.
    ///
    /// This version of constructor accepts an arbitrary "raw" (not
    /// wrapped) MPI communicator.  You are responsible for freeing
    /// the MPI communicator after use, if necessary.
    ///
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters.
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
    MpiPlatform (const RCP<Node> &node, MPI_Comm rawMpiComm)
      : comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (rawMpiComm))),
	node_ (node.is_null () ? Details::getNode<Node> () : node)
    {}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~MpiPlatform() {}

    //@}
    //! @name Methods to access the communicator and Kokkos Node.
    //@{ 

    //! The \c Teuchos::Comm instance with which this object was created.
    RCP<const Comm<int> > getComm() const {
      return comm_; 
    }

    //! The Kokkos Node instance with which this object was created.
    RCP<Node> getNode() const {
      return node_;
    }

    //@}
  protected: 
    //! Teuchos::Comm object instantiated for the platform.
    RCP<Teuchos::MpiComm<int> > comm_;
    //! Kokkos Node object instantiated for the platform.
    RCP<Node> node_;

  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    MpiPlatform (const MpiPlatform<Node> &platform);
    //! Unimplemented assignment operator (syntactically forbidden).
    MpiPlatform& operator= (const MpiPlatform<Node> &platform);
  };


  /// \class MpiPlatform<Kokkos::DefaultNode::DefaultNodeType>
  /// \brief MpiPlatform specialization for \c Kokkos::DefaultNode::DefaultNodeType.
  ///
  /// \warning \c Kokkos::DefaultNode::DefaultNodeType is a typedef,
  ///   and may have a different type, depending on Trilinos' build
  ///   options.  For example, it may be Kokkos::SerialNode if
  ///   Trilinos was built without a threading library, or
  ///   Kokkos::TPINode if Trilinos was built with Pthreads.
  ///
  /// \note In the past (up to and including the 10.8 Trilinos
  ///   release), the specialization of MpiPlatform for the default
  ///   Node type delayed instantiation of the default Node instance
  ///   until getNode() was called.  We have changed this behavior to
  ///   simplify the code and make the specialization of MpiPlatform
  ///   conform more closely to the generic version of MpiPlatform.
  template <>
  class MpiPlatform<Kokkos::DefaultNode::DefaultNodeType> : 
    public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{ 

    //! Kokkos Node type; the template parameter of this class.
    typedef Kokkos::DefaultNode::DefaultNodeType NodeType;

    //@}
    //! @name Constructors and destructor
    //@{ 

    /// \brief Default constructor: uses Kokkos default node and MPI_COMM_WORLD.
    ///
    /// The specialization of MpiPlatform for the default node type
    /// includes a default constructor.  It instantiates a default
    /// Node instance using Kokkos::DefaultNode::getDefaultNode(), and
    /// uses MPI_COMM_WORLD as the communicator.
    MpiPlatform () :
      comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD))),
      node_ (Kokkos::DefaultNode::getDefaultNode ())
    {}

    /// \brief Constructor that accepts a Kokkos Node.
    ///
    /// This version of the constructor uses MPI_COMM_WORLD as the
    /// communicator.  It is declared "explicit" to forbid silent type
    /// conversions from the Node instance to an MpiPlatform.  (A
    /// single-argument constructor that is not declared "explicit"
    /// defines a type conversion method from the input type to the
    /// constructor's class's type.)  The "explicit" declaration does
    /// not affect typical use of this constructor.
    ///
    /// This specialization of MpiPlatform for the default node type
    /// will instantiate a default Node if node.is_null().
    ///
    /// \param node [in/out] The Kokkos Node instance.
    explicit MpiPlatform (const RCP<Kokkos::DefaultNode::DefaultNodeType> &node) :
      comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (MPI_COMM_WORLD))),
      node_ (node.is_null() ? Kokkos::DefaultNode::getDefaultNode() : node)
    {}

    /// \brief Constructor that accepts a Kokkos Node and a wrapped MPI communicator.
    ///
    /// This version of the constructor accepts an arbitrary MPI
    /// communicator.  It requires that you first wrap the MPI
    /// communicator in a Teuchos::OpaqueWrapper.  This is helpful if
    /// you want to "free" the communicator automatically after use.
    ///
    /// If you just have a raw \c MPI_Comm and you want to be
    /// responsible for freeing it after use, use the constructor
    /// version that takes a raw \c MPI_Comm.  Otherwise, see the
    /// documentation of Teuchos::OpaqueWrapper to learn how to wrap
    /// an \c MPI_Comm and how to set the wrapper to free the
    /// communicator automatically after use.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The MPI communicator, wrapped in a
    ///   <tt>Teuchos::OpaqueWrapper</tt>.
    MpiPlatform (const RCP<Kokkos::DefaultNode::DefaultNodeType> &node, 
                 const RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm) 
      : comm_ (Teuchos::createMpiComm<int> (rawMpiComm)),
        node_ (node.is_null() ? Kokkos::DefaultNode::getDefaultNode() : node)
    {}

    /// \brief Constructor that accepts a Kokkos Node and a raw MPI communicator.
    ///
    /// This version of constructor accepts an arbitrary "raw" (not
    /// wrapped) MPI communicator.  You are responsible for freeing
    /// the MPI communicator after use, if necessary.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    /// \param rawMpiComm [in] The "raw" (not wrapped) MPI
    ///   communicator.
    MpiPlatform (const RCP<Kokkos::DefaultNode::DefaultNodeType> &node, 
		 MPI_Comm rawMpiComm)
      : comm_ (Teuchos::createMpiComm<int> (Teuchos::opaqueWrapper<MPI_Comm> (rawMpiComm))),
	node_ (node.is_null() ? Kokkos::DefaultNode::getDefaultNode() : node)
    {}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~MpiPlatform() {}

    //@}
    //! @name Methods to access the communicator and Kokkos Node.
    //@{ 

    //! The \c Teuchos::Comm instance with which this object was created.
    RCP<const Comm<int> > getComm() const {
      return comm_; 
    }

    /// \brief The Kokkos Node instance for this platform to use.
    /// 
    /// This MpiPlatform specialization's constructor may have created
    /// the Node instance, if you invoked the default constructor or
    /// passed in a null Node pointer.
    RCP<Kokkos::DefaultNode::DefaultNodeType> getNode() const {
      return node_;
    }

    //@}

  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    MpiPlatform (const MpiPlatform<Kokkos::DefaultNode::DefaultNodeType> &platform);

    //! Unimplemented assignment operator (syntactically forbidden).
    MpiPlatform& operator= (const MpiPlatform<Kokkos::DefaultNode::DefaultNodeType> &platform);

  protected: 
    //! Teuchos::Comm object instantiated for the platform.
    RCP<Teuchos::MpiComm<int> > comm_;

    //! Kokkos Node object instantiated for the platform.
    RCP<Kokkos::DefaultNode::DefaultNodeType> node_;
  };


} // namespace Tpetra

#endif // TPETRA_MPIPLATFORM_HPP
