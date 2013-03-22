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

#ifndef TPETRA_SERIALPLATFORM_HPP
#define TPETRA_SERIALPLATFORM_HPP

#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_Describable.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Tpetra {

  /// \brief Implementation of the Platform concept for non-MPI platforms.
  ///
  /// SerialPlatform is an implementation of Tpetra's Platform
  /// concept.  Classes implementing the Platform concept are
  /// templated on the Kokkos Node type.  They have at least the
  /// following public interface:
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
  /// SerialPlatform uses a "communicator" containing one process.  It
  /// is available whether or not Trilinos was built with MPI (the
  /// Message-Passing Interface which provides a distributed-memory
  /// parallel programming model).
  template <class Node>
  class SerialPlatform : public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{ 

    //! Kokkos Node type; the template parameter of this class.
    typedef Node NodeType;

    //@}
    //! @name Constructors and destructor
    //@{ 

    /// Constructor that accepts a Kokkos Node.
    ///
    /// \param node [in/out] The Kokkos Node instance.  If null, this
    ///   class will create a Node with default parameters.
    explicit SerialPlatform (const RCP<Node> &node) :
      comm_ (rcp (new Teuchos::SerialComm<int>())),
      node_ (node.is_null () ? Kokkos::Details::getNode<Node> () : node)
    {}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~SerialPlatform() {}

    //@}
    //! @name Methods to access the communicator and Kokkos Node.
    //@{ 

    //! The Teuchos::Comm instance with which this object was created.
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
    RCP<const Teuchos::SerialComm<int> > comm_;
    //! Kokkos Node object instantiated for the platform.
    RCP<Node> node_;

  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    SerialPlatform (const SerialPlatform<Node> &platform);
    //! Unimplemented assignment operator (syntactically forbidden).
    SerialPlatform& operator= (const SerialPlatform<Node> &platform);
  };

  /// \class SerialPlatform<Kokkos::DefaultNode::DefaultNodeType>
  /// \brief SerialPlatform specialization for Kokkos::DefaultNode::DefaultNodeType.
  ///
  /// \warning Kokkos::DefaultNode::DefaultNodeType is a typedef, and
  ///   may have a different type, depending on Trilinos' build
  ///   options.  For example, it may be Kokkos::SerialNode if
  ///   Trilinos was built without a threading library, or
  ///   Kokkos::TPINode if Trilinos was built with Pthreads.
  ///
  /// \note In the past (up to and including the 10.8 Trilinos
  ///   release), the specialization of SerialPlatform for the default
  ///   Node type delayed instantiation of the default Node instance
  ///   until getNode() was called.  We have changed this behavior to
  ///   simplify the code and make the specialization of
  ///   SerialPlatform conform more closely to the generic version of
  ///   SerialPlatform.
  template <>
  class SerialPlatform<Kokkos::DefaultNode::DefaultNodeType> : 
    public Teuchos::Describable {
  public:
    //! @name Typedefs
    //@{ 

    //! Kokkos Node type; the template parameter of this class.
    typedef Kokkos::DefaultNode::DefaultNodeType NodeType;

    //@}
    //! @name Constructors and destructor
    //@{ 

    /// \brief Default constructor: uses Kokkos default node.
    ///
    /// The specialization of SerialPlatform for the default Node type
    /// includes a default constructor.  It instantiates a default
    /// Node instance using Kokkos::DefaultNode::getDefaultNode().
    SerialPlatform() : 
      comm_ (rcp (new Teuchos::SerialComm<int> ())),
      node_ (Kokkos::DefaultNode::getDefaultNode ())
    {}

    /// \brief Constructor that accepts a Kokkos Node.
    ///
    /// This version of the constructor is declared "explicit" to
    /// forbid silent type conversions from the Node instance to a
    /// SerialPlatform.  (A single-argument constructor that is not
    /// declared "explicit" defines a type conversion method from the
    /// input type to the constructor's class's type.)  The "explicit"
    /// declaration does not affect typical use of this constructor.
    ///
    /// This specialization of SerialPlatform for the default node
    /// type will instantiate a default Node if node.is_null().
    ///
    /// \param node [in/out] The Kokkos Node instance.
    explicit SerialPlatform (const RCP<Kokkos::DefaultNode::DefaultNodeType> &node) :
      comm_ (rcp (new Teuchos::SerialComm<int> ())),
      node_ (node.is_null() ? Kokkos::DefaultNode::getDefaultNode() : node)
    {}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~SerialPlatform() {}

    //@}
    //! @name Methods to access the communicator and Kokkos Node.
    //@{ 

    //! The Teuchos::Comm instance with which this object was created.
    RCP<const Teuchos::Comm<int> > getComm() const {
      return comm_;
    }

    //! The Kokkos Node instance for this platform to use.
    RCP<Kokkos::DefaultNode::DefaultNodeType> getNode() const {
      return node_;
    }

    //@}
  private:
    //! Unimplemented copy constructor (syntactically forbidden).
    SerialPlatform (const SerialPlatform<Kokkos::DefaultNode::DefaultNodeType> &platform);

    //! Unimplemented assignment operator (syntactically forbidden).
    SerialPlatform& operator= (const SerialPlatform<Kokkos::DefaultNode::DefaultNodeType> &platform);

  protected: 
    //! Teuchos::Comm object instantiated for the platform.
    RCP<const Teuchos::SerialComm<int> > comm_;

    //! Node object instantiated for the platform.
    RCP<Kokkos::DefaultNode::DefaultNodeType> node_;
  };

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
