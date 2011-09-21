// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_SERIALPLATFORM_HPP
#define TPETRA_SERIALPLATFORM_HPP

#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Tpetra {

  /// \brief Implementation of the Platform concept for MPI-based platforms.
  ///
  /// SerialPlatform is an implementation of Tpetra's Platform
  /// concept.  Classes implementing Tpetra's Platform concept are
  /// templated on the Kokkos Node type.  They have at least the
  /// following public interface:
  /// \code
  /// template<class Node>
  /// class Platform {
  /// public:
  ///   typedef Node NodeType;
  ///
  ///   explicit Platform (const RCP<Node>& node);
  ///
  ///   RCP<const Comm<int> > getComm() const;
  ///
  ///   RCP<Node> getNode() const;
  /// };
  /// \endcode
  /// SerialPlatform uses a "communicator" containing one process.  It
  /// is available whether or not Trilinos was built with MPI.
  ///
  template <class Node>
  class SerialPlatform : public Teuchos::Describable {
  public:
    /// \typedef NodeType 
    /// \brief Kokkos Node type over which the platform is templated.
    typedef Node NodeType;

    //! @name Constructor/Destructor Methods
    //@{ 

    /// Constructor that accepts a Kokkos Node.
    ///
    /// \param node [in/out] The Kokkos Node instance.
    ///
    explicit SerialPlatform (const RCP<Node> &node) :
      comm_ (rcp (new Teuchos::SerialComm<int>())),
      node_ (node)
    {}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~SerialPlatform() {}

    //@}

    //! @name Class Creation and Accessor Methods
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
  /// \brief SerialPlatform specialization for \c Kokkos::DefaultNode::DefaultNodeType.
  ///
  /// \warning \c Kokkos::DefaultNode::DefaultNodeType is a typedef,
  ///   and may have a different type, depending on Trilinos' build
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
    /// \typedef NodeType 
    /// \brief Kokkos Node type over which the platform is templated.
    typedef Kokkos::DefaultNode::DefaultNodeType NodeType;

    //! @name Constructor/Destructor Methods
    //@{ 

    /// \brief Default constructor: uses Kokkos default node.
    ///
    /// The specialization of SerialPlatform for the default node type
    /// includes a default constructor.  It instantiates a default
    /// Node instance using Kokkos::DefaultNode::getDefaultNode().
    ///
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
    ///
    explicit SerialPlatform (const RCP<Kokkos::DefaultNode::DefaultNodeType> &node) :
      comm_ (rcp (new Teuchos::SerialComm<int> ())),
      node_ (node.is_null() ? Kokkos::DefaultNode::getDefaultNode() : node)
    {}

    //! Destructor (virtual for memory safety of derived classes).
    virtual ~SerialPlatform() {}

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! The \c Teuchos::Comm instance with which this object was created.
    RCP<const Teuchos::Comm<int> > getComm() const {
      return comm_;
    }

    /// \brief The Kokkos Node instance for this platform to use.
    /// 
    /// This SerialPlatform specialization's constructor may have
    /// created the Node instance, if you invoked the default
    /// constructor or passed in a null Node pointer.
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
