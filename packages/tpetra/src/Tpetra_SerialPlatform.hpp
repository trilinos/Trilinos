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

	//! \brief A implementation of the Platform class for serial platforms.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template <class Node>
	class SerialPlatform : public Teuchos::Describable {
  public:
    //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
    typedef Node NodeType;
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Node-accepting constructor
    explicit SerialPlatform(const RCP<Node> &node) {
      node_ = node;
      comm_ = rcp(new Teuchos::SerialComm<int>() );
    }

    //! Destructor
    ~SerialPlatform() {}

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Comm Instance
    RCP< const Comm<int> > getComm() const {
      return comm_;
    }

    //! Get Get a node for parallel computation.
    RCP<Node> getNode() const {
      return node_; 
    }

    //@}
    private:
    SerialPlatform(const SerialPlatform<Node> &platform);

    protected: 
    //! Teuchos::Comm object instantiated for the platform.
    RCP<const Teuchos::SerialComm<int> > comm_;
    //! Node object instantiated for the platform.
    RCP<Node> node_;
  };


  template <>
	class SerialPlatform<Kokkos::DefaultNode::DefaultNodeType> : public Teuchos::Describable {
  private:
    SerialPlatform(const SerialPlatform<Kokkos::DefaultNode::DefaultNodeType> &platform);

  protected: 
    //! Teuchos::Comm object instantiated for the platform.
    RCP<const Teuchos::SerialComm<int> > comm_;
    //! Node object instantiated for the platform.
    RCP<Kokkos::DefaultNode::DefaultNodeType> dnode_;

  public:
    //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
    typedef Kokkos::DefaultNode::DefaultNodeType NodeType;
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Default constructor uses Kokkos default node
    SerialPlatform() {
      // will construct the node late
      dnode_ = null;
      comm_ = rcp(new Teuchos::SerialComm<int>() );
    }

    //! Node-accepting constructor
    explicit SerialPlatform(const RCP<Kokkos::DefaultNode::DefaultNodeType> &node) {
      dnode_ = node;
      comm_ = rcp(new Teuchos::SerialComm<int>() );
    }

    //! Destructor
    ~SerialPlatform() {}

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Comm Instance
    RCP< const Teuchos::SerialComm<int> > getComm() const {
      return comm_;
    }

    //! Get Get a node for parallel computation.
    RCP<Kokkos::DefaultNode::DefaultNodeType> getNode() const {
      RCP<Kokkos::DefaultNode::DefaultNodeType> def = dnode_;
      if (def == null) {
        def = Kokkos::DefaultNode::getDefaultNode(); 
      }
      return def;
    }

    //@}

  };

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
