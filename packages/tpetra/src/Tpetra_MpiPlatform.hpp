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

#ifndef TPETRA_MPIPLATFORM_HPP
#define TPETRA_MPIPLATFORM_HPP

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Tpetra {

	//! \brief A implementation of the Platform class for MPI-based platforms.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template <class Node>
	class MpiPlatform : public Teuchos::Describable {
  private:
    MpiPlatform(const MpiPlatform<Node> &platform);

  protected: 
    //! Teuchos::Comm object instantiated for the platform.
    RCP<Teuchos::MpiComm<int> > comm_;
    //! Node object instantiated for the platform.
    RCP<Node> node_;
  public:
    //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
    typedef Node NodeType;
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Node-accepting constructor uses MPI_COMM_WORLD
    explicit MpiPlatform(const RCP<Node> &node) {
      node_ = node;
      comm_ = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
    }

    //! Node and MPI_Comm accepting constructor
    MpiPlatform(const RCP<Node> &node, const RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm) {
      node_ = node;
      comm_ = Teuchos::createMpiComm<int>(rawMpiComm);
    }

    //! Destructor
    ~MpiPlatform() {}

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
  };


  template <>
	class MpiPlatform<Kokkos::DefaultNode::DefaultNodeType> : public Teuchos::Describable {
  private:
    MpiPlatform(const MpiPlatform<Kokkos::DefaultNode::DefaultNodeType> &platform);

  protected: 
    //! Teuchos::Comm object instantiated for the platform.
    RCP<Teuchos::MpiComm<int> > comm_;
    //! Node object instantiated for the platform.
    RCP<Kokkos::DefaultNode::DefaultNodeType> dnode_;
  public:
    //! Typedef indicating the node type over which the platform is templated. This default to the Kokkos default node type.
    typedef Kokkos::DefaultNode::DefaultNodeType NodeType;
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Default constructor uses Kokkos default node and MPI_COMM_WORLD
    MpiPlatform() {
      dnode_ = null;
      comm_ = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
    }

    //! Node and MPI_Comm accepting constructor
    MpiPlatform(const RCP<Kokkos::DefaultNode::DefaultNodeType> &node, const RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm) {
      dnode_ = node;
      comm_ = Teuchos::createMpiComm<int>(rawMpiComm);
    }

    //! Destructor
    ~MpiPlatform() {}

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Comm Instance
    RCP< const Comm<int> > getComm() const {
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

#endif // TPETRA_MPIPLATFORM_HPP
