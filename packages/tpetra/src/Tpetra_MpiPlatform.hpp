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
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Tpetra {

	//! \brief A implementation of the Platform class for MPI-based platforms.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template<class Node=Kokkos::DefaultNode::DefaultNodeType>
	class MpiPlatform : public Teuchos::Describable {
    public:
      typedef Node NodeType;
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor
      explicit MpiPlatform(Teuchos::RCP<Node> node);

      //! Constructor
      MpiPlatform(Teuchos::RCP<Node> node, const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm);

      //! Destructor
      ~MpiPlatform();

      //@}

      //! @name Class Creation and Accessor Methods
      //@{ 

      //! Comm Instance
      Teuchos::RCP< const Teuchos::Comm<int> > getComm() const;

      //! Get Get a node for parallel computation.
      Teuchos::RCP<Node> getNode() const;

      //@}

    protected: 
      Teuchos::RCP<Node> node_;

    private:
      Teuchos::RCP<Teuchos::MpiComm<int> > comm_;
      MpiPlatform(const MpiPlatform<Node> &platform);
  };

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(Teuchos::RCP<Node> node, const Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > &rawMpiComm)
  : node_(node)
  {
    comm_ = Teuchos::createMpiComm<int>(rawMpiComm);
  }

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(Teuchos::RCP<Node> node)
  : node_(node)
  {
    comm_ = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  } 

  template <class Node>
  MpiPlatform<Node>::~MpiPlatform() {}

  template <class Node>
  MpiPlatform<Node>::MpiPlatform(const MpiPlatform<Node> &platform)
  {
    comm_ = platform.comm_;
  }

  template <class Node>
  Teuchos::RCP< const Teuchos::Comm<int> > 
  MpiPlatform<Node>::getComm() const 
  {
    return comm_;
  }

  template <class Node>
  Teuchos::RCP<Node> MpiPlatform<Node>::getNode() const 
  { return node_; }

} // namespace Tpetra

#endif // TPETRA_MPIPLATFORM_HPP

