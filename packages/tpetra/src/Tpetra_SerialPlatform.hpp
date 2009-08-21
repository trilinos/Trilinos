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
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Tpetra {

	//! \brief A implementation of the Platform class for serial platforms.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
     type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template<class Node=Kokkos::DefaultNode::DefaultNodeType>
	class SerialPlatform : public Teuchos::Describable {
  public:
    typedef Node NodeType;
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor
    SerialPlatform(Teuchos::RCP<Node> node);

    //! Destructor
    ~SerialPlatform();

    //@}

    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Comm Instance
    Teuchos::RCP< Teuchos::Comm<int> > getComm() const;

    //! Get Get a node for parallel computation.
    Teuchos::RCP<Node> getNode() const;

    //@}
    private:
    SerialPlatform(const SerialPlatform<Node> &platform);

    protected: 
    Teuchos::RCP<Node> node_;
  };

  template<class Node>
  SerialPlatform<Node>::SerialPlatform(Teuchos::RCP<Node> node) 
  : node_(node) {}

  template<class Node>
  SerialPlatform<Node>::~SerialPlatform() {}

  template<class Node>
  Teuchos::RCP< Teuchos::Comm<int> > 
  SerialPlatform<Node>::getComm() const 
  {
    return Teuchos::rcp(new Teuchos::SerialComm<int>() );
  }

  template<class Node>
  Teuchos::RCP<Node> SerialPlatform<Node>::getNode() const 
  { return node_; }

} // namespace Tpetra

#endif // TPETRA_SERIALPLATFORM_HPP
