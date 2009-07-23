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

#ifndef TPETRA_PLATFORM_HPP
#define TPETRA_PLATFORM_HPP

#include "Tpetra_ConfigDefs.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_Describable.hpp>
#include <Kokkos_DefaultNode.hpp>

namespace Tpetra {

  //! \brief The Tpetra platform abstract base class.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
     type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class Platform : public Teuchos::Describable {
  public:
    typedef Scalar        ScalarType;
    typedef LocalOrdinal  LocalOrdinalType;
    typedef GlobalOrdinal GlobalOrdinalType;
    typedef Node          NodeType;
  
    //! @name Constructor/Destructor Methods
    //@{ 

    //! Constructor
    Platform(Node &node);

    //! Constructor with object label
    Platform(const std::string &label);

    //! Destructor
    virtual ~Platform() {};

    //! Clone method
    /*! Returns a copy of this Platform instance. It is dynamically allocated and 
        encapsulated in a Teuchos RCP.
    */
    virtual Teuchos::RCP< Platform<Scalar, LocalOrdinal, GlobalOrdinal, Node> > clone() const = 0;

    //@}
  
    //! @name Class Creation and Accessor Methods
    //@{ 

    //! Create a Comm instance for global communication between nodes.
    virtual Teuchos::RCP< Teuchos::Comm<int> > getComm() const = 0;

    //! Get Get a node for parallel computation.
    virtual Node & getNode() const;

    //@}

    protected: 
    Node &node_;
  }; // Platform class

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Platform<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Platform(Node &node) 
  : node_(node) {}
  
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Platform<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Platform(const std::string &str)
  : Teuchos::LabeledObject(str) {}

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Node &
  Platform<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getNode() const 
  { return node_; }

  
} // namespace Tpetra

#endif // TPETRA_PLATFORM_HPP

