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

#ifndef TPETRA_HYBRIDPLATFORM_HPP
#define TPETRA_HYBRIDPLATFORM_HPP

#include <Teuchos_Describable.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Comm.hpp> 
#include <Teuchos_TypeNameTraits.hpp>
#include <string>

#include <Kokkos_SerialNode.hpp>
#ifdef HAVE_KOKKOS_TBB
#include <Kokkos_TBBNode.hpp>
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include <Kokkos_TPINode.hpp>
#endif
#ifdef HAVE_KOKKOS_THRUST
#include <Kokkos_ThrustGPUNode.hpp>
#endif

namespace Tpetra {

	//! \brief A platform class for hybrid nodes.
  /*!
    This class is templated on two types, those of the two underlying Nodes.
    In this way, the HybridPlatform is compiled with support for a particular 
    hybrid architecture.
   */
  class HybridPlatform : public Teuchos::Describable {
    public:
      //! @name Constructor/Destructor Methods
      //@{ 

      //! Constructor
      HybridPlatform(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, Teuchos::ParameterList &pl);

      //! Destructor
      ~HybridPlatform();

      //@}

      //! @name Class Query, Creation and Accessor Methods
      //@{ 

      //! Comm Instance
      Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

      void createNode();

      //! Run user code with the runtime-selected Node type.
      template <template <class Node> class UserCode> 
      void runUserCode();

      //@}

    private:
      HybridPlatform(const HybridPlatform &platform); // not supported
      const Teuchos::RCP<const Teuchos::Comm<int> > comm_;
      mutable Teuchos::ParameterList sublist;
      Teuchos::RCP<Kokkos::SerialNode>    serialNode_;
#ifdef HAVE_KOKKOS_TBB
      Teuchos::RCP<Kokkos::TBBNode>       tbbNode_;
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
      Teuchos::RCP<Kokkos::TPINode>       tpiNode_;
#endif
#ifdef HAVE_KOKKOS_THRUST
      Teuchos::RCP<Kokkos::ThrustGPUNode> thrustNode_;
#endif

      enum NodeType {
        SERIALNODE
#ifdef HAVE_KOKKOS_TBB
        , TBBNODE
#endif        
#ifdef HAVE_KOKKOS_THREADPOOL
        , TPINODE
#endif        
#ifdef HAVE_KOKKOS_THRUST
        , THRUSTNODE
#endif        
      } nodeType_;
  };

  HybridPlatform::HybridPlatform(const Teuchos::RCP<const Teuchos::Comm<int> > &comm, Teuchos::ParameterList &pl)
  : comm_(comm)
  , nodeType_(SERIALNODE)
  {
    // FINISH: process the parameter list
    switch (comm_->getRank()) {
    case 0:
      nodeType_ = SERIALNODE;
      break;
    case 1:
      nodeType_ = TPINODE;
      break;
    }
  } 

  HybridPlatform::~HybridPlatform() 
  {}

  Teuchos::RCP<const Teuchos::Comm<int> > 
  HybridPlatform::getComm() const {
    return comm_;
  }

  template<template<class Node> class UserCode>
  void HybridPlatform::runUserCode() {
    switch (nodeType_) {
      case SERIALNODE:
        UserCode<Kokkos::SerialNode>::run(comm_, serialNode_);
        break;
#ifdef HAVE_KOKKOS_TBB
      case TBBNODE:
        UserCode<Kokkos::TBBNode>::run(comm_, tbbNode_);
        break;
#endif        
#ifdef HAVE_KOKKOS_THREADPOOL
      case TPINODE:
        UserCode<Kokkos::TPINode>::run(comm_, tpiNode_);
        break;
#endif        
#ifdef HAVE_KOKKOS_THRUST
      case THRUSTNODE:
        UserCode<Kokkos::ThrustGPUNode>::run(comm_, thrustNode_);
        break;
#endif        
      default:
        TEST_FOR_EXCEPTION(true, std::runtime_error, 
            Teuchos::typeName(*this) << "::runUserCode(): Invalid node type." << std::endl);
    } // end of switch
  }

} // namespace Tpetra

#endif // TPETRA_HYBRIDPLATFORM_HPP
