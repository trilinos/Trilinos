// @HEADER
// ***********************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
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

#ifndef KOKKOS_DEFAULT_NODE_HPP_
#define KOKKOS_DEFAULT_NODE_HPP_

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_SerialNode.hpp"
#ifdef HAVE_KOKKOS_TBB
#include "Kokkos_TBBNode.hpp"
#endif
#ifdef HAVE_KOKKOS_THREADPOOL
#include "Kokkos_TPINode.hpp"
#endif

#include <Teuchos_RCP.hpp>

namespace Kokkos {

  /** \brief Class to specify %Kokkos default node type and instantiate the default node.
      \ingroup kokkos_node_api
    */
  class DefaultNode {
    public:
#if defined(HAVE_KOKKOS_THREADPOOL)
      typedef TPINode DefaultNodeType;
#elif defined(HAVE_KOKKOS_TBB)
      typedef TBBNode DefaultNodeType;
#else
      //! Typedef specifying the default node type.
      typedef SerialNode DefaultNodeType;
#endif

      //! \brief Return a pointer to the default node.
      static RCP<DefaultNodeType> getDefaultNode();

    private:
      static RCP<DefaultNodeType> node_;
  };

}

#endif
