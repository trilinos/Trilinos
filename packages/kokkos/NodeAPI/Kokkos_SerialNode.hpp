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

#ifndef KOKKOS_SERIALNODE_HPP_
#define KOKKOS_SERIALNODE_HPP_

#include <Teuchos_ParameterList.hpp>
#include <Kokkos_StandardNodeMemoryModel.hpp>
#include "Kokkos_NodeHelpers.hpp"

namespace Kokkos {

  /** \brief %Kokkos node interface for a serial, CPU node.
      \ingroup kokkos_node_api
   */
  class SerialNode : public StandardNodeMemoryModel {
    public:
      /*! \brief Default constructor, accepts a parameter list but reads no parameters. */
      SerialNode(Teuchos::ParameterList &pl) {}

      //! \begin parallel for skeleton, with a trivial serial implementation. See \ref kokkos_node_api "Kokkos Node API"
      template <class WDP>
      static void parallel_for(int beg, int end, WDP wd) {
        for (int i=beg; i != end; ++i) {
          wd.execute(i);
        }
      }

      //! \begin parallel reduction skeleton, with a trivial serial implementation. See \ref kokkos_node_api "Kokkos Node API"
      template <class WDP>
      static typename WDP::ReductionType
      parallel_reduce(int begin, int end, WDP wd) {
        typename WDP::ReductionType result = wd.identity();
        for (int i=begin; i != end; ++i) {
          result = wd.reduce( result, wd.generate(i) );
        }
        return result;
      }

      //! \begin No-op for SerialNode.
      inline void sync() const {};

  };

  template <> class ArrayOfViewsHelper<SerialNode> : public ArrayOfViewsHelperTrivialImpl<SerialNode> {};

} // end of Kokkos namespace

#endif
