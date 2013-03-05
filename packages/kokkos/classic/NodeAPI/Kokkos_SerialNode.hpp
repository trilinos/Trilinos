//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER

#ifndef KOKKOS_SERIALNODE_HPP_
#define KOKKOS_SERIALNODE_HPP_

#include <Kokkos_StandardNodeMemoryModel.hpp>
#include "Kokkos_NodeHelpers.hpp"

namespace Kokkos {

  /** \brief %Kokkos node interface for a serial, CPU node.
      \ingroup kokkos_node_api
   */
  class SerialNode : public StandardNodeMemoryModel {
  public:
    //! Constructor; sets default parameters.
    SerialNode () {}

    /// Constructor that takes a list of parameters.
    /// This class doesn't currently use any parameters.
    SerialNode (ParameterList &pl) {
      (void) pl; 
    }

    //! Get default parameters for this class.
    static ParameterList getDefaultParameters() {
      ParameterList params;
      return params;
    }

    /// Skeleton for a parallel for, with a trivial serial implementation. 
    /// See \ref kokkos_node_api "Kokkos Node API"
    template <class WDP>
    static void parallel_for(int beg, int end, WDP wd) {
      for (int i=beg; i != end; ++i) {
	wd.execute(i);
      }
    }

    /// Skeleton for a parallel reduction, with a trivial serial implementation. 
    /// See \ref kokkos_node_api "Kokkos Node API"
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

  template <> class ArrayOfViewsHelper<SerialNode> : 
    public ArrayOfViewsHelperTrivialImpl<SerialNode> {};

} // end of Kokkos namespace

#endif
