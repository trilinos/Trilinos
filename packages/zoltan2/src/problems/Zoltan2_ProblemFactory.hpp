// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

/*! \file Zoltan2_ProblemFactory.hpp
    \brief Returns a pointer to a new problem class.
    Is not responsible for memory management!
*/
#ifndef ZOLTAN2_PROBLEM_FACTORY_HPP
#define ZOLTAN2_PROBLEM_FACTORY_HPP
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>
#include <Zoltan2_OrderingProblem.hpp>

namespace Zoltan2 {
/// \brief ProblemFactory class contains 1 static factory method
///        templated the input adapter type
  class ProblemFactory {
  public:

    /// \brif Zoltan2::Problem factory method
    ///
    /// @param kind A string equal to the type of the problem (paritioning, ordering, coloring)
    /// @param input Input adapter used to construct the problem
    /// @param params Zolta2 parameter list
    /// @param (MPI) MPI world communicator
    ///
    /// @return returns a pointer to new Zoltan2::Problem or a nullptr if kind was not known.
    template<typename base_adapter, typename cast_adapter>
#ifdef HAVE_ZOLTAN2_MPI
    static Problem<base_adapter> * newProblem( const std::string &kind,
                                          base_adapter *input,
                                          ParameterList *params,
                                          MPI_Comm comm)
#else
    static Problem<base_adapter> * newProblem( const std::string & kind,
                                          base_adapter *input,
                                          ParameterList *params)
#endif
    {
      if(kind == "partitioning") {
#ifdef HAVE_ZOLTAN2_MPI
        return reinterpret_cast< Problem<base_adapter> *>( new Zoltan2::PartitioningProblem<cast_adapter>(reinterpret_cast<cast_adapter *>(input), params, comm));
#else
        return reinterpret_cast< Problem<base_adapter> *>( new Zoltan2::PartitioningProblem<cast_adapter>(reinterpret_cast<cast_adapter *>(input), params));
#endif
      } else if (kind == "ordering") {

#ifdef HAVE_ZOLTAN2_MPI
        return reinterpret_cast< Problem<base_adapter> *>( new Zoltan2::OrderingProblem<cast_adapter>(reinterpret_cast<cast_adapter *>(input), params, comm));
#else
        return reinterpret_cast< Problem<base_adapter> *>( new Zoltan2::OrderingProblem<cast_adapter>(reinterpret_cast<cast_adapter *>(input), params));
#endif
      } else if (kind == "coloring") {

#ifdef HAVE_ZOLTAN2_MPI
        return reinterpret_cast< Problem<base_adapter> *>( new Zoltan2::ColoringProblem<cast_adapter>(reinterpret_cast<cast_adapter *>(input), params, comm));
#else
        return reinterpret_cast< Problem<base_adapter> *>( new Zoltan2::ColoringProblem<cast_adapter>(reinterpret_cast<cast_adapter *>(input), params));
#endif
      }
     return nullptr; // problem type not known 
    }
  };
}
#endif // ZOLTAN2_PROBLEM_FACTORY_HPP

