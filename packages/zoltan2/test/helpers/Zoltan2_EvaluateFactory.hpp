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

/*! \file Zoltan2_EvaluateFactory.hpp
    \brief Returns a pointer to new test classes.
    Is not responsible for memory management!
*/
#ifndef ZOLTAN2_EVALUATE_FACTORY_HPP
#define ZOLTAN2_EVALUATE_FACTORY_HPP
#include <Zoltan2_Typedefs.hpp>
#include <Zoltan2_EvaluatePartitionFactory.hpp>
#include <Zoltan2_EvaluateOrderingFactory.hpp>
// #include <Zoltan2_EvaluateColoringFactory.hpp>

#include <assert.h>
using namespace Zoltan2_TestingFramework;
using namespace Zoltan2;

namespace Zoltan2_TestingFramework {
/// \brief ProblemFactory class contains 1 static factory method
  class EvaluateFactory {
  public:

    /// \brif Zoltan2::EvaluateBaseClass factory method
    ///
    /// @param kind A string equal to the type of the problem (paritioning, ordering, coloring)
    /// @param kind A string equal to the type of the adapter
    /// @param input Input adapter used to construct the problem
    /// @param params Zolta2 parameter list
    /// @param (MPI) MPI world communicator
    ///
    /// @return returns a pointer to new Zoltan2::Problem or a nullptr if kind was not known.

    static base_evaluate_t * newEvaluate(const std::string & kind,
                                          const std::string &adapter_name,
                                          base_adapter_t *input,
                                          ParameterList *params,
                                          base_problem_t * problem) {
      if (kind == "partitioning") {
        return reinterpret_cast<base_evaluate_t *>(
          EvaluatePartitionFactory::newEvaluatePartition(
                   reinterpret_cast<partitioning_problem_t*> (problem),
                   adapter_name, input, params));
      }
      else if(kind == "ordering") {
        return reinterpret_cast<base_evaluate_t *>(
          EvaluateOrderingFactory::newEvaluateOrdering(
                   reinterpret_cast<ordering_problem_t*> (problem),
                   adapter_name, input, params));
      }
      else if(kind == "coloring") {
        // Coloring code here... EvaluateColoringFactory not created yet
        // return reinterpret_cast<base_evaluate_t *>(
        //  EvaluateColoringFactory::newEvaluatColoring(
        //           reinterpret_cast<coloring_problem_t*> (problem),
        //           adapter_name, input, params));
      }

      return nullptr;
    }
  };
}
#endif // ZOLTAN2_EVALUATE_FACTORY_HPP

