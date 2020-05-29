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

/*! \file Zoltan2_TestFactory.hpp
    \brief Returns a pointer to new test classes.
    Is not responsible for memory management!
*/
#ifndef ZOLTAN2_PROBLEM_FACTORY_HPP
#define ZOLTAN2_PROBLEM_FACTORY_HPP
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_OrderingProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>
#include <Zoltan2_Typedefs.hpp>

using namespace Zoltan2_TestingFramework;
using namespace Zoltan2;

namespace Zoltan2_TestingFramework {
/// \brief ProblemFactory class contains 1 static factory method
  class ProblemFactory {
  public:
    /// \brif Zoltan2::Problem factory method
    ///
    /// @param kind A string equal to the type of the problem (paritioning, ordering, coloring)
    /// @param kind A string equal to the type of the adapter
    /// @param params Zolta2 parameter list
    /// @param (MPI) MPI world communicator
    ///
    /// @return returns a pointer to new Zoltan2::Problem or a nullptr if
    /// problem_name was not known.
    ProblemFactory(const std::string & problemName,
                   RCP<AdapterFactory> adapterFactory,
                   ParameterList *params
                   #ifdef HAVE_ZOLTAN2_MPI
                    , MPI_Comm comm
                   #endif
                   ) {

      problem_name = problemName;
      adapterType = adapterFactory->getMainAdapterType();

      #ifdef HAVE_ZOLTAN2_MPI
        #define CREATE_PRBLM(problemClass, adapterClass)                       \
          adapterClass * pCast = dynamic_cast<adapterClass *>                  \
            (adapterFactory->getMainAdapter());                                \
          if(!pCast) { throw std::logic_error(                                 \
            "ProblemFactory adapter dynamic_cast failed for problem name "     \
              + problem_name + " and adapterClass " + #adapterClass ); }       \
          problem = rcp(new problemClass<adapterClass>(pCast, params, comm));
      #else
        #define CREATE_PRBLM(problemClass, adapterClass)                       \
          adapterClass * pCast = dynamic_cast<adapterClass *>                  \
            (adapterFactory->getMainAdapter());                                \
          if(!pCast) { throw std::logic_error(                                 \
            "ProblemFactory adapter dynamic_cast failed for problem name "     \
              + problem_name + " and adapterClass " + #adapterClass ); }       \
          problem = rcp(new problemClass<adapterClass>(pCast, params));
      #endif

      #define MAKE_PARTITION_PROBLEM(adapterClass)  \
        CREATE_PRBLM(PartitioningProblem, adapterClass);

      #define MAKE_ORDERING_PROBLEM(adapterClass)  \
        CREATE_PRBLM(OrderingProblem, adapterClass);

      // PartitioningProblem
      if(problem_name == "partitioning") {
        Z2_TEST_UPCAST(adapterType, MAKE_PARTITION_PROBLEM)
      }
      else if(problem_name == "ordering") {
        Z2_TEST_UPCAST(adapterType, MAKE_ORDERING_PROBLEM)
      }

      if(problem == Teuchos::null) {
        throw std::logic_error("ProblemFactory failed to create Problem!");
      }
    }

    RCP<ProblemRoot> getProblem() { return problem; }
    const std::string & getProblemName() const { return problem_name; }
    EAdapterType getAdapterType() const { return adapterType; }

    private:
      std::string problem_name; // string converts to a problem type
      EAdapterType adapterType; // converts to an adapter type
      RCP<ProblemRoot> problem;
  };
}
#endif // ZOLTAN2_PROBLEM_FACTORY_HPP

